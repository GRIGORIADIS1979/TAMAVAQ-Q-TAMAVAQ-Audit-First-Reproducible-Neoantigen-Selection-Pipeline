"""Reproducible TAMAVAQ / Q-TAMAVAQ pipeline.

This module is intentionally dependency-light and deterministic. It implements
classical simulations of the graph/CTQW and Grover-style weighting policies
specified in the manuscript. It is not a clinical decision tool.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from dataclasses import dataclass
from datetime import datetime, timezone
from itertools import permutations
from pathlib import Path
from typing import Any, Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd
import yaml
from scipy.linalg import expm

AA = set("ACDEFGHIKLMNPQRSTVWY")
PROPERTY_GROUPS = {
    "hydrophobic": set("AILMFWYV"),
    "polar": set("STNQCY"),
    "positive": set("KRH"),
    "negative": set("DE"),
    "special": set("GP"),
}


@dataclass(frozen=True)
class RunArtifacts:
    selected_panel: pd.DataFrame
    scores: pd.DataFrame
    batch_state: Dict[str, Any]
    construct_sequence: str
    teleport_events: List[Dict[str, Any]]


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def canonical_id(row: pd.Series) -> str:
    payload = "|".join(str(row.get(c, "")).strip() for c in ["peptide", "allele", "gene", "source_variant", "anchor_type"])
    return hashlib.sha256(payload.encode()).hexdigest()[:16]


def parse_bool(x: Any) -> bool:
    if isinstance(x, bool):
        return x
    return str(x).strip().lower() in {"true", "1", "yes", "y"}


def validate_candidates(df: pd.DataFrame) -> pd.DataFrame:
    required = [
        "peptide", "allele", "gene", "source_variant", "calibrated_immunogenicity",
        "delta_g_kcal_mol", "manufacturable", "provenance_ok", "expression_tpm",
        "vaf", "clonality", "anchor_type",
    ]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
    out = df.copy()
    out["peptide"] = out["peptide"].astype(str).str.upper().str.strip()
    out["manufacturable"] = out["manufacturable"].map(parse_bool)
    out["provenance_ok"] = out["provenance_ok"].map(parse_bool)
    numeric = ["calibrated_immunogenicity", "delta_g_kcal_mol", "expression_tpm", "vaf"]
    if "tumor_cell_fraction" in out.columns:
        numeric.append("tumor_cell_fraction")
    for col in numeric:
        out[col] = pd.to_numeric(out[col], errors="raise")
    out["valid_peptide_alphabet"] = out["peptide"].map(lambda p: bool(p) and set(p).issubset(AA))
    out["canonical_id"] = out.apply(canonical_id, axis=1)
    out = out.sort_values(["canonical_id", "peptide", "allele"], kind="mergesort").reset_index(drop=True)
    out["candidate_index"] = np.arange(len(out), dtype=int)
    return out


def property_fingerprint(peptide: str) -> np.ndarray:
    counts = []
    length = max(len(peptide), 1)
    for residues in PROPERTY_GROUPS.values():
        counts.append(sum(1 for aa in peptide if aa in residues) / length)
    # crude but deterministic length and terminal descriptors
    counts.append(length / 30.0)
    counts.append(1.0 if peptide and peptide[-1] in PROPERTY_GROUPS["hydrophobic"] else 0.0)
    counts.append(1.0 if peptide and peptide[0] in PROPERTY_GROUPS["positive"] else 0.0)
    return np.asarray(counts, dtype=float)


def build_feature_matrix(df: pd.DataFrame) -> np.ndarray:
    prop = np.vstack([property_fingerprint(p) for p in df["peptide"]])
    scalar = df[["calibrated_immunogenicity", "expression_tpm", "vaf"]].to_numpy(float)
    scalar[:, 1] = np.log1p(scalar[:, 1]) / (1.0 + np.log1p(max(df["expression_tpm"].max(), 1.0)))
    dg = df[["delta_g_kcal_mol"]].to_numpy(float)
    dg = (dg - dg.mean()) / (dg.std() + 1e-9)
    X = np.hstack([prop, scalar, dg])
    norms = np.linalg.norm(X, axis=1, keepdims=True)
    return X / np.maximum(norms, 1e-12)


def fs_graph(X: np.ndarray, sigma: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    K = np.clip((X @ X.T) ** 2, 0.0, 1.0)
    Dfs = np.arccos(np.sqrt(K))
    W = np.exp(-(Dfs ** 2) / (2.0 * sigma ** 2))
    np.fill_diagonal(W, 0.0)
    Deg = np.diag(W.sum(axis=1))
    L = Deg - W
    return K, Dfs, W, L


def initial_prior(df: pd.DataFrame, mode: str) -> np.ndarray:
    n = len(df)
    if mode == "uniform":
        return np.ones(n) / n
    bio = df["calibrated_immunogenicity"].to_numpy(float)
    phys = -df["delta_g_kcal_mol"].to_numpy(float)
    phys = phys - phys.min() + 1e-6
    expr = np.log1p(df["expression_tpm"].to_numpy(float)) + 1e-6
    prior = np.maximum(bio, 1e-6) * phys * expr
    return prior / prior.sum()


def ctqw_scores(L: np.ndarray, prior: np.ndarray, t: float) -> Tuple[np.ndarray, Dict[str, Any]]:
    psi0 = np.sqrt(prior).astype(complex)
    U = expm(-1j * L * t)
    psit = U @ psi0
    probs = np.abs(psit) ** 2
    probs = probs / probs.sum()
    eig = np.linalg.eigvalsh(L)
    entropy = float(-(probs * np.log2(np.maximum(probs, 1e-12))).sum())
    return probs.real, {"laplacian_eigenvalues": eig.round(10).tolist(), "walk_entropy_bits": entropy}


def predicate_table(df: pd.DataFrame, ctqw: np.ndarray, cfg: Dict[str, Any]) -> pd.DataFrame:
    th = cfg["thresholds"]
    pred = pd.DataFrame(index=df.index)
    pred["BIO"] = df["calibrated_immunogenicity"] >= float(th["bio_min"])
    pred["PHYS"] = df["delta_g_kcal_mol"] <= float(th["delta_g_max_kcal_mol"])
    pred["MFG"] = df["manufacturable"].astype(bool)
    pred["PROV"] = df["provenance_ok"].astype(bool)
    pred["VALID"] = df["valid_peptide_alphabet"].astype(bool)
    pred["GEOM"] = ctqw >= float(th.get("ctqw_min", 0.0))
    pred["MARK"] = pred[["BIO", "PHYS", "MFG", "PROV", "VALID", "GEOM"]].all(axis=1)
    return pred


def grover_amplified_scores(ctqw: np.ndarray, mark: np.ndarray, iterations: int) -> np.ndarray:
    n = len(ctqw)
    m = int(mark.sum())
    if m == 0 or m == n:
        return ctqw / ctqw.sum()
    theta = math.asin(math.sqrt(m / n))
    marked_mass = math.sin((2 * iterations + 1) * theta) ** 2
    out = np.zeros(n, dtype=float)
    marked_base = ctqw * mark
    unmarked_base = ctqw * (~mark)
    out[mark] = marked_mass * marked_base[mark] / marked_base.sum()
    out[~mark] = (1.0 - marked_mass) * unmarked_base[~mark] / unmarked_base.sum()
    return out / out.sum()


def quota_ok(current: pd.DataFrame, cand: pd.Series, cfg: Dict[str, Any]) -> bool:
    panel_cfg = cfg["panel"]
    if len(current) == 0:
        allele_count = 0
    else:
        allele_count = int((current["allele"] == cand["allele"]).sum())
    if allele_count >= int(panel_cfg.get("max_per_allele", 10**9)):
        return False
    return True


def select_panel(df: pd.DataFrame, scores: np.ndarray, pred: pd.DataFrame, cfg: Dict[str, Any]) -> pd.DataFrame:
    work = df.copy()
    work["score"] = scores
    work["MARK"] = pred["MARK"].values
    candidates = work[work["MARK"]].copy()
    candidates = candidates.sort_values(
        ["score", "expression_tpm", "vaf", "canonical_id"],
        ascending=[False, False, False, True],
        kind="mergesort",
    )
    selected: List[pd.Series] = []
    k = int(cfg["panel"]["k"])
    seen_peptides = set()
    for _, row in candidates.iterrows():
        if row["peptide"] in seen_peptides:
            continue
        current = pd.DataFrame(selected) if selected else pd.DataFrame(columns=candidates.columns)
        if not quota_ok(current, row, cfg):
            continue
        selected.append(row)
        seen_peptides.add(row["peptide"])
        if len(selected) >= k:
            break
    panel = pd.DataFrame(selected) if selected else candidates.head(0)
    min_long = int(cfg["panel"].get("min_long_anchors", 0))
    if min_long and not panel.empty and (panel["anchor_type"].str.lower().str.contains("long").sum() < min_long):
        long_candidates = candidates[candidates["anchor_type"].str.lower().str.contains("long")]
        if not long_candidates.empty:
            replacement = long_candidates.iloc[0]
            panel = pd.concat([pd.DataFrame([replacement]), panel[panel["canonical_id"] != replacement["canonical_id"]]], ignore_index=True).head(k)
    return panel.reset_index(drop=True)


def junction_risk(left: str, right: str, window: int) -> float:
    junction = (left + right)
    center = len(left)
    start = max(0, center - window + 1)
    end = min(len(junction), center + window - 1)
    segment = junction[start:end]
    # deterministic proxy: penalize hydrophobic runs and repeated amino acids at junction
    hydrophobic = sum(1 for aa in segment if aa in PROPERTY_GROUPS["hydrophobic"]) / max(len(segment), 1)
    repeats = sum(1 for a, b in zip(segment, segment[1:]) if a == b) / max(len(segment) - 1, 1)
    return min(1.0, 0.7 * hydrophobic + 0.3 * repeats)


def assemble_construct(panel: pd.DataFrame, cfg: Dict[str, Any]) -> Tuple[str, Dict[str, Any]]:
    peptides = panel["peptide"].tolist()
    if not peptides:
        return "", {"order": [], "linkers": []}
    window = int(cfg["construct"].get("junction_window", 9))
    risk_max = float(cfg["construct"].get("direct_junction_risk_max", 0.35))
    linkers = list(cfg["construct"].get("linker_library", ["AAY", "GPGPG", "GGGS"]))

    def order_cost(order: Tuple[str, ...]) -> float:
        return sum(junction_risk(a, b, window) for a, b in zip(order, order[1:]))

    if len(peptides) <= 8:
        order = min(permutations(peptides), key=lambda p: (order_cost(p), p))
    else:
        remaining = sorted(peptides)
        order_list = [remaining.pop(0)]
        while remaining:
            nxt = min(remaining, key=lambda p: (junction_risk(order_list[-1], p, window), p))
            remaining.remove(nxt)
            order_list.append(nxt)
        order = tuple(order_list)

    seq = order[0]
    linker_manifest = []
    for left, right in zip(order, order[1:]):
        direct = junction_risk(left, right, window)
        chosen = ""
        if direct > risk_max:
            chosen = min(linkers, key=lambda l: (junction_risk(left, l + right, window), len(l), l))
        seq += chosen + right
        linker_manifest.append({"left": left, "right": right, "direct_risk": direct, "linker": chosen})
    return seq, {"order": list(order), "linkers": linker_manifest, "lambda_evaluation": "run neoDesign/LinearDesign lambda evaluation downstream"}


def make_teleport_event(batch_id: str, context: str, predicate_snapshot: Dict[str, Any], prior_checksum: str) -> Dict[str, Any]:
    payload = json.dumps({"batch_id": batch_id, "context": context, "predicate_snapshot": predicate_snapshot, "prior": prior_checksum}, sort_keys=True)
    digest = sha256_bytes(payload.encode())
    bits = [int(x, 16) % 2 for x in digest[:4]]
    frame = "".join("X" if b else "I" for b in bits[:2]) + ":" + "".join("Z" if b else "I" for b in bits[2:])
    return {
        "event_id": digest[:16],
        "batch_id": batch_id,
        "context": context,
        "bell_outcomes": bits,
        "pauli_frame": frame,
        "predicate_snapshot": predicate_snapshot,
        "checksum": digest,
    }


def run_pipeline(candidates_path: Path, config_path: Path, outdir: Path) -> RunArtifacts:
    cfg = yaml.safe_load(config_path.read_text())
    batch_id = str(cfg.get("batch_id", "tamavaq_batch"))
    df = validate_candidates(pd.read_csv(candidates_path))
    index_map = df[["candidate_index", "canonical_id", "peptide", "allele"]].to_dict(orient="records")
    index_hash = sha256_bytes(json.dumps(index_map, sort_keys=True).encode())

    X = build_feature_matrix(df)
    K, Dfs, W, L = fs_graph(X, sigma=float(cfg["ctqw"].get("sigma", 0.85)))
    prior = initial_prior(df, cfg["ctqw"].get("prior", "biology_physics"))
    ctqw, diag = ctqw_scores(L, prior, t=float(cfg["ctqw"].get("time", 1.0)))
    pred = predicate_table(df, ctqw, cfg)
    amplified = grover_amplified_scores(ctqw, pred["MARK"].to_numpy(bool), int(cfg["amplification"].get("grover_iterations", 1)))
    panel = select_panel(df, amplified, pred, cfg)
    construct, construct_info = assemble_construct(panel, cfg)

    scores = df[["candidate_index", "canonical_id", "peptide", "allele", "gene", "anchor_type"]].copy()
    scores["prior"] = prior
    scores["ctqw_score"] = ctqw
    scores["amplified_score"] = amplified
    for col in pred.columns:
        scores[col] = pred[col].values

    lap_hash = sha256_bytes(np.ascontiguousarray(L).tobytes())
    selected_ids = panel.get("canonical_id", pd.Series(dtype=str)).tolist()
    batch_state = {
        "batch_id": batch_id,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "input_sha256": sha256_file(candidates_path),
        "config_sha256": sha256_file(config_path),
        "candidate_index_map_sha256": index_hash,
        "software": {"tamavaq_algorithm_package": "0.1.0", "python_module": "tamavaq.pipeline"},
        "thresholds": cfg["thresholds"],
        "ctqw": {**cfg["ctqw"], **diag, "laplacian_sha256": lap_hash},
        "amplification": cfg["amplification"],
        "predicate_counts": {c: int(pred[c].sum()) for c in pred.columns},
        "selected_ids": selected_ids,
        "construct": construct_info,
        "scope": "Computational selector/audit package; not a clinical efficacy predictor.",
    }
    state_checksum = sha256_bytes(json.dumps(batch_state, sort_keys=True).encode())
    pred_snapshot = {c: int(pred[c].sum()) for c in pred.columns}
    events = [
        make_teleport_event(batch_id, "post_ctqw_to_predicate", pred_snapshot, lap_hash),
        make_teleport_event(batch_id, "post_predicate_to_panel", {"selected_n": len(selected_ids), **pred_snapshot}, state_checksum),
    ]

    outdir.mkdir(parents=True, exist_ok=True)
    scores.to_csv(outdir / "scores.csv", index=False)
    panel.to_csv(outdir / "selected_panel.csv", index=False)
    (outdir / "optimized_construct.txt").write_text(construct + "\n")
    (outdir / "batch_state.json").write_text(json.dumps(batch_state, indent=2, sort_keys=True) + "\n")
    with (outdir / "teleport_log.jsonl").open("w") as fh:
        for ev in events:
            fh.write(json.dumps(ev, sort_keys=True) + "\n")
    return RunArtifacts(panel, scores, batch_state, construct, events)


def main(argv: Iterable[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Run reproducible TAMAVAQ/Q-TAMAVAQ selector")
    parser.add_argument("--candidates", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args(argv)
    run_pipeline(args.candidates, args.config, args.outdir)


if __name__ == "__main__":
    main()
