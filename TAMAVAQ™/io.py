"""Input/output helpers for CSV ledgers and JSON configs."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

from .models import Candidate, Config


REQUIRED_COLUMNS = {
    "id",
    "peptide",
    "role",
    "bio",
    "delta_g_kcal_mol",
    "contact_persistence",
    "rmsd_star",
    "rmsf_star",
    "geom",
    "cmc",
    "prov",
    "valid",
}


def load_config(path: str | Path) -> Config:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    return Config(
        panel_size=int(payload["panel_size"]),
        thresholds={k: float(v) for k, v in payload["thresholds"].items()},
        weights={k: {kk: float(vv) for kk, vv in v.items()} for k, v in payload["weights"].items()},
        policy_version=str(payload["policy_version"]),
        claim_boundary=str(payload["claim_boundary"]),
    )


def load_ledger(path: str | Path) -> list[Candidate]:
    with Path(path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        columns = set(reader.fieldnames or [])
        missing = REQUIRED_COLUMNS - columns
        if missing:
            raise ValueError(f"Ledger is missing required columns: {sorted(missing)}")
        return [_row_to_candidate(row) for row in reader]


def _row_to_candidate(row: dict[str, str]) -> Candidate:
    return Candidate(
        id=row["id"].strip(),
        peptide=row["peptide"].strip().upper(),
        role=row["role"].strip(),
        bio=float(row["bio"]),
        delta_g_kcal_mol=float(row["delta_g_kcal_mol"]),
        contact_persistence=float(row["contact_persistence"]),
        rmsd_star=float(row["rmsd_star"]),
        rmsf_star=float(row["rmsf_star"]),
        geom=float(row["geom"]),
        cmc=int(row["cmc"]),
        prov=int(row["prov"]),
        valid=int(row["valid"]),
        reported_status=(row.get("reported_status") or None),
        reported_q_score=float(row["reported_q_score"]) if row.get("reported_q_score") else None,
    )


def write_csv(rows: list[dict[str, Any]], path: str | Path) -> None:
    if not rows:
        raise ValueError("Cannot write an empty result table.")
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
