"""High-level deterministic replay procedure."""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict

from .models import Candidate, CandidateResult, Config, ReplayResult
from .scoring import margins, pass_bit, s_ctqw, s_phys, s_q
from .validation import validate_candidates


def replay(candidates: list[Candidate], config: Config) -> ReplayResult:
    """Replay candidate validation, scoring, panel selection, and checksums."""
    validate_candidates(candidates)

    interim: list[tuple[Candidate, float, float, float, int, dict[str, float]]] = []
    for candidate in candidates:
        phys = s_phys(candidate, candidates)
        pass_value = pass_bit(candidate, config)
        ctqw = s_ctqw(candidate, phys, config)
        q_score = s_q(candidate, phys, ctqw, pass_value, config)
        interim.append((candidate, phys, ctqw, q_score, pass_value, margins(candidate, config)))

    ranked = sorted(interim, key=lambda row: (row[4], row[3], row[0].id), reverse=True)
    selected_ids = [row[0].id for row in ranked if row[4] == 1][: config.panel_size]
    reserve_ids = [row[0].id for row in ranked if row[4] == 1][config.panel_size :]
    fail_ids = [row[0].id for row in ranked if row[4] == 0]

    results: list[CandidateResult] = []
    for rank, (candidate, phys, ctqw, q_score, pass_value, margin_values) in enumerate(ranked, start=1):
        if candidate.id in selected_ids:
            state = "selected"
        elif candidate.id in reserve_ids:
            state = "reserve"
        else:
            state = "fail"
        digest = _candidate_hash(candidate, q_score, pass_value, state, margin_values, config)
        results.append(
            CandidateResult(
                candidate=candidate,
                s_phys=phys,
                s_ctqw=ctqw,
                s_q=q_score,
                pass_bit=pass_value,
                state=state,
                rank=rank,
                margins=margin_values,
                hash_id=digest,
            )
        )

    delta_k = _delta_k(results, config.panel_size)
    hr_at_k = sum(item.pass_bit for item in results[: config.panel_size]) / config.panel_size
    positives = sum(item.pass_bit for item in results)
    ef_at_k = hr_at_k / (positives / len(results)) if positives else 0.0
    manifest_hash = _manifest_hash(results, config)
    return ReplayResult(
        results=results,
        selected_ids=selected_ids,
        reserve_ids=reserve_ids,
        fail_ids=fail_ids,
        delta_k=delta_k,
        hr_at_k=hr_at_k,
        ef_at_k=ef_at_k,
        manifest_hash=manifest_hash,
    )


def _delta_k(results: list[CandidateResult], panel_size: int) -> float:
    if len(results) <= panel_size:
        return 0.0
    return results[panel_size - 1].s_q - results[panel_size].s_q


def _candidate_hash(
    candidate: Candidate,
    q_score: float,
    pass_value: int,
    state: str,
    margin_values: dict[str, float],
    config: Config,
) -> str:
    payload = {
        "candidate": asdict(candidate),
        "q_score": round(q_score, 12),
        "pass": pass_value,
        "state": state,
        "margins": {key: round(value, 12) for key, value in margin_values.items()},
        "policy_version": config.policy_version,
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode("utf-8")).hexdigest()


def _manifest_hash(results: list[CandidateResult], config: Config) -> str:
    payload = {
        "config": asdict(config),
        "state_vector": [
            {
                "id": item.candidate.id,
                "rank": item.rank,
                "state": item.state,
                "s_q": round(item.s_q, 12),
                "hash": item.hash_id,
            }
            for item in results
        ],
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode("utf-8")).hexdigest()
