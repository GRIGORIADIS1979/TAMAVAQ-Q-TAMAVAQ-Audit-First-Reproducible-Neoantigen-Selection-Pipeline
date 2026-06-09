"""Typed data structures for the TAMAVAQ replay pipeline."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any


@dataclass(frozen=True)
class Candidate:
    """Immutable input dossier for one peptide candidate."""

    id: str
    peptide: str
    role: str
    bio: float
    delta_g_kcal_mol: float
    contact_persistence: float
    rmsd_star: float
    rmsf_star: float
    geom: float
    cmc: int
    prov: int
    valid: int
    reported_status: str | None = None
    reported_q_score: float | None = None


@dataclass(frozen=True)
class Config:
    """Frozen replay configuration."""

    panel_size: int
    thresholds: dict[str, float]
    weights: dict[str, dict[str, float]]
    policy_version: str
    claim_boundary: str


@dataclass(frozen=True)
class CandidateResult:
    """Computed replay state for one candidate."""

    candidate: Candidate
    s_phys: float
    s_ctqw: float
    s_q: float
    pass_bit: int
    state: str
    rank: int
    margins: dict[str, float]
    hash_id: str

    def as_dict(self) -> dict[str, Any]:
        data: dict[str, Any] = {
            "id": self.candidate.id,
            "peptide": self.candidate.peptide,
            "role": self.candidate.role,
            "BIO": round(self.candidate.bio, 6),
            "S_PHYS": round(self.s_phys, 6),
            "GEOM": round(self.candidate.geom, 6),
            "S_CTQW": round(self.s_ctqw, 6),
            "S_Q": round(self.s_q, 6),
            "PASS": self.pass_bit,
            "rank": self.rank,
            "state": self.state,
            "hash": self.hash_id,
        }
        data.update({f"margin_{key}": round(value, 6) for key, value in self.margins.items()})
        return data


@dataclass(frozen=True)
class ReplayResult:
    """Full replay result and manifest certificate."""

    results: list[CandidateResult]
    selected_ids: list[str]
    reserve_ids: list[str]
    fail_ids: list[str]
    delta_k: float
    hr_at_k: float
    ef_at_k: float
    manifest_hash: str
