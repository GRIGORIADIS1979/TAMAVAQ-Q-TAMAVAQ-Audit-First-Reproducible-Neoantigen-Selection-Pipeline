"""Validation rules for TAMAVAQ candidate ledgers."""

from __future__ import annotations

import re

from .models import Candidate

AMINO_ACID_PATTERN = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")


def validate_candidates(candidates: list[Candidate]) -> None:
    if not candidates:
        raise ValueError("Ledger is empty.")
    ids = [candidate.id for candidate in candidates]
    duplicates = sorted({value for value in ids if ids.count(value) > 1})
    if duplicates:
        raise ValueError(f"Duplicate candidate IDs found: {duplicates}")
    for candidate in candidates:
        _validate_candidate(candidate)


def _validate_candidate(candidate: Candidate) -> None:
    if not candidate.id:
        raise ValueError("Candidate has an empty id.")
    if not AMINO_ACID_PATTERN.fullmatch(candidate.peptide):
        raise ValueError(f"Candidate {candidate.id} has invalid peptide syntax: {candidate.peptide!r}")
    if not 5 <= len(candidate.peptide) <= 35:
        raise ValueError(f"Candidate {candidate.id} length outside supported range 5-35.")
    for name, value in {
        "bio": candidate.bio,
        "contact_persistence": candidate.contact_persistence,
        "rmsd_star": candidate.rmsd_star,
        "rmsf_star": candidate.rmsf_star,
        "geom": candidate.geom,
    }.items():
        if not 0.0 <= value <= 1.0:
            raise ValueError(f"Candidate {candidate.id} field {name} must be in [0, 1].")
    for name, value in {"cmc": candidate.cmc, "prov": candidate.prov, "valid": candidate.valid}.items():
        if value not in (0, 1):
            raise ValueError(f"Candidate {candidate.id} field {name} must be 0 or 1.")
