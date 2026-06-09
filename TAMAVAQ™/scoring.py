"""Deterministic scoring functions for TAMAVAQ replay."""

from __future__ import annotations

from .models import Candidate, Config


def normalized_energy(candidate: Candidate, candidates: list[Candidate]) -> float:
    abs_values = [abs(item.delta_g_kcal_mol) for item in candidates]
    lower = min(abs_values)
    upper = max(abs_values)
    if upper == lower:
        return 0.0
    return (abs(candidate.delta_g_kcal_mol) - lower) / (upper - lower)


def s_phys(candidate: Candidate, candidates: list[Candidate]) -> float:
    e_norm = normalized_energy(candidate, candidates)
    return (
        0.35 * e_norm
        + 0.35 * candidate.contact_persistence
        + 0.15 * (1.0 - candidate.rmsd_star)
        + 0.15 * (1.0 - candidate.rmsf_star)
    )


def s_ctqw(candidate: Candidate, phys_score: float, config: Config) -> float:
    weights = config.weights["s_ctqw"]
    return (
        weights["geom"] * candidate.geom
        + weights["s_phys"] * phys_score
        + weights["bio"] * candidate.bio
    )


def pass_bit(candidate: Candidate, config: Config) -> int:
    t = config.thresholds
    checks = [
        candidate.bio >= t["bio"],
        candidate.delta_g_kcal_mol <= t["delta_g_kcal_mol_max"],
        candidate.contact_persistence >= t["contact_persistence"],
        candidate.rmsd_star <= t["rmsd_star_max"],
        candidate.rmsf_star <= t["rmsf_star_max"],
        candidate.geom >= t["geom"],
        candidate.cmc >= int(t["cmc"]),
        candidate.prov >= int(t["prov"]),
        candidate.valid >= int(t["valid"]),
    ]
    return int(all(checks))


def s_q(candidate: Candidate, phys_score: float, ctqw_score: float, pass_value: int, config: Config) -> float:
    weights = config.weights["s_q"]
    return (
        weights["bio"] * candidate.bio
        + weights["s_phys"] * phys_score
        + weights["geom"] * candidate.geom
        + weights["s_ctqw"] * ctqw_score
        + weights["pass"] * pass_value
    )


def margins(candidate: Candidate, config: Config) -> dict[str, float]:
    t = config.thresholds
    return {
        "bio": candidate.bio - t["bio"],
        "delta_g": t["delta_g_kcal_mol_max"] - candidate.delta_g_kcal_mol,
        "contact_persistence": candidate.contact_persistence - t["contact_persistence"],
        "rmsd_star": t["rmsd_star_max"] - candidate.rmsd_star,
        "rmsf_star": t["rmsf_star_max"] - candidate.rmsf_star,
        "geom": candidate.geom - t["geom"],
        "cmc": candidate.cmc - t["cmc"],
        "prov": candidate.prov - t["prov"],
        "valid": candidate.valid - t["valid"],
    }
