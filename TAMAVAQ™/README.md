# TAMAVAQ/Q-TAMAVAQ Audit Replay

This repository contains a cleaned, executable, and test-covered Python implementation of the deterministic TAMAVAQ/Q-TAMAVAQ audit-first neoantigen candidate replay described in the manuscript and supplementary appendices.

The package reproduces the frozen demonstrator endpoint:

- selected panel: `C01, C02, C03, C04`
- reserve shell: `C05, C06, C07`
- PASS-negative background: `C08-C14`
- `Delta_4 = 0.040`, `HR@4 = 1.00`, `EF@4 = 2.00`

## Claim boundary

This code is a computational replay and reproducibility package. It is **not** clinical, diagnostic, or treatment advice. It does not claim animal, clinical, survival, or patient-efficacy outcomes.

## Repository layout

```text
configs/default.json              Frozen thresholds, scoring weights, panel size, policy label
data/candidate_ledger.csv         Ordered candidate ledger for replay
src/tamavaq/                      Python package
tests/                            Pytest validation suite
scripts/run_replay.sh             One-command replay helper
.github/workflows/ci.yml          GitHub Actions CI
```

## Quick start

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements-dev.txt
pip install -e .
pytest -q
python -m tamavaq.cli --ledger data/candidate_ledger.csv --config configs/default.json --out outputs
```

After execution, outputs are written to:

```text
outputs/replay_results.csv
outputs/replay_summary.json
```

## Method implemented

The code implements:

1. peptide syntax and ledger validation;
2. non-compensatory predicate gates for BIO, PHYS, GEOM, CMC, PROV, and VALID;
3. physical support `S_PHYS` from energy normalization, contact persistence, RMSD*, and RMSF*;
4. CTQW-style transported support proxy `S_CTQW` from GEOM, PHYS, and BIO;
5. final score `S_Q` and selected/reserve/fail state assignment;
6. SHA-256 candidate and manifest hashes for replay governance.

The numerical defaults are intentionally explicit in `configs/default.json` so that any threshold, score-weight, or policy change is visible in version control.

## GitHub submission commands

```bash
git init
git add .
git commit -m "Add TAMAVAQ audit replay pipeline"
git branch -M main
git remote add origin <YOUR_GITHUB_REPOSITORY_URL>
git push -u origin main
```
