#!/usr/bin/env bash
set -euo pipefail
python -m tamavaq.cli --ledger data/candidate_ledger.csv --config configs/default.json --out outputs
