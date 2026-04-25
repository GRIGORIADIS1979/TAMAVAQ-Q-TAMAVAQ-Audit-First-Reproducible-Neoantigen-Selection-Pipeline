#!/usr/bin/env bash
set -euo pipefail
PYTHONPATH=src python -m tamavaq.pipeline --candidates examples/candidates.csv --config examples/config.yaml --outdir runs/example
