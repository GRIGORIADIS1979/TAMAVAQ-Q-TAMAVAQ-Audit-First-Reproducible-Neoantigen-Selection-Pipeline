"""Command-line entry point for TAMAVAQ replay."""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path

from .io import load_config, load_ledger, write_csv
from .replay import replay


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Replay TAMAVAQ/Q-TAMAVAQ deterministic panel selection.")
    parser.add_argument("--ledger", default="data/candidate_ledger.csv", help="Input candidate CSV ledger.")
    parser.add_argument("--config", default="configs/default.json", help="Frozen replay configuration JSON.")
    parser.add_argument("--out", default="outputs", help="Output directory.")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    config = load_config(args.config)
    candidates = load_ledger(args.ledger)
    result = replay(candidates, config)

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    write_csv([item.as_dict() for item in result.results], out_dir / "replay_results.csv")
    summary = {
        "selected_ids": result.selected_ids,
        "reserve_ids": result.reserve_ids,
        "fail_ids": result.fail_ids,
        "delta_k": round(result.delta_k, 6),
        "hr_at_k": round(result.hr_at_k, 6),
        "ef_at_k": round(result.ef_at_k, 6),
        "manifest_hash": result.manifest_hash,
        "claim_boundary": config.claim_boundary,
        "config": asdict(config),
    }
    (out_dir / "replay_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
