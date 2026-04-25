from pathlib import Path
import json

from tamavaq.pipeline import run_pipeline


def test_pipeline_outputs(tmp_path):
    root = Path(__file__).resolve().parents[1]
    artifacts = run_pipeline(root / "examples" / "candidates.csv", root / "examples" / "config.yaml", tmp_path)
    assert (tmp_path / "scores.csv").exists()
    assert (tmp_path / "selected_panel.csv").exists()
    assert (tmp_path / "batch_state.json").exists()
    assert (tmp_path / "teleport_log.jsonl").exists()
    assert len(artifacts.selected_panel) <= 4
    assert "scope" in artifacts.batch_state
    assert artifacts.batch_state["predicate_counts"]["MARK"] >= len(artifacts.selected_panel)
