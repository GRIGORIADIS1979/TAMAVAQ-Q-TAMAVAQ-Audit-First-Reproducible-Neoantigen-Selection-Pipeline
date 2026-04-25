from pathlib import Path

from tamavaq.pipeline import run_pipeline


def test_reproducible_selection(tmp_path):
    root = Path(__file__).resolve().parents[1]
    a = run_pipeline(root / "examples" / "candidates.csv", root / "examples" / "config.yaml", tmp_path / "a")
    b = run_pipeline(root / "examples" / "candidates.csv", root / "examples" / "config.yaml", tmp_path / "b")
    assert a.batch_state["candidate_index_map_sha256"] == b.batch_state["candidate_index_map_sha256"]
    assert a.batch_state["selected_ids"] == b.batch_state["selected_ids"]
    assert a.construct_sequence == b.construct_sequence
