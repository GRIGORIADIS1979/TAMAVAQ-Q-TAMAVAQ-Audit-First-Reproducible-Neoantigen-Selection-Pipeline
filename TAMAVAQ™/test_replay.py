from pathlib import Path

from tamavaq.io import load_config, load_ledger
from tamavaq.replay import replay

ROOT = Path(__file__).resolve().parents[1]


def test_replay_reproduces_selected_and_reserve_panel() -> None:
    result = replay(
        load_ledger(ROOT / "data" / "candidate_ledger.csv"),
        load_config(ROOT / "configs" / "default.json"),
    )
    assert result.selected_ids == ["C01", "C02", "C03", "C04"]
    assert result.reserve_ids == ["C05", "C06", "C07"]
    assert len(result.fail_ids) == 7
    assert round(result.delta_k, 3) == 0.040
    assert result.hr_at_k == 1.0
    assert result.ef_at_k == 2.0


def test_replay_q_scores_match_reported_rounded_values() -> None:
    result = replay(
        load_ledger(ROOT / "data" / "candidate_ledger.csv"),
        load_config(ROOT / "configs" / "default.json"),
    )
    by_id = {item.candidate.id: item for item in result.results}
    for candidate_id in ["C01", "C02", "C03", "C04", "C05", "C06", "C07"]:
        reported = by_id[candidate_id].candidate.reported_q_score
        assert reported is not None
        assert round(by_id[candidate_id].s_q, 3) == reported


def test_background_records_fail_non_compensatory_predicate() -> None:
    result = replay(
        load_ledger(ROOT / "data" / "candidate_ledger.csv"),
        load_config(ROOT / "configs" / "default.json"),
    )
    by_id = {item.candidate.id: item for item in result.results}
    assert all(by_id[f"C{i:02d}"].pass_bit == 0 for i in range(8, 15))
    assert by_id["C08"].margins["bio"] < 0
    assert by_id["C08"].margins["geom"] < 0
