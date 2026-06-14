from __future__ import annotations

from types import SimpleNamespace

from primer_cli.cli.commands import pipeline


def test_build_diverse_pair_order_preserves_pre_scored_order() -> None:
    pair_a = SimpleNamespace(forward_seq="AAAA", reverse_seq="TTTT", marker="a")
    pair_b = SimpleNamespace(forward_seq="CCCC", reverse_seq="GGGG", marker="b")
    pair_cov_by_key = {
        ("AAAA", "TTTT"): pair_a,
        ("CCCC", "GGGG"): pair_b,
    }
    pre_scored = [
        SimpleNamespace(forward_seq="CCCC", reverse_seq="GGGG", final_score=99.0),
        SimpleNamespace(forward_seq="AAAA", reverse_seq="TTTT", final_score=98.0),
    ]

    ordered = pipeline._build_diverse_pair_order(pre_scored, pair_cov_by_key)

    assert ordered == [pair_b, pair_a]
