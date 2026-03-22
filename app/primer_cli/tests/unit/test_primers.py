from __future__ import annotations

from types import SimpleNamespace

import pytest

from primer_cli.services.primers.pair_candidates import (
    PrimerPairingConfig,
    build_candidate_primer_pairs,
)


class _FakeHeterodimer:
    def __init__(self, tm: float) -> None:
        self.tm = tm


def test_pair_candidates_not_empty_when_constraints_allow(monkeypatch: pytest.MonkeyPatch) -> None:
    import primer_cli.services.primers.pair_candidates as pc

    monkeypatch.setattr(pc.primer3.bindings, "calc_heterodimer", lambda *_: _FakeHeterodimer(10.0))

    fwd = SimpleNamespace(
        sequence="ACGTACGTACGTACGTAC",
        orientation="forward",
        msa_start=100,
        msa_end=119,
        tm=60.0,
        gc_percent=50.0,
    )
    rev = SimpleNamespace(
        sequence="TGCATGCATGCATGCATG",
        orientation="reverse",
        msa_start=150,
        msa_end=169,
        tm=60.5,
        gc_percent=52.0,
    )

    pairs = build_candidate_primer_pairs(
        [fwd, rev],
        cfg=PrimerPairingConfig(
            min_amplicon_len=40,
            max_amplicon_len=220,
            preferred_min_amplicon_len=40,
            preferred_max_amplicon_len=120,
            max_tm_diff=2.0,
            max_heterodimer_tm=47.0,
        ),
    )

    assert len(pairs) == 1
    assert pairs[0].amplicon_length == 69
