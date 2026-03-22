from __future__ import annotations

from primer_cli.services.aligners.mafft import _MafftProgress


def test_mafft_progress_phase_progress_increases() -> None:
    p = _MafftProgress(phase_idx=1, phase_total=2)

    start = p._phase_progress(0, 10)
    mid = p._phase_progress(5, 10)
    end = p._phase_progress(10, 10)

    assert start < mid < end
