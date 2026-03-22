# src/primer_cli/cli/commands/align.py
from __future__ import annotations

from pathlib import Path
import shutil
import sys
import time

from primer_cli.core.exceptions import PrimerCliError
from primer_cli.services.aligners.mafft import MafftAligner


def cmd_align(args) -> int:
    """
    Align FASTA sequences using MAFFT.
    CLI orchestration only: validation + service invocation.
    """
    in_path = Path(args.inp)
    out_path = Path(args.out)

    if not in_path.exists():
        raise PrimerCliError(f"Input FASTA does not exist: {in_path}")

    if out_path.exists():
        raise PrimerCliError(f"Output file already exists: {out_path}")

    mafft_bin = args.mafft
    if shutil.which(mafft_bin) is None:
        raise PrimerCliError(f"MAFFT binary not found in PATH: {mafft_bin}")

    out_path.parent.mkdir(parents=True, exist_ok=True)

    aligner = MafftAligner(binary=mafft_bin)
    print(
        f"Starting alignment with MAFFT: {in_path} -> {out_path}",
        file=sys.stderr,
        flush=True,
    )
    t0 = time.monotonic()
    aligner.align_fasta(
        input_path=str(in_path),
        output_path=str(out_path),
        extra_args=args.mafft_args,
    )
    dt = time.monotonic() - t0
    print(f"Alignment completed in {dt:.1f}s", file=sys.stderr, flush=True)

    return 0
