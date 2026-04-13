# src/primer_cli/cli/commands/align.py
from __future__ import annotations

from pathlib import Path
import shutil
import sys
import time

from primer_cli.core.validation import require_file_exists, require_not_directory, validation_error
from primer_cli.services.aligners.mafft import MafftAligner


def cmd_align(args) -> int:
    in_path = Path(args.inp)
    out_path = Path(args.out)

    require_file_exists(in_path, where="align --input", arg_name="--input")
    require_not_directory(out_path, where="align --output", arg_name="--output")

    mafft_bin = args.mafft
    if shutil.which(mafft_bin) is None:
        raise validation_error(
            what=f"MAFFT executable not found: {mafft_bin}",
            where="align --mafft",
            fix="Install MAFFT or pass a valid executable path via --mafft.",
        )

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
