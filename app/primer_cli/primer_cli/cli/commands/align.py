# src/primer_cli/cli/commands/align.py
from __future__ import annotations

import logging
from pathlib import Path
import time

from primer_cli.core.validation import require_file_exists, require_not_directory, validation_error
from primer_cli.services.aligners.mafft import MafftAligner

logger = logging.getLogger(__name__)


def cmd_align(args) -> int:
    in_path = Path(args.inp)
    out_path = Path(args.out)

    require_file_exists(in_path, where="align --input", arg_name="--input")
    require_not_directory(out_path, where="align --output", arg_name="--output")

    mafft_bin = args.mafft

    out_path.parent.mkdir(parents=True, exist_ok=True)

    aligner = MafftAligner(binary=mafft_bin)
    logger.info(
        "Starting alignment with MAFFT: input=%s output=%s binary=%s args=%s",
        in_path,
        out_path,
        mafft_bin,
        args.mafft_args,
    )
    t0 = time.monotonic()
    aligner.align_fasta(
        input_path=str(in_path),
        output_path=str(out_path),
        extra_args=args.mafft_args,
    )
    dt = time.monotonic() - t0
    logger.info("Alignment completed in %.1fs: %s", dt, out_path)

    return 0
