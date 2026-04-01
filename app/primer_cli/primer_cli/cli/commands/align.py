from __future__ import annotations

import logging
import shutil
import time
from pathlib import Path

from primer_cli.core.exceptions import PrimerCliError
from primer_cli.services.aligners.mafft import MafftAligner

logger = logging.getLogger(__name__)


def cmd_align(args) -> int:
    input_path = getattr(args, "input_path", None) or getattr(args, "inp", None)
    output_path = getattr(args, "output", None) or getattr(args, "out", None)
    if not input_path:
        raise PrimerCliError("--input is required")
    if not output_path:
        raise PrimerCliError("--output is required")
    in_path = Path(input_path)
    out_path = Path(output_path)

    if not in_path.exists():
        raise PrimerCliError(f"Input FASTA does not exist: {in_path}")

    if out_path.exists() and out_path.is_dir():
        raise PrimerCliError(f"Output path is a directory, expected file: {out_path}")

    mafft_bin = args.mafft
    if shutil.which(mafft_bin) is None:
        raise PrimerCliError(f"MAFFT binary not found in PATH: {mafft_bin}")

    out_path.parent.mkdir(parents=True, exist_ok=True)

    aligner = MafftAligner(binary=mafft_bin)
    logger.info("Starting alignment with MAFFT: %s -> %s", in_path, out_path)
    t0 = time.monotonic()
    aligner.align_fasta(
        input_path=str(in_path),
        output_path=str(out_path),
        extra_args=args.mafft_args,
    )
    dt = time.monotonic() - t0
    logger.info("Alignment completed in %.1fs", dt)

    return 0
