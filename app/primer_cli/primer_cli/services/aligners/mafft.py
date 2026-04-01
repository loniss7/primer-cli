from __future__ import annotations

import re
import subprocess
import sys
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence

from primer_cli.core.exceptions import PrimerCliError
from primer_cli.utils.subprocess import run_cmd


_COUNT_RE = re.compile(r"^\s*(\d+)\s*/\s*(\d+)\s*$")
_STEP_RE = re.compile(r"STEP\s+(\d+)\s*/\s*(\d+)")
_PHASE_RE = re.compile(r"Progressive alignment\s+(\d+)\s*/\s*(\d+)\.\.\.")


@dataclass
class _MafftProgress:
    width: int = 34
    phase_idx: int = 1
    phase_total: int = 1
    last_render_len: int = 0

    def _render(self, percent: float, status: str) -> None:
        p = min(max(percent, 0.0), 100.0)
        filled = int(self.width * (p / 100.0))
        bar = "#" * filled + "-" * (self.width - filled)
        line = f"[mafft] [{bar}] {p:5.1f}% | {status}"
        pad = " " * max(0, self.last_render_len - len(line))
        print(f"\r{line}{pad}", end="", file=sys.stderr, flush=True)
        self.last_render_len = len(line)

    def _phase_progress(self, step: int, total: int) -> float:
        # 30% reserved for preprocessing, 70% for progressive alignment.
        phase_span = 70.0 / max(self.phase_total, 1)
        phase_base = 30.0 + (self.phase_idx - 1) * phase_span
        return phase_base + phase_span * (step / max(total, 1))

    def handle(self, line: str) -> None:
        if not line:
            return

        phase_match = _PHASE_RE.search(line)
        if phase_match:
            self.phase_idx = int(phase_match.group(1))
            self.phase_total = int(phase_match.group(2))
            self._render(self._phase_progress(0, 1), f"Progressive phase {self.phase_idx}/{self.phase_total}")
            return

        if "generating a scoring matrix" in line:
            self._render(5.0, "Generating scoring matrix")
            return
        if "Gap Penalty" in line:
            self._render(10.0, "Applying gap penalties")
            return
        if "Making a distance matrix" in line:
            self._render(12.0, "Building distance matrix")
            return
        if "Constructing a UPGMA tree" in line:
            self._render(20.0, "Constructing guide tree")
            return

        step_match = _STEP_RE.search(line)
        if step_match:
            step = int(step_match.group(1))
            total = int(step_match.group(2))
            pct = self._phase_progress(step, total)
            self._render(pct, f"Progressive alignment step {step}/{total}")
            return

        count_match = _COUNT_RE.match(line)
        if count_match:
            cur = int(count_match.group(1))
            total = int(count_match.group(2))
            if total > 0:
                if cur <= total and cur < max(3, total // 2):
                    pct = 12.0 + 8.0 * (cur / total)
                    self._render(pct, f"Distance matrix {cur}/{total}")
                else:
                    pct = 20.0 + 10.0 * (cur / total)
                    self._render(pct, f"Guide tree {cur}/{total}")
            return

        if line == "done.":
            self._render(30.0, "Preprocessing done")
            return

    def finish(self) -> None:
        self._render(100.0, "Alignment completed")
        print(file=sys.stderr, flush=True)


@dataclass
class MafftAligner:

    binary: str = "mafft"

    def _ensure_binary(self) -> None:
        if shutil.which(self.binary) is None:
            raise PrimerCliError(f"MAFFT binary not found in PATH: {self.binary}")

    def version(self) -> str:
        self._ensure_binary()
        res = run_cmd([self.binary, "--version"], capture_stdout=True)
        return res.stdout.strip()

    def align_fasta(
        self,
        input_path: str | Path,
        output_path: str | Path,
        extra_args: Optional[Sequence[str] | str] = None,
    ) -> None:

        self._ensure_binary()

        in_path = Path(input_path)
        out_path = Path(output_path)

        if not in_path.exists():
            raise PrimerCliError(f"Input FASTA does not exist: {in_path}")

        out_path.parent.mkdir(parents=True, exist_ok=True)

        if extra_args is None:
            args: list[str] = []
        elif isinstance(extra_args, str):
            args = extra_args.split()
        else:
            args = list(extra_args)

        cmd = [self.binary, *args, str(in_path)]
        stderr_lines: list[str] = []
        progress = _MafftProgress()
        try:
            with out_path.open("w", encoding="utf-8") as out_f:
                proc = subprocess.Popen(
                    cmd,
                    text=True,
                    stdout=out_f,
                    stderr=subprocess.PIPE,
                )
                assert proc.stderr is not None
                for raw in proc.stderr:
                    line = raw.rstrip()
                    if not line:
                        continue
                    stderr_lines.append(line)
                    progress.handle(line)
                rc = proc.wait()
        except FileNotFoundError as e:
            raise PrimerCliError(f"Executable not found: {cmd[0]}") from e
        except Exception as e:
            raise PrimerCliError(f"Failed to run command: {' '.join(cmd)}") from e

        if rc != 0:
            print(file=sys.stderr, flush=True)
            tail = "\n".join(stderr_lines[-10:])
            raise PrimerCliError(f"Command failed ({rc}): {' '.join(cmd)}\n{tail}".rstrip())

        progress.finish()

        out_raw = out_path.read_text(encoding="utf-8")
        if not out_raw.strip():
            raise PrimerCliError("MAFFT produced empty output")

        lines = out_raw.splitlines()
        processed_lines: list[str] = []

        for line in lines:
            if line.startswith(">"):
                processed_lines.append(line)
            else:
                processed_lines.append(line.upper())

        out_text = "\n".join(processed_lines) + "\n"

        out_path.write_text(out_text, encoding="utf-8")
