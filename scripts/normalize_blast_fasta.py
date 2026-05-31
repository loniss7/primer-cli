from __future__ import annotations

import argparse
import re
from pathlib import Path


def clean_id(raw: str, idx: int) -> str:
    first = raw.strip().split()[0] if raw.strip() else f"seq_{idx}"
    first = first.replace("|", "_")
    first = re.sub(r"[^A-Za-z0-9_.:-]+", "_", first)
    return first or f"seq_{idx}"


def normalize_fasta(input_path: Path, output_path: Path, manifest_path: Path) -> None:
    seen: dict[str, int] = {}
    seq_idx = 0

    output_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)

    with (
        input_path.open("r", encoding="utf-8") as fin,
        output_path.open("w", encoding="utf-8") as fout,
        manifest_path.open("w", encoding="utf-8") as mf,
    ):
        mf.write("new_id\told_header\n")

        for line in fin:
            line = line.rstrip("\n")
            if not line:
                continue

            if line.startswith(">"):
                seq_idx += 1
                old_header = line[1:].strip()
                base_id = clean_id(old_header, seq_idx)

                n = seen.get(base_id, 0)
                seen[base_id] = n + 1
                new_id = base_id if n == 0 else f"{base_id}_{n + 1}"

                fout.write(f">{new_id}\n")
                mf.write(f"{new_id}\t{old_header}\n")
            else:
                fout.write(line.strip().upper() + "\n")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    args = parser.parse_args()

    normalize_fasta(args.input, args.output, args.manifest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
