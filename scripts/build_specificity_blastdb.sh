#!/usr/bin/env bash
set -euo pipefail

DB_VERSION="${1:-amr_background_v2026_05}"
ROOT="data/blastdb/${DB_VERSION}"

RAW="${ROOT}/work/subjects.raw.fna"
CLEAN="${ROOT}/work/subjects.clean.fna"
MANIFEST="${ROOT}/work/manifest.tsv"
DB_OUT="${ROOT}/db/specificity_subjects"

mkdir -p "${ROOT}/source" "${ROOT}/work" "${ROOT}/db"

if [ ! -s "${RAW}" ]; then
  echo "ERROR: ${RAW} not found or empty."
  echo "Put source FASTA files into ${ROOT}/source and run:"
  echo "cat ${ROOT}/source/*.fna > ${RAW}"
  exit 1
fi

python scripts/normalize_blast_fasta.py \
  --input "${RAW}" \
  --output "${CLEAN}" \
  --manifest "${MANIFEST}"

makeblastdb \
  -in "${CLEAN}" \
  -dbtype nucl \
  -parse_seqids \
  -blastdb_version 5 \
  -title "${DB_VERSION}" \
  -out "${DB_OUT}"

blastdbcmd -db "${DB_OUT}" -info

sha256sum "${CLEAN}" > "${ROOT}/work/subjects.clean.fna.sha256"

cat > "${ROOT}/build_report.json" <<EOF
{
  "db_version": "${DB_VERSION}",
  "clean_fasta": "${CLEAN}",
  "manifest": "${MANIFEST}",
  "blast_db": "${DB_OUT}",
  "dbtype": "nucl"
}
EOF

echo "Built BLAST DB: ${DB_OUT}"
