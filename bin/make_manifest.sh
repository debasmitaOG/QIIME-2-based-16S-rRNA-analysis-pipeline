#!/usr/bin/env bash
# Builds a QIIME 2 paired-end manifest from the FASTQ files in the current directory.
# Expects Illumina-style names: <sample>_<anything>_R1_001.fastq / ..._R2_001.fastq
# Sample ID = text before the first underscore (F3D0_S188_L001_R1_001.fastq -> F3D0).
# gzips copies of the reads so the manifest points at real files inside the task folder.
set -euo pipefail

out="${1:-manifest.tsv}"
printf 'sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n' > "$out"

for r1 in *_R1_001.fastq; do
    r2="${r1/_R1_/_R2_}"
    sid="${r1%%_*}"
    gzip -c "$r1" > "${r1}.gz"
    gzip -c "$r2" > "${r2}.gz"
    printf '%s\t%s/%s.gz\t%s/%s.gz\n' "$sid" "$PWD" "$r1" "$PWD" "$r2" >> "$out"
done
