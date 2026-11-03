#!/usr/bin/env bash
# Safe, explicitly gated quality trimming of the three original FASTQ pairs.
# Usage: QC_REVIEWED=yes bash scripts/07_trim_three_pairs.sh quality-only [THREADS]
# Or:    QC_REVIEWED=yes bash scripts/07_trim_three_pairs.sh adapters [THREADS] PATH_TO_ADAPTER_FASTA
set -euo pipefail
cd "$(dirname "$0")/.."
[[ "${QC_REVIEWED:-no}" == yes ]] || { echo 'STOP: review docs/ADAPTER_REVIEW.md first; set QC_REVIEWED=yes only after review.' >&2; exit 2; }
MODE=${1:-}; THREADS=${2:-2}; ADAPTERS=${3:-}
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || { echo 'THREADS must be positive integer' >&2; exit 2; }
[[ "$MODE" == 'quality-only' || "$MODE" == 'adapters' ]] || { echo 'MODE must be quality-only or adapters' >&2; exit 2; }
if [[ "$MODE" == 'adapters' ]]; then
  [[ -s "$ADAPTERS" ]] || { echo 'Missing adapter FASTA' >&2; exit 2; }
fi
command -v trimmomatic >/dev/null || { echo 'Activate ngs-variant-portfolio environment' >&2; exit 1; }
command -v fastqc >/dev/null || { echo 'fastqc missing' >&2; exit 1; }
command -v multiqc >/dev/null || { echo 'multiqc missing' >&2; exit 1; }
mkdir -p results/trimmed results/qc/trimmed results/reports logs
shopt -s nullglob
pairs=(data/source/NA12878_20k/*_reads_1.fastq)
[[ "${#pairs[@]}" == 3 ]] || { echo "Expected three R1 FASTQ files, found ${#pairs[@]}" >&2; exit 1; }
for r1 in "${pairs[@]}"; do
  r2=${r1/_reads_1.fastq/_reads_2.fastq}
  [[ -s "$r2" ]] || { echo "Missing mate: $r2" >&2; exit 1; }
  group=$(basename "$r1" _reads_1.fastq)
  r1p="results/trimmed/${group}_R1_paired.fastq.gz"
  r2p="results/trimmed/${group}_R2_paired.fastq.gz"
  r1u="results/trimmed/${group}_R1_unpaired.fastq.gz"
  r2u="results/trimmed/${group}_R2_unpaired.fastq.gz"
  for f in "$r1p" "$r2p" "$r1u" "$r2u"; do
    [[ ! -e "$f" ]] || { echo "Refusing to overwrite existing output: $f" >&2; exit 2; }
  done
  clip=()
  if [[ "$MODE" == 'adapters' ]]; then clip=("ILLUMINACLIP:${ADAPTERS}:2:30:10"); fi
  echo "Trimming group $group ($MODE)"
  trimmomatic PE -threads "$THREADS" -phred33 \
    "$r1" "$r2" "$r1p" "$r1u" "$r2p" "$r2u" \
    "${clip[@]}" LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50 \
    2>&1 | tee "logs/trimmomatic_${group}.log"
  fastqc --threads "$THREADS" --outdir results/qc/trimmed "$r1p" "$r2p"
done
multiqc results/qc/trimmed -o results/reports -n multiqc_trimmed.html -f
printf 'Done. Next run python scripts/08_summarize_trimming.py\n'
