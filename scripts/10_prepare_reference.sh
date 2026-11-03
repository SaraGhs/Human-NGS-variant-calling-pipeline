#!/usr/bin/env bash
# Validate an existing putative GRCh38/hg38 reference and link its index files
# into the project WITHOUT altering the original data on /mnt/d.
# Usage: bash scripts/10_prepare_reference.sh [REFERENCE_FASTA]
set -euo pipefail
cd "$(dirname "$0")/.."
ORIGINAL=${1:-/mnt/d/Udemy/NGS/hg38_human_genome/hg38.fa}
[[ -s "$ORIGINAL" ]] || { echo "Missing reference: $ORIGINAL" >&2; exit 1; }
for suffix in fai amb ann bwt pac sa; do
  [[ -s "${ORIGINAL}.${suffix}" ]] || { echo "Missing BWA/SAMtools index: ${ORIGINAL}.${suffix}" >&2; exit 1; }
done
for cmd in awk gatk readlink; do command -v "$cmd" >/dev/null || { echo "Missing: $cmd" >&2; exit 1; }; done
chr1=$(awk -F '\t' '$1=="chr1" {print $2; exit}' "${ORIGINAL}.fai")
chr19=$(awk -F '\t' '$1=="chr19" {print $2; exit}' "${ORIGINAL}.fai")
[[ "$chr1" == 248956422 && "$chr19" == 58617616 ]] || {
  echo "STOP: Reference contig lengths do not match GRCh38 primary chromosome 1/19." >&2
  echo "Observed chr1=${chr1:-absent} chr19=${chr19:-absent}" >&2
  exit 2
}
mkdir -p data/reference docs
REF=data/reference/hg38.fa
for suffix in '' .fai .amb .ann .bwt .pac .sa; do
  src="${ORIGINAL}${suffix}"; dst="${REF}${suffix}"
  if [[ -e "$dst" || -L "$dst" ]]; then
    [[ -L "$dst" && "$(readlink -f "$dst")" == "$(readlink -f "$src")" ]] || {
      echo "STOP: Existing destination differs: $dst" >&2; exit 2;
    }
  else
    ln -s "$(readlink -f "$src")" "$dst"
  fi
done
DICT=data/reference/hg38.dict
if [[ ! -s "$DICT" ]]; then
  [[ ! -e "$DICT" ]] || { echo "Dictionary exists but is empty: $DICT" >&2; exit 2; }
  gatk --java-options '-Xmx2g' CreateSequenceDictionary -R "$REF" -O "$DICT"
fi
printf 'reference_original\t%s\nreference_project\t%s\nchr1_length\t%s\nchr19_length\t%s\nreference_provenance\t%s\n' \
  "$ORIGINAL" "$REF" "$chr1" "$chr19" 'Source URL/hash of original course reference still to be documented' \
  > docs/reference_check.tsv
printf 'Reference checks passed. Project symlinks and local dictionary ready.\n'
ls -lh "$REF" "$REF.fai" "$DICT"
