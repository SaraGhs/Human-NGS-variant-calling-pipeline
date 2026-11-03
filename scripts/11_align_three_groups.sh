#!/usr/bin/env bash
# Map the three verified Trimmomatic paired-end outputs to the complete reference.
# BWA uses one thread to reduce RAM demand on the 7.5 GiB WSL setup.
# Usage: bash scripts/11_align_three_groups.sh
set -euo pipefail
cd "$(dirname "$0")/.."
REF=data/reference/hg38.fa
[[ -s "$REF" && -s "$REF.fai" && -s data/reference/hg38.dict ]] || {
  echo 'STOP: run scripts/10_prepare_reference.sh and verify provenance first' >&2; exit 2;
}
for cmd in bwa samtools; do command -v "$cmd" >/dev/null || { echo "Missing $cmd" >&2; exit 2; }; done
mkdir -p results/alignments results/reports logs docs
shopt -s nullglob
r1s=(results/trimmed/*_R1_paired.fastq.gz)
[[ ${#r1s[@]} -eq 3 ]] || { echo "STOP: expected three paired R1 FASTQ files, found ${#r1s[@]}" >&2; exit 2; }
FINAL=results/alignments/NA12878.markdup.bam
# Never silently overwrite existing alignments or intermediary outputs.
for f in results/alignments/*.bam "$FINAL.bai" results/reports/alignment_flagstat.txt; do
  [[ ! -e "$f" ]] || { echo "STOP: output already exists: $f. Review it before rerunning." >&2; exit 2; }
done
coordinate_bams=()
for r1 in "${r1s[@]}"; do
  r2=${r1/_R1_paired.fastq.gz/_R2_paired.fastq.gz}
  [[ -s "$r2" ]] || { echo "STOP: missing paired mate: $r2" >&2; exit 2; }
  group=$(basename "$r1" _R1_paired.fastq.gz)
  bam="results/alignments/${group}.coordinate.bam"
  rg="@RG\\tID:${group}\\tSM:NA12878\\tLB:NA12878_WGS\\tPL:ILLUMINA\\tPU:${group}"
  echo "[align] $group"
  bwa mem -t 1 -R "$rg" "$REF" "$r1" "$r2" \
    2> "logs/bwa_${group}.log" | \
    samtools sort -@ 1 -m 256M -o "$bam" -
  samtools quickcheck -v "$bam"
  coordinate_bams+=("$bam")
done
samtools merge -f -@ 1 -o results/alignments/merged.coordinate.bam "${coordinate_bams[@]}"
samtools sort -n -@ 1 -m 256M -o results/alignments/merged.name.bam results/alignments/merged.coordinate.bam
samtools fixmate -m results/alignments/merged.name.bam results/alignments/merged.fixmate.bam
samtools sort -@ 1 -m 256M -o results/alignments/merged.fixmate.coordinate.bam results/alignments/merged.fixmate.bam
samtools markdup -s -f results/reports/duplicate_statistics.txt \
  results/alignments/merged.fixmate.coordinate.bam "$FINAL" 2> logs/markdup.log
samtools index "$FINAL"
samtools quickcheck -v "$FINAL"
samtools flagstat "$FINAL" > results/reports/alignment_flagstat.txt
samtools idxstats "$FINAL" > results/reports/alignment_idxstats.tsv
samtools stats "$FINAL" > results/reports/alignment_stats.txt
samtools view -H "$FINAL" | grep '^@RG' > results/reports/alignment_read_groups.sam
printf 'Alignment finished. Review results/reports/alignment_flagstat.txt and alignment_read_groups.sam\n'
