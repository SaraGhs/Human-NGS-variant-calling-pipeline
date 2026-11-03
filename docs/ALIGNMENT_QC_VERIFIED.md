# Alignment Quality Control — NA12878

## Methods
Trimmed paired-end reads were aligned to an hg38
reference using BWA-MEM. BAM processing included
sorting, merging, mate correction, duplicate marking,
and indexing using SAMtools.

Alignment statistics were evaluated using samtools
flagstat, view, and coverage.

## Alignment Results
- Input read pairs: 23,848
- Primary alignment records: 47,696
- Primary mapped reads: 47,659 (99.92%)
- Mapped primary reads with MAPQ >=20: 44,470
- Fraction of mapped reads with MAPQ >=20: 93.31%
- Properly paired records: 46,770 (98.06%)
- Primary duplicate records: 11
- Supplementary alignment records: 478
- Read groups: 3

## Chromosome 1 Coverage
- Reference length: 248,956,422 bp
- Covered bases: 644,601
- Coverage breadth: 0.258921%
- Mean depth: 0.00284488x
- Mean base quality: 33.3
- Mean mapping quality: 58.8

## Chromosome 19 Coverage
- Reference length: 58,617,616 bp
- Covered bases: 135,613
- Coverage breadth: 0.231352%
- Mean depth: 0.00249563x
- Mean base quality: 33.1
- Mean mapping quality: 59.7

## Interpretation
Reference-guided alignment was successful, with
99.92% of primary reads mapped and 93.31% of mapped
primary reads having MAPQ >=20.

However, extremely low genomic coverage prevents
reliable genome-wide germline variant interpretation.

This dataset is suitable for pipeline validation
and exploratory variant-calling exercises.

## Remaining Work
- Evaluate cross-chromosome mate alignments
- Perform exploratory GATK HaplotypeCaller analysis
- Inspect and summarize VCF output
- Investigate a higher-coverage benchmarking dataset
- Finalize reference provenance and reproducibility

## Supporting Outputs
- results/reports/alignment_flagstat.txt
- results/reports/alignment_read_groups.sam
- results/reports/alignment_idxstats.tsv
- results/reports/alignment_additional_qc.txt

## Status
Alignment QC complete.
Exploratory variant calling pending.
