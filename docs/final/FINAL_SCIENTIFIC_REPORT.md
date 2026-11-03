# Human NGS Variant Calling — Scientific Report

## Project Overview

This project analyzed publicly available NA12878/HG001
sequencing data using two complementary workflows:

1. Raw FASTQ-to-VCF analysis using a small test dataset.
2. High-coverage regional variant calling and GIAB benchmarking.

## Methods

- Raw sequencing QC: FastQC and MultiQC
- Read preprocessing: Trimmomatic
- Reference alignment: BWA-MEM
- BAM processing: SAMtools
- Germline variant discovery: GATK HaplotypeCaller
- Variant analysis: BCFtools
- Benchmarking: hap.py
- Reference: GRCh38/hg38
- Truth dataset: GIAB HG001 v4.2.1

## Preprocessing

Original read pairs: 30,000
Retained paired reads: 23,848
Paired-read retention: 79.49%

## Low-Coverage Experiment

Chromosome 19:
- Raw candidates: 25
- SNPs: 22
- Indels: 3
- Read depths: 1–2

These calls were not considered reliable genotypes.

## High-Coverage Regional Benchmarking

| Region | Type | Truth total | TP | FN | FP | Unassessed | Precision | Recall |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| 30 kb | SNP | 15 | 15 | 0 | 0 | 3 | 1.00 | 1.00 |
| 30 kb | INDEL | 6 | 6 | 0 | 0 | 4 | 1.00 | 1.00 |
| 200 kb | SNP | 331 | 331 | 0 | 0 | 18 | 1.00 | 1.00 |
| 200 kb | INDEL | 53 | 53 | 0 | 0 | 31 | 1.00 | 1.00 |

### GIAB High-Confidence Regions

30 kb interval: 22,888 confident bases.

200 kb interval: 194,617 confident bases.

### Interpretation

Both regional benchmarks showed high agreement within
the GIAB-evaluable regions.

The 30 kb analysis contained 7 unassessed
query units, while the 200 kb analysis contained
49 unassessed query units.

These are not confirmed true or false positives.

## Limitations

- Both regional analyses used the same biological sample.
- Results are not representative of the entire genome.
- Original variant FILTER values were unset.
- The 30 kb region overlaps structural complexity.
- Reference provenance requires final verification.
- No clinical conclusions can be drawn.

## Supporting Documentation

See the existing preprocessing, alignment, variant
quality, and GIAB benchmarking reports in docs/.

## Status

Scientific analyses complete for the defined project.
Reproducibility and publication audit pending.
