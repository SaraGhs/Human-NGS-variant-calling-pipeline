# High-Coverage NA12878 — Regional Alignment QC

## Objective
Evaluate alignment integrity, sequencing depth,
and coverage of a high-coverage human genomic region.

## Dataset
- Sample: NA12878 / HG001
- Source: 1000 Genomes high-coverage release
- Original format: CRAM
- Reference: GRCh38/hg38
- Region: chr20:32925000-32955000
- Region length: 30,001 bp
- Software: SAMtools 1.25

## Reference Validation
Chromosome 20 reference MD5:
b18e6c531b0bd70e949a7fc20859cb01

Remote and local reference checksums matched.

## Alignment Results

| Metric | Value |
|---|---:|
| Total alignment records | 7170 |
| Primary records | 7140 |
| Primary mapped records | 7133 |
| Primary mapping rate | 99.90% |
| Properly paired records | 6994 |
| Properly paired percentage | 97.96% |
| Duplicate-marked records | 683 |
| Supplementary alignments | 30 |

## Coverage Assessment

Coverage calculated using minimum mapping
quality 20 and minimum base quality 20.

| Metric | Value |
|---|---:|
| Covered bases | 30001 |
| Coverage breadth | 100% |
| Average depth | 31.1887x |
| Mean base quality | 29.9 |
| Mean mapping quality | 59.7 |

## Interpretation
The regional extraction and alignment QC completed
successfully.

The interval has substantially greater sequencing
depth than the original low-coverage training dataset.

All 30,001 positions have at least one qualifying
base under the coverage assessment settings.

However, average depth and coverage breadth alone
do not establish uniform depth or variant accuracy.

## Limitations
The region includes a documented structural
variation associated with altered read depth.

This interval is suitable for technical analysis
but will not automatically be treated as an
uncomplicated small-variant benchmark region.

## Next Steps
1. Calculate per-position depth statistics.
2. Run regional GATK HaplotypeCaller.
3. Summarize SNPs, indels, DP, GQ and FILTER.
4. Select an appropriate GIAB benchmarking region.
5. Compare calls against matched HG001 truth data.

## Supporting Files
- docs/HIGH_COVERAGE_REFERENCE_CHECK.txt
- docs/HIGH_COVERAGE_DATA_PROVENANCE.md
- results/highcov/chr20_flagstat.txt
- results/highcov/chr20_coverage.txt
- results/highcov/NA12878.chr20_test.bam

## Status
High-coverage regional extraction and QC completed.
High-coverage variant calling pending.
