# High-Coverage Variant Quality Evaluation

## Dataset
- Sample: NA12878 (HG001)
- Source: 1000 Genomes high-coverage WGS
- Reference: GRCh38/hg38
- Interval: chr20:32925000-32955000
- Caller: GATK HaplotypeCaller

## Regional Coverage
- Positions assessed: 30,001
- Mean depth: 30.99x
- Positions with >=10x: 98.71%
- Positions with >=20x: 87.91%
- Positions with >=30x: 63.56%

## Variant Results
- Raw candidates: 28
- SNPs: 18
- Indels: 10
- DP min/median/max: 7/30/42
- GQ min/median/max: 21/99/99
- QUAL min/median/max: 141.6/390.14/797.6
- DP >=10 and GQ >=20: 26 of 28
- FILTER unset: 28
- Explicit FILTER PASS: 0

## GIAB Confidence Assessment
- Total interval: 30,001 bp
- Overlapping high-confidence bases: 22,888
- High-confidence overlap: 76.29%
- Intersecting confidence segments: 8

## Interpretation
The higher-coverage dataset provides substantially
stronger genotype support than the original low-depth
demonstration.

The 26 calls meeting exploratory DP/GQ thresholds
are not yet validated true-positive variants.

All 28 calls have an unset FILTER status.

Formal variant accuracy assessment requires comparison
against the matching HG001 benchmark VCF using a
haplotype-aware comparison method.

## Limitations
The test interval overlaps a region associated with
structural variation and is not wholly covered by
the GIAB high-confidence benchmark.

This small interval cannot establish genome-wide
variant-calling performance.

## Supporting Evidence
- results/highcov/chr20_depth_thresholds.txt
- results/highcov/chr20_variant_metrics.tsv
- docs/CHR20_GIAB_CONFIDENCE_CHECK.md
- docs/HIGH_COVERAGE_VARIANT_REVIEW.md

## Status
High-coverage candidate variant evaluation completed.
GIAB truth comparison pending.
