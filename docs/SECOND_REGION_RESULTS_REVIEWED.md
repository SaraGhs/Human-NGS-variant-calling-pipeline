# Second regional benchmark — reviewed results

## Scope and provenance

- Sample: NA12878 / HG001, 1000 Genomes high-coverage GRCh38-aligned CRAM.
- Candidate method: GATK HaplotypeCaller on a preselected 200,000 bp chromosome 20 region (`chr20:20000001-20200000`).
- Truth source: NIST Genome in a Bottle HG001 v4.2.1 GRCh38 small-variant VCF and high-confidence BED.
- Comparison: Illumina hap.py; `ALL` rows of SNP and INDEL summary; confident-region restricted scoring.
- GIAB confidence overlap: **97.31%** of this 200 kb interval (percent rounded from the user's local output).
- Regional mean depth: **31.8787×**, coverage breadth **99.9995%** under prior SAMtools coverage criteria.

## Variant results

The original candidate VCF was summarized with BCFtools as **428 variant records**, with 349 SNP and 80 indel classifications and 10 multiallelic sites. Summing these category counts does not always recover the number of VCF records, because variant types and multiallelic alleles can be counted differently.

| hap.py `ALL` measurement | SNP | INDEL | Sum |
|---|---:|---:|---:|
| Truth units (`TRUTH.TOTAL`) | 331 | 53 | 384 |
| True-positive truth units (`TRUTH.TP`) | 331 | 53 | 384 |
| False negatives (`TRUTH.FN`) | 0 | 0 | 0 |
| Query false positives (`QUERY.FP`) | 0 | 0 | 0 |
| Unassessed query units (`QUERY.UNK`) | 18 | 31 | 49 |
| Total query units (`QUERY.TOTAL`) | 349 | 85 | 434 |
| Precision | 1.00 | 1.00 | 1.00 (within assessed calls) |
| Recall | 1.00 | 1.00 | 1.00 (within assessed calls) |
| F1 | 1.00 | 1.00 | 1.00 (within assessed calls) |

**Interpretation:** All 384 truth-side units considered in the benchmark were recovered, and no false-positive query units were reported within GIAB-evaluable regions. However, **49 of 434 query benchmark units (11.29%) were unassessed**, not verified as true positives. Haplotype-aware normalization and representation can make hap.py units differ from raw BCFtools record counts; do not present those totals as directly interchangeable.

## Comparison with prior 30 kb interval

In the first test interval, GIAB confidence BED covered 76.29% of the 30,001-bp target. The benchmark reported 21 truth-side TP units, zero false positives and zero false negatives, and seven unassessed query units, with precision, recall and F1 of 1.00 within assessed regions.

The second preselected interval was ~6.7 times larger and had higher GIAB confident-region breadth (97.31%). Both intervals used the same biological sample and do not constitute independent biological replication.

## Limitations and follow-up

- The benchmark was performed on a single sample and two selected genomic intervals, so perfect regional metrics do not establish genome-wide accuracy.
- Raw variant-call FILTER status should be documented separately; `PASS` in the benchmarking summary is not equivalent to a validated filtering protocol for the raw GATK VCF.
- Calls outside the confident BED are not independently verified. Report `QUERY.UNK`, rather than treating these as errors or correct variants.
- Small-variant benchmarking does not validate structural variants or long complex genomic regions.
- Review the full hap.py log and preserve software versions, reference checksum evidence, input provenance, and the exact command before publishing.

## Supporting local files

- `results/highcov/chr20_validation_stats.txt`
- `results/highcov/chr20_validation_coverage.txt`
- `results/benchmark/NA12878_chr20_validation_happy.summary.csv`
- `results/benchmark/NA12878_chr20_happy.summary.csv`
- `results/benchmark/chr20_independent_confident.bed`
- `logs/happy_chr20_validation.log`

**Status:** Second-region hap.py assessment complete. Final reproducibility audit and GitHub packaging pending. This report describes results supplied by the user; the automated script in this package should be run locally to regenerate the comparative table from the original outputs.
