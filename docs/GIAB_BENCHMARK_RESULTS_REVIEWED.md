# HG001/NA12878 chr20 — GIAB Benchmark Review

## Scope and provenance
- Query: GATK HaplotypeCaller raw candidates, `results/highcov/NA12878.chr20.highcov.raw.vcf.gz`.
- Truth: GIAB HG001 GRCh38 v4.2.1 small-variant benchmark.
- Comparator: hap.py; evaluation region `chr20:32925000-32955000` (1-based inclusive).
- Evaluable confident bases: **22,888 / 30,001 (76.29%)**, across 8 merged segments.
- Original candidate VCF has unset `FILTER=.`; no validated quality filters were applied. The hap.py `PASS` summary row is not proof that these input calls passed a scientific quality filter.

## hap.py results — ALL rows

| Variant class | GIAB truth | Truth TP | Truth FN | Query calls | Query FP | Query UNK | Recall | Precision | F1 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| SNP | 15 | 15 | 0 | 18 | 0 | 3 | 100.00% | 100.00% | 1.0000 |
| INDEL | 6 | 6 | 0 | 10 | 0 | 4 | 100.00% | 100.00% | 1.0000 |
| **Combined counts** | **21** | **21** | **0** | **28** | **0** | **7** | — | — | — |

## Interpretation
- The benchmark matched 21/21 GIAB truth-side variants in the evaluable region; 0 false negatives and 0 query false positives were reported.
- **7/28 query calls (25.0%)** were **unassessed** (`QUERY.UNK`), not classified as true or false positives.
- Precision, recall, and F1 are 100%/100%/1.00 **only within this small evaluated interval and confidence mask**; do not generalize to the genome.
- The selected interval intersects a known structural-variation example; performance can differ markedly elsewhere and across challenging genomic contexts.
- These metrics describe unfiltered exploratory variant calls. A quality-filtering strategy was not independently validated.
- Differences between truth and query VCF representation may yield different truth-side vs query-side true-positive counts; use hap.py output rather than naive coordinate matches.

## Reproducibility and next work
- Preserve `results/benchmark/NA12878_chr20_happy.summary.csv`, `.extended.csv`, and `.metrics.json` as source results (if present).
- Inspect annotated hap.py comparison VCF for TP/FP/UNK annotations and confirm sample/reference identity.
- Record exact commands, software versions, and SHA-256 checksums of input VCF and BED files.
- Extend the analysis to a larger, separate, well-covered and uncomplicated high-confidence region before making broader performance claims.
- Publish scripts, small summaries and figures only after the full scientific and reproducibility review; exclude large sequencing data.

**Status:** Small-interval hap.py benchmark completed; broader validation and final project audit remain pending.
