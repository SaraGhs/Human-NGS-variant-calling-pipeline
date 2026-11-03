# Chromosome 19 — Exploratory Variant Quality Evaluation

**Dataset:** NA12878 (HG001), public Broad Institute GATK WGS test subset.

**Method:** Single-sample GATK HaplotypeCaller v4.6.2.0, hg38, chr19 interval. Raw VCF evaluated using BCFtools; tabulated here from the extracted per-sample variant fields.

## Raw variant records

| Classification | Count |
|---|---:|
| SNPs | 22 |
| Indels | 3 |
| Other | 0 |
| **Total raw candidate records** | **25** |

## Read depth and quality

| Metric | Value |
|---|---:|
| Depth (DP) present in records | 25 |
| DP min / median / max | 1 / 2 / 2 |
| Genotype quality (GQ) present | 25 |
| GQ min / median / max | 3 / 6 / 6 |
| Site quality (QUAL) min / median / max | 30.48 / 37.32 / 78.32 |
| Candidates with DP < 3 | 25 |
| Candidates with DP < 10 | 25 |
| Candidates with GQ < 20 | 25 |
| Candidates meeting exploratory DP >= 10 AND GQ >= 20 screen | 0 |

## Genotype distribution

| Genotype | Records |
|---|---:|
| `1/1` | 25 |

## FILTER field distribution

| FILTER value | Records |
|---|---:|
| `.` | 25 |

## Interpretation and limitations

- This dataset has very low chromosome-19 average depth (approximately 0.0025x), so genotype calls cannot be treated as reliable without further evidence.
- `GT=1/1` in a site with DP 1–2 and low GQ is **not** proof of a true homozygous alternate genotype.
- Site QUAL does not replace per-sample depth or GQ checks.
- A `FILTER` value of `.` means filtering status was not set. It does **not** mean the call passed quality filters.
- The DP >= 10 and GQ >= 20 rule above is a **descriptive sensitivity check**, not a validated clinical or research filter.
- Do not treat these 25 raw records as independently confirmed variants or use them to make clinical conclusions.
- Variant precision/recall and biological interpretation require a more suitable high-coverage dataset and an assembly-matched truth set.

## Supporting files

- `results/variants/NA12878.chr19.exploratory.raw.vcf.gz` (local only; not intended for GitHub)
- `results/reports/chr19_exploratory_bcftools_stats.txt`
- `results/reports/chr19_exploratory_variants.tsv`
- `results/reports/chr19_variant_qc_summary.tsv`
- `logs/gatk_haplotypecaller_chr19.log`

**Status:** Exploratory variant-level summary complete; benchmark validation not performed.
