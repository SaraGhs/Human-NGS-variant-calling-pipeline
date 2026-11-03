#!/usr/bin/env python3
"""Summarize research-only NA12878 chr20 variant calls from the local TSV."""
import csv
import math
import statistics
from collections import Counter
from pathlib import Path

ROOT=Path(__file__).resolve().parent.parent
SRC=ROOT/'results/highcov/chr20_highcov_variants.tsv'
OUT=ROOT/'docs/HIGH_COVERAGE_VARIANT_REVIEW.md'
TABLE=ROOT/'results/highcov/chr20_variant_metrics.tsv'

def number(s):
    try:
        x=float(s)
        return x if math.isfinite(x) else None
    except (ValueError,TypeError):
        return None

def describe(vals):
    if not vals: return 'not available'
    return f'{min(vals):g} / {statistics.median(vals):g} / {max(vals):g}'

def classify(ref, alt):
    if ',' in alt: return 'multiallelic'
    if len(ref)==1 and len(alt)==1: return 'SNP'
    if len(ref)!=len(alt): return 'indel'
    return 'other'

def run():
    if not SRC.is_file(): raise SystemExit(f'Missing input: {SRC}')
    with SRC.open(newline='') as f:
        records=list(csv.DictReader(f,delimiter='\t'))
    if not records: raise SystemExit('No variant records in TSV. Do not generate conclusions.')
    names={'CHROM','POS','REF','ALT','QUAL','FILTER','GT','DP','GQ'}
    if not names.issubset(records[0]): raise SystemExit('Missing expected fields in variant TSV.')
    n=len(records)
    counts=Counter(classify(r['REF'],r['ALT']) for r in records)
    gt=Counter(r['GT'] for r in records)
    filt=Counter(r['FILTER'] for r in records)
    dps=[number(r['DP']) for r in records]
    gqs=[number(r['GQ']) for r in records]
    quals=[number(r['QUAL']) for r in records]
    screened=sum(d is not None and q is not None and d>=10 and q>=20 for d,q in zip(dps,gqs))
    metric={
        'raw_candidate_records':n,
        'snps':counts['SNP'],
        'indels':counts['indel'],
        'other_or_multiallelic':n-counts['SNP']-counts['indel'],
        'dp_min_median_max':describe([v for v in dps if v is not None]),
        'gq_min_median_max':describe([v for v in gqs if v is not None]),
        'qual_min_median_max':describe([v for v in quals if v is not None]),
        'records_dp_ge_10_and_gq_ge_20':screened,
        'records_without_DP':sum(v is None for v in dps),
        'records_without_GQ':sum(v is None for v in gqs),
        'filter_unset_dot':filt['.'],
        'filter_pass':filt['PASS']
    }
    TABLE.parent.mkdir(parents=True,exist_ok=True)
    with TABLE.open('w') as f:
        f.write('metric\tvalue\n')
        for k,v in metric.items(): f.write(f'{k}\t{v}\n')
    md=f'''# High-coverage chromosome 20 variant review — NA12878 (HG001)\n\n## Methods\n\nSingle-sample GATK HaplotypeCaller v4.6.2.0 on a 30,001-bp chromosome 20 interval (`chr20:32925000-32955000`) using the verified GRCh38-compatible reference and indexed high-coverage BAM. Variant-level measurements were exported from the VCF using BCFtools. This document is generated from the local per-record TSV.\n\n## Independent QC evidence\n\n- Interval positions assessed: **30,001**.\n- Mean overlap-aware depth (base Q >=20, MAPQ >=20): **30.99×**.\n- Depth >=10×: **98.71%** of positions.\n- Depth >=20×: **87.91%** of positions.\n- Depth >=30×: **63.56%** of positions.\n- Source: `results/highcov/chr20_depth_thresholds.txt`.\n\n## Variant observations\n\n- Raw candidate records: **{n}**.\n- SNP records: **{counts['SNP']}**.\n- Indel records: **{counts['indel']}**.\n- Other or multiallelic records: **{n-counts['SNP']-counts['indel']}**.\n- DP min / median / max: **{metric['dp_min_median_max']}**.\n- GQ min / median / max: **{metric['gq_min_median_max']}**.\n- QUAL min / median / max: **{metric['qual_min_median_max']}**.\n- DP >=10 and GQ >=20 (illustrative screen only): **{screened}/{n}**.\n- FILTER field unset (`.`): **{filt['.']}**; FILTER explicitly `PASS`: **{filt['PASS']}**.\n\n## Genotype distribution\n\n'''
    for k,v in sorted(gt.items()): md+=f'- `{k}`: {v}\n'
    md+='''\n## Interpretation and limitations\n\nA higher average depth and stronger sample-level quality evidence than the preceding low-coverage test support evaluating raw regional calls. However, **the candidates are not yet validated**. An unset FILTER (`.`) is not equivalent to `PASS`. Read trimming, selection of a 30 kb interval, region-only CRAM extraction, and an interval containing known structural variation may affect the callset, including edge positions and complex loci. The illustrative DP/GQ threshold is not a validated filter. Formal benchmarking requires the matched HG001/NA12878 GIAB GRCh38 VCF, its high-confidence BED, consistent reference alleles and variant representation, and haplotype-aware comparison. Do not report sensitivity, precision or accuracy until comparison is run.\n\n## Evidence\n\n- `results/highcov/chr20_coverage.txt`\n- `results/highcov/chr20_depth_thresholds.txt`\n- `results/highcov/chr20_highcov_bcftools_stats.txt`\n- `results/highcov/chr20_highcov_variants.tsv`\n- `results/highcov/chr20_variant_metrics.tsv`\n\n**Status:** Preliminary high-coverage variant review completed; GIAB benchmarking pending.\n'''
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.write_text(md)
    print(f'Wrote {OUT} and {TABLE}')
    print(f'{n} records, {counts["SNP"]} SNPs, {counts["indel"]} indels; illustrative screen {screened}/{n}')

if __name__=='__main__': run()
