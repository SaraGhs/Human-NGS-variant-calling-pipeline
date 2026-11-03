#!/usr/bin/env python3
"""Summarize chr19 candidate variants from bcftools-query TSV (no extra dependencies)."""
import argparse
import csv
from collections import Counter
from pathlib import Path
from statistics import median


def parse_num(value):
    try:
        if value in ('', '.', 'NA', None):
            return None
        return float(value)
    except (ValueError, TypeError):
        return None


def variant_class(ref, alt):
    alts = alt.split(',')
    if any(a.startswith('<') or a == '*' for a in alts):
        return 'other'
    if all(len(a) == len(ref) == 1 for a in alts):
        return 'SNP'
    if any(len(a) != len(ref) for a in alts):
        return 'indel'
    return 'other'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', default='results/reports/chr19_exploratory_variants.tsv')
    ap.add_argument('--out', default='docs/CHR19_VARIANT_QC.md')
    ap.add_argument('--summary', default='results/reports/chr19_variant_qc_summary.tsv')
    args = ap.parse_args()
    path = Path(args.input)
    if not path.is_file():
        raise SystemExit(f'Missing input: {path}')
    with path.open(newline='') as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    expected = {'CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'FILTER', 'GT', 'DP', 'GQ'}
    if rows and not expected.issubset(rows[0]):
        raise SystemExit('Unexpected columns in input; required: ' + ', '.join(sorted(expected)))
    if not rows:
        raise SystemExit('No variant rows found: cannot summarize an empty table')
    counts = Counter(variant_class(r['REF'], r['ALT']) for r in rows)
    genotypes = Counter(r['GT'] for r in rows)
    filters = Counter(r['FILTER'] for r in rows)
    depths = [x for r in rows if (x := parse_num(r['DP'])) is not None]
    gqs = [x for r in rows if (x := parse_num(r['GQ'])) is not None]
    quals = [x for r in rows if (x := parse_num(r['QUAL'])) is not None]
    dp_under3 = sum(1 for r in rows if (x := parse_num(r['DP'])) is not None and x < 3)
    dp_under10 = sum(1 for r in rows if (x := parse_num(r['DP'])) is not None and x < 10)
    gq_under20 = sum(1 for r in rows if (x := parse_num(r['GQ'])) is not None and x < 20)
    pass_screen = sum(1 for r in rows if (dp := parse_num(r['DP'])) is not None and (gq := parse_num(r['GQ'])) is not None and dp >= 10 and gq >= 20)
    def stats(values):
        return f'{min(values):g} / {median(values):g} / {max(values):g}' if values else 'not available'
    summary = [
        ('total_candidate_records', len(rows)),
        ('SNPs', counts['SNP']),
        ('indels', counts['indel']),
        ('other', counts['other']),
        ('DP_available', len(depths)),
        ('DP_less_than_3', dp_under3),
        ('DP_less_than_10', dp_under10),
        ('GQ_available', len(gqs)),
        ('GQ_less_than_20', gq_under20),
        ('DP_at_least_10_and_GQ_at_least_20', pass_screen),
        ('FILTER_missing_dot', filters['.']),
    ]
    out = Path(args.out)
    sumout = Path(args.summary)
    out.parent.mkdir(parents=True, exist_ok=True)
    sumout.parent.mkdir(parents=True, exist_ok=True)
    with sumout.open('w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['metric', 'value'])
        w.writerows(summary)
        for k, v in sorted(genotypes.items()):
            w.writerow(['genotype_' + k, v])
        for k, v in sorted(filters.items()):
            w.writerow(['FILTER_' + k, v])
    lines = [
        '# Chromosome 19 — Exploratory Variant Quality Evaluation',
        '',
        '**Dataset:** NA12878 (HG001), public Broad Institute GATK WGS test subset.',
        '',
        '**Method:** Single-sample GATK HaplotypeCaller v4.6.2.0, hg38, chr19 interval. Raw VCF evaluated using BCFtools; tabulated here from the extracted per-sample variant fields.',
        '',
        '## Raw variant records',
        '',
        '| Classification | Count |',
        '|---|---:|',
        f'| SNPs | {counts["SNP"]} |',
        f'| Indels | {counts["indel"]} |',
        f'| Other | {counts["other"]} |',
        f'| **Total raw candidate records** | **{len(rows)}** |',
        '',
        '## Read depth and quality',
        '',
        '| Metric | Value |',
        '|---|---:|',
        f'| Depth (DP) present in records | {len(depths)} |',
        f'| DP min / median / max | {stats(depths)} |',
        f'| Genotype quality (GQ) present | {len(gqs)} |',
        f'| GQ min / median / max | {stats(gqs)} |',
        f'| Site quality (QUAL) min / median / max | {stats(quals)} |',
        f'| Candidates with DP < 3 | {dp_under3} |',
        f'| Candidates with DP < 10 | {dp_under10} |',
        f'| Candidates with GQ < 20 | {gq_under20} |',
        f'| Candidates meeting exploratory DP >= 10 AND GQ >= 20 screen | {pass_screen} |',
        '',
        '## Genotype distribution', '', '| Genotype | Records |', '|---|---:|',
    ]
    lines.extend(f'| `{k}` | {v} |' for k, v in sorted(genotypes.items()))
    lines.extend(['', '## FILTER field distribution', '', '| FILTER value | Records |', '|---|---:|'])
    lines.extend(f'| `{k}` | {v} |' for k, v in sorted(filters.items()))
    lines.extend([
        '',
        '## Interpretation and limitations',
        '',
        '- This dataset has very low chromosome-19 average depth (approximately 0.0025x), so genotype calls cannot be treated as reliable without further evidence.',
        '- `GT=1/1` in a site with DP 1–2 and low GQ is **not** proof of a true homozygous alternate genotype.',
        '- Site QUAL does not replace per-sample depth or GQ checks.',
        '- A `FILTER` value of `.` means filtering status was not set. It does **not** mean the call passed quality filters.',
        '- The DP >= 10 and GQ >= 20 rule above is a **descriptive sensitivity check**, not a validated clinical or research filter.',
        '- Do not treat these 25 raw records as independently confirmed variants or use them to make clinical conclusions.',
        '- Variant precision/recall and biological interpretation require a more suitable high-coverage dataset and an assembly-matched truth set.',
        '',
        '## Supporting files',
        '',
        '- `results/variants/NA12878.chr19.exploratory.raw.vcf.gz` (local only; not intended for GitHub)',
        '- `results/reports/chr19_exploratory_bcftools_stats.txt`',
        '- `results/reports/chr19_exploratory_variants.tsv`',
        '- `results/reports/chr19_variant_qc_summary.tsv`',
        '- `logs/gatk_haplotypecaller_chr19.log`',
        '',
        '**Status:** Exploratory variant-level summary complete; benchmark validation not performed.',
    ])
    out.write_text('\n'.join(lines) + '\n')
    print(f'Wrote {out}')
    print(f'Wrote {sumout}')
    for k, v in summary:
        print(f'{k}: {v}')

if __name__ == '__main__':
    main()
