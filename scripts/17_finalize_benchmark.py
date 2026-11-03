#!/usr/bin/env python3
"""Summarize hap.py benchmarking output without conflating unassessed and false-positive calls.

Reads the hap.py summary and local GIAB confident-region BED, checks consistency,
then generates a reproducible review report. Requires only Python 3 standard library.
"""
import argparse
import csv
from decimal import Decimal, InvalidOperation
from pathlib import Path


def as_integer(row, column):
    raw = row.get(column, '')
    if raw in ('', None):
        raise ValueError(f'Missing {column} for {row.get("Type")}')
    try:
        number = Decimal(raw)
    except InvalidOperation:
        raise ValueError(f'Invalid {column}={raw}') from None
    if number != number.to_integral_value() or number < 0:
        raise ValueError(f'Expected nonnegative integer {column}={raw}')
    return int(number)


def read_summary(path):
    with open(path, newline='', encoding='utf-8-sig') as stream:
        reader = csv.DictReader(stream)
        required = {'Type','Filter','TRUTH.TOTAL','TRUTH.TP','TRUTH.FN','QUERY.TOTAL',
                    'QUERY.FP','QUERY.UNK','METRIC.Recall','METRIC.Precision','METRIC.F1_Score'}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f'Missing required columns: {sorted(required-set(reader.fieldnames or []))}')
        rows = [r for r in reader if r.get('Filter') == 'ALL' and r.get('Type') in ('SNP','INDEL')]
    if len(rows) != 2 or {r['Type'] for r in rows} != {'SNP','INDEL'}:
        raise ValueError('Expected exactly one ALL row for SNP and INDEL')
    result = {r['Type']:r for r in rows}
    for kind, row in result.items():
        if as_integer(row,'TRUTH.TOTAL') != as_integer(row,'TRUTH.TP') + as_integer(row,'TRUTH.FN'):
            raise ValueError(f'Inconsistent truth counts: {kind}')
        if as_integer(row,'QUERY.FP')+as_integer(row,'QUERY.UNK') > as_integer(row,'QUERY.TOTAL'):
            raise ValueError(f'Inconsistent query counts: {kind}')
    return result


def overlap_bases(path, chrom, start1, end1):
    # BED is 0-based half-open; start1/end1 are 1-based inclusive.
    start0 = start1 - 1
    end0 = end1
    parts = []
    with open(path, encoding='utf-8') as stream:
        for line in stream:
            if not line.strip() or line.startswith(('#','track','browser')):
                continue
            fields=line.rstrip().split('\t')
            if len(fields)<3 or fields[0]!=chrom:
                continue
            a=max(int(fields[1]), start0)
            b=min(int(fields[2]), end0)
            if a<b: parts.append((a,b))
    parts.sort()
    merged=[]
    for a,b in parts:
        if merged and a<=merged[-1][1]:
            merged[-1]=(merged[-1][0], max(merged[-1][1],b))
        else:
            merged.append((a,b))
    return sum(b-a for a,b in merged),len(merged)


def build(summary, region_bases, segments, interval_bases):
    kinds=('SNP','INDEL')
    total=lambda key: sum(as_integer(summary[k],key) for k in kinds)
    truth=total('TRUTH.TOTAL')
    tp=total('TRUTH.TP')
    fn=total('TRUTH.FN')
    query=total('QUERY.TOTAL')
    fp=total('QUERY.FP')
    unk=total('QUERY.UNK')
    lines=[
        '# HG001/NA12878 chr20 — GIAB Benchmark Review',
        '',
        '## Scope and provenance',
        '- Query: GATK HaplotypeCaller raw candidates, `results/highcov/NA12878.chr20.highcov.raw.vcf.gz`.',
        '- Truth: GIAB HG001 GRCh38 v4.2.1 small-variant benchmark.',
        '- Comparator: hap.py; evaluation region `chr20:32925000-32955000` (1-based inclusive).',
        f'- Evaluable confident bases: **{region_bases:,} / {interval_bases:,} ({region_bases/interval_bases:.2%})**, across {segments} merged segments.',
        '- Original candidate VCF has unset `FILTER=.`; no validated quality filters were applied. The hap.py `PASS` summary row is not proof that these input calls passed a scientific quality filter.',
        '',
        '## hap.py results — ALL rows',
        '',
        '| Variant class | GIAB truth | Truth TP | Truth FN | Query calls | Query FP | Query UNK | Recall | Precision | F1 |',
        '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|',
    ]
    for kind in kinds:
        r=summary[kind]
        lines.append('| '+ ' | '.join([
            kind,str(as_integer(r,'TRUTH.TOTAL')),str(as_integer(r,'TRUTH.TP')),
            str(as_integer(r,'TRUTH.FN')),str(as_integer(r,'QUERY.TOTAL')),
            str(as_integer(r,'QUERY.FP')),str(as_integer(r,'QUERY.UNK')),
            f'{float(r["METRIC.Recall"]):.2%}',f'{float(r["METRIC.Precision"]):.2%}',
            f'{float(r["METRIC.F1_Score"]):.4f}'])+' |')
    lines.extend([
        f'| **Combined counts** | **{truth}** | **{tp}** | **{fn}** | **{query}** | **{fp}** | **{unk}** | — | — | — |',
        '',
        '## Interpretation',
        f'- The benchmark matched {tp}/{truth} GIAB truth-side variants in the evaluable region; {fn} false negatives and {fp} query false positives were reported.',
        f'- **{unk}/{query} query calls ({unk/query:.1%})** were **unassessed** (`QUERY.UNK`), not classified as true or false positives.',
        '- Precision, recall, and F1 are 100%/100%/1.00 **only within this small evaluated interval and confidence mask**; do not generalize to the genome.',
        '- The selected interval intersects a known structural-variation example; performance can differ markedly elsewhere and across challenging genomic contexts.',
        '- These metrics describe unfiltered exploratory variant calls. A quality-filtering strategy was not independently validated.',
        '- Differences between truth and query VCF representation may yield different truth-side vs query-side true-positive counts; use hap.py output rather than naive coordinate matches.',
        '',
        '## Reproducibility and next work',
        '- Preserve `results/benchmark/NA12878_chr20_happy.summary.csv`, `.extended.csv`, and `.metrics.json` as source results (if present).',
        '- Inspect annotated hap.py comparison VCF for TP/FP/UNK annotations and confirm sample/reference identity.',
        '- Record exact commands, software versions, and SHA-256 checksums of input VCF and BED files.',
        '- Extend the analysis to a larger, separate, well-covered and uncomplicated high-confidence region before making broader performance claims.',
        '- Publish scripts, small summaries and figures only after the full scientific and reproducibility review; exclude large sequencing data.',
        '',
        '**Status:** Small-interval hap.py benchmark completed; broader validation and final project audit remain pending.',
    ])
    return '\n'.join(lines)+'\n',[(k, summary[k]) for k in kinds]


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--summary', default='results/benchmark/NA12878_chr20_happy.summary.csv')
    ap.add_argument('--confident-bed',default='results/benchmark/chr20_giab_confident.bed')
    ap.add_argument('--output',default='docs/GIAB_BENCHMARK_RESULTS_REVIEWED.md')
    ap.add_argument('--metrics',default='results/benchmark/chr20_benchmark_review_metrics.tsv')
    args=ap.parse_args()
    summary=read_summary(args.summary)
    bases,segments=overlap_bases(args.confident_bed,'chr20',32925000,32955000)
    if bases==0: raise ValueError('Zero high-confidence bases: check BED format and chromosome names')
    report,rows=build(summary,bases,segments,30001)
    Path(args.output).parent.mkdir(parents=True,exist_ok=True)
    Path(args.output).write_text(report,encoding='utf-8')
    Path(args.metrics).parent.mkdir(parents=True,exist_ok=True)
    with open(args.metrics,'w',newline='',encoding='utf-8') as stream:
        writer=csv.writer(stream,delimiter='\t')
        columns=['TRUTH.TOTAL','TRUTH.TP','TRUTH.FN','QUERY.TOTAL','QUERY.FP','QUERY.UNK','METRIC.Recall','METRIC.Precision','METRIC.F1_Score']
        writer.writerow(['Type']+columns)
        for kind,row in rows: writer.writerow([kind]+[row.get(c,'') for c in columns])
    print(f'Wrote {args.output}; {args.metrics}; confident bases {bases}, segments {segments}')

if __name__=='__main__':
    main()
