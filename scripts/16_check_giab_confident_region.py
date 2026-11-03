#!/usr/bin/env python3
"""Check fraction of a 1-based inclusive interval inside HG001 GIAB BED confidence regions."""
import argparse
from pathlib import Path

p=argparse.ArgumentParser()
p.add_argument('--bed',default='data/giab/HG001_GRCh38_1_22_v4.2.1_benchmark.bed')
p.add_argument('--chrom',default='chr20')
p.add_argument('--start',type=int,default=32925000)
p.add_argument('--end',type=int,default=32955000)
a=p.parse_args()
if a.end<a.start or a.start<1: raise SystemExit('Invalid region')
path=Path(a.bed)
if not path.exists(): raise SystemExit(f'Missing BED: {path}')
start,end=a.start-1,a.end # 1-based inclusive to BED half-open
intervals=[]
with path.open() as fh:
 for line in fh:
  if not line.strip() or line.startswith(('#','track','browser')): continue
  fields=line.strip().split('\t')
  if len(fields)<3 or fields[0] not in {a.chrom,a.chrom.removeprefix('chr')}: continue
  try: s,e=int(fields[1]),int(fields[2])
  except ValueError: continue
  left,right=max(s,start),min(e,end)
  if right>left: intervals.append((left,right))
intervals.sort()
merged=[]
for s,e in intervals:
 if merged and s<=merged[-1][1]: merged[-1]=(merged[-1][0],max(merged[-1][1],e))
 else: merged.append((s,e))
bases=sum(e-s for s,e in merged)
total=end-start
out=Path('docs/CHR20_GIAB_CONFIDENCE_CHECK.md')
out.parent.mkdir(parents=True,exist_ok=True)
report=f'''# HG001 GIAB confident-region overlap\n\n- BED file: `{path}`\n- Candidate interval: `{a.chrom}:{a.start}-{a.end}` (1-based, inclusive; {total:,} bases).\n- Overlap with GIAB v4.2.1 high-confidence BED: **{bases:,} bases ({100*bases/total:.2f}%)**.\n- Number of intersecting merged confident segments: **{len(merged)}**.\n\n**Warning:** High-confidence BED overlap does not by itself validate any calls. Variants outside the confidence BED cannot be classified reliably as false positives for a standard benchmark. Any formal comparison must use a matched HG001 GRCh38 truth VCF and a haplotype-aware comparator (e.g., hap.py).\n\n## Overlapping BED coordinates (0-based, half-open)\n\n'''
for s,e in merged: report+=f'- `{a.chrom}:{s}-{e}` ({e-s:,} bases)\n'
if not merged: report+='- None. Select another region for benchmarking.\n'
out.write_text(report)
print(f'High-confidence overlap: {bases}/{total} bases ({100*bases/total:.2f}%)')
print(f'Wrote {out}')
