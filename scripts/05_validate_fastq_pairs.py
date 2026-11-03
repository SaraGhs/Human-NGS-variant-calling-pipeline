#!/usr/bin/env python3
"""Read-only record-by-record validation of the 3 NA12878 paired-end FASTQ groups.
Writes a small per-group TSV; does not modify sequencing data.
"""
from __future__ import annotations
import argparse, csv, gzip
from pathlib import Path

def op(path):
    return gzip.open(path,'rt',encoding='ascii') if str(path).endswith('.gz') else path.open('rt',encoding='ascii')

def record_iter(path):
    with op(path) as handle:
        n=0
        while True:
            header=handle.readline()
            if not header: return
            seq=handle.readline(); plus=handle.readline(); qual=handle.readline()
            n+=1
            if not seq or not plus or not qual:
                raise ValueError(f'Truncated FASTQ record {n} in {path}')
            if not header.startswith('@') or not plus.startswith('+'):
                raise ValueError(f'Invalid FASTQ header/separator record {n} in {path}')
            s,q=seq.rstrip('\r\n'),qual.rstrip('\r\n')
            if len(s)!=len(q) or not s:
                raise ValueError(f'Sequence/quality length mismatch or empty read at record {n} in {path}')
            ident=header[1:].split()[0].removesuffix('/1').removesuffix('/2')
            yield ident,len(s)

def examine_pair(r1,r2):
    reads=0; minlen=None; maxlen=0; total=0
    from itertools import zip_longest
    for rec1,rec2 in zip_longest(record_iter(r1),record_iter(r2)):
        reads+=1
        if rec1 is None or rec2 is None:
            raise ValueError(f'Unequal mate counts for {r1.name} and {r2.name} at record {reads}')
        if rec1[0] != rec2[0]:
            raise ValueError(f'Mismatched read IDs in pair {r1.name} and {r2.name} at record {reads}: {rec1[0]} != {rec2[0]}')
        for l in (rec1[1],rec2[1]):
            total+=l; minlen=l if minlen is None else min(minlen,l); maxlen=max(maxlen,l)
    if reads==0: raise ValueError(f'Empty FASTQ pair: {r1.name}, {r2.name}')
    return reads,minlen,maxlen,round(total/(2*reads),2)

def inventory(folder):
    r1s=sorted(folder.glob('*_reads_1.fastq'))+sorted(folder.glob('*_reads_1.fastq.gz'))
    r2s=sorted(folder.glob('*_reads_2.fastq'))+sorted(folder.glob('*_reads_2.fastq.gz'))
    if len(r1s)!=3 or len(r2s)!=3:
        raise ValueError(f'Expected three R1 and three R2 files; found {len(r1s)} and {len(r2s)}')
    rows=[]
    for r1 in r1s:
        group=r1.name.split('_reads_1.fastq')[0]
        r2=folder / r1.name.replace('_reads_1.fastq','_reads_2.fastq')
        if not r2.is_file(): raise ValueError(f'Missing paired FASTQ for {r1.name}: {r2}')
        count,mn,mx,mean=examine_pair(r1,r2)
        rows.append(dict(group=group,r1=r1.name,r2=r2.name,read_pairs=count,min_length_bp=mn,max_length_bp=mx,mean_length_bp=mean,status='PASS'))
    return rows

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--fastq-dir',type=Path,default=Path('data/source/NA12878_20k'))
    ap.add_argument('--out',type=Path,default=Path('docs/fastq_pair_inventory.tsv'))
    a=ap.parse_args()
    rows=inventory(a.fastq_dir)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    with a.out.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t'); w.writeheader(); w.writerows(rows)
    print(f'Validated {len(rows)} paired sequencing groups; output: {a.out}')
    for r in rows: print(f"  {r['group']}: {r['read_pairs']} read pairs; read lengths {r['min_length_bp']}-{r['max_length_bp']} bp")
if __name__=='__main__': main()
