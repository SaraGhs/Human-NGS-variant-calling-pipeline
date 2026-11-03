#!/usr/bin/env python3
"""Export exact raw FastQC metrics from original FastQC ZIP files.

Requires only the Python standard library. No raw reads are copied.
"""
import argparse
import csv
from pathlib import Path
import zipfile

CHECKS = {
    'per_base_quality': 'Per base sequence quality',
    'per_sequence_quality': 'Per sequence quality scores',
    'per_base_n_content': 'Per base N content',
    'adapter_content': 'Adapter Content',
    'overrepresented_sequences': 'Overrepresented sequences',
    'sequence_duplication_levels': 'Sequence Duplication Levels',
}
FIELDS = [
    'fastqc_zip','filename','total_sequences','sequence_length','gc_percent',
    'deduplicated_percent','max_per_base_n_percent','max_n_position',
    *CHECKS.keys()
]


def parse_fastqc_text(text):
    modules = {}
    module = None
    for line in text.splitlines():
        if line.startswith('>>END_MODULE'):
            module = None
            continue
        if line.startswith('>>'):
            parts = line[2:].split('\t')
            module = parts[0]
            modules[module] = {'status': parts[1] if len(parts)>1 else '', 'lines': []}
            continue
        if module is not None:
            modules[module]['lines'].append(line)
    basic = {}
    for line in modules.get('Basic Statistics', {}).get('lines', []):
        if '\t' in line and not line.startswith('#'):
            key,value = line.split('\t', 1)
            basic[key] = value
    dedup = ''
    for line in modules.get('Sequence Duplication Levels',{}).get('lines',[]):
        if line.startswith('#Total Deduplicated Percentage\t'):
            dedup = line.split('\t',1)[1]
            break
    max_n, max_pos = None, ''
    for line in modules.get('Per base N content',{}).get('lines',[]):
        if line.startswith('#') or not line.strip():
            continue
        cols = line.split('\t')
        if len(cols)>1:
            try:
                val=float(cols[1])
            except ValueError:
                continue
            if max_n is None or val>max_n:
                max_n=val
                max_pos=cols[0]
    row = {
        'filename': basic.get('Filename',''),
        'total_sequences': basic.get('Total Sequences',''),
        'sequence_length': basic.get('Sequence length',''),
        'gc_percent': basic.get('%GC',''),
        'deduplicated_percent': dedup,
        'max_per_base_n_percent': '' if max_n is None else str(max_n),
        'max_n_position': max_pos,
    }
    for key,modname in CHECKS.items():
        row[key]=modules.get(modname,{}).get('status','missing')
    return row


def extract_one(path):
    with zipfile.ZipFile(path) as z:
        names=[name for name in z.namelist() if name.endswith('/fastqc_data.txt') or name=='fastqc_data.txt']
        if len(names)!=1:
            raise ValueError(f'Expected exactly one fastqc_data.txt in {path}; got {len(names)}')
        text=z.read(names[0]).decode('utf-8-sig')
    row=parse_fastqc_text(text)
    row['fastqc_zip']=path.name
    if not row['total_sequences']:
        raise ValueError(f'Missing Total Sequences: {path}')
    return row


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fastqc-dir',type=Path,default=Path('results/qc/raw'))
    parser.add_argument('--out',type=Path,default=Path('docs/qc_metrics.tsv'))
    args=parser.parse_args()
    files=sorted(args.fastqc_dir.glob('*_fastqc.zip'))
    if len(files)!=6:
        raise SystemExit(f'Expected six raw FastQC ZIP files; found {len(files)} in {args.fastqc_dir}.')
    rows=[extract_one(p) for p in files]
    args.out.parent.mkdir(parents=True,exist_ok=True)
    with args.out.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=FIELDS,delimiter='\t')
        w.writeheader()
        w.writerows(rows)
    print(f'Wrote {len(rows)} per-file QC rows to {args.out}')
    print('Review with: column -s $\'\\t\' -t '+str(args.out)+' | less -S')

if __name__=='__main__':
    main()
