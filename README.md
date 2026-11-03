# Human NGS Variant Calling and GIAB Benchmarking

A bioinformatics portfolio project analyzing real human
NA12878/HG001 next-generation sequencing data.

## Tools

FastQC, MultiQC, Trimmomatic, BWA-MEM, SAMtools,
GATK HaplotypeCaller, BCFtools, Python, and hap.py.

## Key Results

- 23,848 of 30,000 read pairs retained after trimming.
- 30 kb GIAB benchmark: 21 truth-side TP.
- 200 kb GIAB benchmark: 384 truth-side TP.
- Both regional evaluations reported zero FP and FN.
- Unassessed query units: 7 and 49.

Accuracy measurements apply only to evaluable
GIAB high-confidence regions.

## Documentation

[Scientific Report](docs/final/FINAL_SCIENTIFIC_REPORT.md)

[Benchmark Metrics](docs/final/BENCHMARK_METRICS.tsv)

[Data Provenance](docs/final/REPRODUCIBILITY_AND_DATA.md)

## Status

Final publication review pending.
