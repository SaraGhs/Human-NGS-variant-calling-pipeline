# Reproducibility and Data Provenance

## Public Datasets

Raw FASTQ:
s3://gatk-test-data/wgs_fastq/NA12878_20k/

High-coverage CRAM:
https://s3.amazonaws.com/1000genomes/1000G_2504_high_coverage/data/ERR3239334/NA12878.final.cram

Truth benchmark:
GIAB HG001 v4.2.1, GRCh38.

## Reference
Chromosome 20 MD5 compatibility confirmed.
Original full-reference provenance requires review.

## Environment
Separate Conda environments were used for the
NGS pipeline and hap.py benchmarking.

## Publication
Large raw reads, BAM, CRAM and VCF datasets
should not be committed to GitHub.
