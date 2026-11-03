from pathlib import Path
import csv

root = Path.cwd()
out = root / "docs/final"
fig = root / "docs/figures/final"

files = {
    "30 kb": root / "results/benchmark/NA12878_chr20_happy.summary.csv",
    "200 kb": root / "results/benchmark/NA12878_chr20_validation_happy.summary.csv"
}

beds = {
    "30 kb": root / "results/benchmark/chr20_giab_confident.bed",
    "200 kb": root / "results/benchmark/chr20_independent_confident.bed"
}

for path in list(files.values()) + list(beds.values()):
    if not path.is_file():
        raise SystemExit(f"Missing required file: {path}")

rows = []

for region, path in files.items():
    records = list(csv.DictReader(path.open()))

    for variant_type in ("SNP", "INDEL"):
        matches = [
            r for r in records
            if r["Type"] == variant_type and r["Filter"] == "ALL"
        ]

        if len(matches) != 1:
            raise SystemExit(f"Unexpected benchmark summary: {path}")

        r = matches[0]

        values = [
            int(r[x]) for x in (
                "TRUTH.TOTAL", "TRUTH.TP", "TRUTH.FN",
                "QUERY.TOTAL", "QUERY.FP", "QUERY.UNK"
            )
        ]

        scores = [
            float(r[x]) for x in (
                "METRIC.Precision",
                "METRIC.Recall",
                "METRIC.F1_Score"
            )
        ]

        if values[0] != values[1] + values[2]:
            raise SystemExit("Inconsistent truth variant counts")

        rows.append([region, variant_type, *values, *scores])

confidence = {}

for region, path in beds.items():
    confidence[region] = sum(
        int(fields[2]) - int(fields[1])
        for line in path.read_text().splitlines()
        if line.strip() and not line.startswith("#")
        for fields in [line.split()]
    )

out.mkdir(parents=True, exist_ok=True)
fig.mkdir(parents=True, exist_ok=True)

header = [
    "region", "type", "truth_total", "truth_tp",
    "truth_fn", "query_total", "query_fp",
    "query_unknown", "precision", "recall", "f1"
]

with (out / "BENCHMARK_METRICS.tsv").open("w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(header)
    writer.writerows(rows)

def totals(region):
    selected = [r for r in rows if r[0] == region]
    return {
        "tp": sum(r[3] for r in selected),
        "fn": sum(r[4] for r in selected),
        "fp": sum(r[6] for r in selected),
        "unk": sum(r[7] for r in selected)
    }

small = totals("30 kb")
large = totals("200 kb")

table = "\n".join(
    f"| {r[0]} | {r[1]} | {r[2]} | {r[3]} | "
    f"{r[4]} | {r[6]} | {r[7]} | "
    f"{r[8]:.2f} | {r[9]:.2f} |"
    for r in rows
)

report = f"""# Human NGS Variant Calling — Scientific Report

## Project Overview

This project analyzed publicly available NA12878/HG001
sequencing data using two complementary workflows:

1. Raw FASTQ-to-VCF analysis using a small test dataset.
2. High-coverage regional variant calling and GIAB benchmarking.

## Methods

- Raw sequencing QC: FastQC and MultiQC
- Read preprocessing: Trimmomatic
- Reference alignment: BWA-MEM
- BAM processing: SAMtools
- Germline variant discovery: GATK HaplotypeCaller
- Variant analysis: BCFtools
- Benchmarking: hap.py
- Reference: GRCh38/hg38
- Truth dataset: GIAB HG001 v4.2.1

## Preprocessing

Original read pairs: 30,000
Retained paired reads: 23,848
Paired-read retention: 79.49%

## Low-Coverage Experiment

Chromosome 19:
- Raw candidates: 25
- SNPs: 22
- Indels: 3
- Read depths: 1–2

These calls were not considered reliable genotypes.

## High-Coverage Regional Benchmarking

| Region | Type | Truth total | TP | FN | FP | Unassessed | Precision | Recall |
|---|---|---:|---:|---:|---:|---:|---:|---:|
{table}

### GIAB High-Confidence Regions

30 kb interval: {confidence["30 kb"]:,} confident bases.

200 kb interval: {confidence["200 kb"]:,} confident bases.

### Interpretation

Both regional benchmarks showed high agreement within
the GIAB-evaluable regions.

The 30 kb analysis contained {small["unk"]} unassessed
query units, while the 200 kb analysis contained
{large["unk"]} unassessed query units.

These are not confirmed true or false positives.

## Limitations

- Both regional analyses used the same biological sample.
- Results are not representative of the entire genome.
- Original variant FILTER values were unset.
- The 30 kb region overlaps structural complexity.
- Reference provenance requires final verification.
- No clinical conclusions can be drawn.

## Supporting Documentation

See the existing preprocessing, alignment, variant
quality, and GIAB benchmarking reports in docs/.

## Status

Scientific analyses complete for the defined project.
Reproducibility and publication audit pending.
"""

(out / "FINAL_SCIENTIFIC_REPORT.md").write_text(report)

readme = f"""# Human NGS Variant Calling and GIAB Benchmarking

A bioinformatics portfolio project analyzing real human
NA12878/HG001 next-generation sequencing data.

## Tools

FastQC, MultiQC, Trimmomatic, BWA-MEM, SAMtools,
GATK HaplotypeCaller, BCFtools, Python, and hap.py.

## Key Results

- 23,848 of 30,000 read pairs retained after trimming.
- 30 kb GIAB benchmark: {small["tp"]} truth-side TP.
- 200 kb GIAB benchmark: {large["tp"]} truth-side TP.
- Both regional evaluations reported zero FP and FN.
- Unassessed query units: {small["unk"]} and {large["unk"]}.

Accuracy measurements apply only to evaluable
GIAB high-confidence regions.

## Documentation

[Scientific Report](docs/final/FINAL_SCIENTIFIC_REPORT.md)

[Benchmark Metrics](docs/final/BENCHMARK_METRICS.tsv)

[Data Provenance](docs/final/REPRODUCIBILITY_AND_DATA.md)

## Status

Final publication review pending.
"""

(out / "GITHUB_README_CANDIDATE.md").write_text(readme)

provenance = """# Reproducibility and Data Provenance

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
"""

(out / "REPRODUCIBILITY_AND_DATA.md").write_text(provenance)

audit = """# GitHub Publication Checklist

- [ ] Verify all scientific numbers and interpretation.
- [ ] Review original reference-genome provenance.
- [ ] Check outdated scripts and duplicate reports.
- [ ] Verify script dependencies and relative paths.
- [ ] Inspect all files before Git staging.
- [ ] Exclude large genomic datasets and logs.
- [ ] Preserve benchmark limitations and warnings.
- [ ] Finalize the GitHub README.
"""

(out / "PUBLICATION_AUDIT.md").write_text(audit)

max_tp = max(small["tp"], large["tp"])

svg = f"""<svg xmlns="http://www.w3.org/2000/svg"
width="850" height="260" viewBox="0 0 850 260">
<rect width="850" height="260" fill="white"/>
<g font-family="Arial" fill="#183047">
<text x="25" y="40" font-size="23">
GIAB Truth Variants Recovered
</text>
<text x="25" y="100" font-size="18">30 kb</text>
<text x="25" y="175" font-size="18">200 kb</text>
<rect x="145" y="75" width="{560*small['tp']/max_tp:.1f}"
height="30" fill="#287f87"/>
<rect x="145" y="150" width="{560*large['tp']/max_tp:.1f}"
height="30" fill="#287f87"/>
<text x="730" y="100" font-size="18">{small['tp']}</text>
<text x="730" y="175" font-size="18">{large['tp']}</text>
<text x="25" y="225" font-size="13">
GIAB high-confidence regions only; not genome-wide accuracy.
</text>
</g></svg>"""

(fig / "giab_truth_tp_by_region.svg").write_text(svg)

print("SUCCESS: Final report drafts and figure generated.")
