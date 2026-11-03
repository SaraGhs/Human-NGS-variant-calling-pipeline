#!/usr/bin/env python3
"""Summarize two hap.py regional benchmarks without changing input VCFs.

Run from the root of the NA12878 portfolio project:
  python scripts/18_finalize_benchmarks.py
Python standard library only.
"""
import csv
from pathlib import Path
from datetime import date
from xml.sax.saxutils import escape

ROOT = Path.cwd()
OUT = ROOT / "docs"
RES = ROOT / "results" / "benchmark"
REGIONS = [
    ("Initial 30 kb", 30001, RES / "NA12878_chr20_happy.summary.csv", RES / "chr20_giab_confident.bed"),
    ("Second 200 kb", 200000, RES / "NA12878_chr20_validation_happy.summary.csv", RES / "chr20_independent_confident.bed"),
]

def numeric(row, key, integer=True):
    value = row.get(key, "")
    if value in (None, "", "NaN"):
        raise ValueError(f"Missing metric {key} in {row}")
    return int(value) if integer else float(value)

def load_rows(path):
    if not path.is_file():
        raise FileNotFoundError(f"Missing hap.py results: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    out = {}
    for kind in ("SNP", "INDEL"):
        match = [r for r in rows if r.get("Type") == kind and r.get("Filter") == "ALL"]
        if len(match) != 1:
            raise ValueError(f"Expected one {kind}/ALL row in {path}, found {len(match)}")
        row = match[0]
        out[kind] = {key: numeric(row, key) for key in (
            "TRUTH.TOTAL", "TRUTH.TP", "TRUTH.FN", "QUERY.TOTAL", "QUERY.FP", "QUERY.UNK")}
        out[kind].update({key: numeric(row, key, integer=False) for key in (
            "METRIC.Recall", "METRIC.Precision", "METRIC.F1_Score")})
        assert out[kind]["TRUTH.TOTAL"] == out[kind]["TRUTH.TP"] + out[kind]["TRUTH.FN"]
    return out

def confident_bases(path):
    if not path.is_file():
        raise FileNotFoundError(f"Missing confident-region BED: {path}")
    intervals = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#"):
            continue
        chrom, start, end, *_ = line.split("\t")
        intervals.setdefault(chrom, []).append((int(start), int(end)))
    total = 0
    for blocks in intervals.values():
        current_end = -1
        for start, end in sorted(blocks):
            if start >= end:
                raise ValueError(f"Invalid interval: {start}-{end}")
            total += max(0, end - max(start, current_end))
            current_end = max(current_end, end)
    return total

report_data = []
for title, size, csv_file, bed_file in REGIONS:
    rows = load_rows(csv_file)
    confident = confident_bases(bed_file)
    if not 0 <= confident <= size:
        raise ValueError(f"Invalid confident-region size for {title}: {confident}/{size}")
    totals = {field: sum(rows[k][field] for k in ("SNP", "INDEL")) for field in (
        "TRUTH.TOTAL", "TRUTH.TP", "TRUTH.FN", "QUERY.TOTAL", "QUERY.FP", "QUERY.UNK")}
    report_data.append(dict(title=title, size=size, confident=confident,
                            percentage=100*confident/size, rows=rows, totals=totals,
                            csv_path=csv_file.relative_to(ROOT)))

OUT.mkdir(exist_ok=True)
(OUT / "figures").mkdir(exist_ok=True)
RES.mkdir(parents=True, exist_ok=True)

# Audit-friendly long-format table.
tsv = RES / "benchmark_metrics_comparison.tsv"
with tsv.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["region", "type", "truth_total", "truth_tp", "truth_fn",
                     "query_total", "query_fp", "query_unk", "recall", "precision", "f1", "giab_confident_pct"])
    for entry in report_data:
        for kind, row in entry["rows"].items():
            writer.writerow([entry["title"], kind, row["TRUTH.TOTAL"], row["TRUTH.TP"],
                             row["TRUTH.FN"], row["QUERY.TOTAL"], row["QUERY.FP"],
                             row["QUERY.UNK"], row["METRIC.Recall"], row["METRIC.Precision"],
                             row["METRIC.F1_Score"], f"{entry['percentage']:.2f}"])

first, second = report_data
lines = [
    "# GIAB HG001 — Two-region small-variant benchmarking results",
    "",
    "## Objective and data sources",
    "",
    "Assess exploratory GATK HaplotypeCaller small-variant results for NA12878/HG001 against "
    "NIST Genome in a Bottle HG001 v4.2.1 GRCh38 truth variants within its high-confidence BED, "
    "using Illumina hap.py. These are two regions from the same biological sample, **not** independent samples.",
    "",
    "## Benchmark scope",
    "",
    "| Region | Interval length | GIAB high-confidence bases | Confident fraction |",
    "|---|---:|---:|---:|",
]
for r in report_data:
    lines.append(f"| {r['title']} | {r['size']:,} | {r['confident']:,} | {r['percentage']:.2f}% |")
lines += ["", "## hap.py ALL-row results", "",
          "| Region | Type | Truth TP / total | FN | Query FP | Query UNK / total | Precision | Recall | F1 |",
          "|---|---|---:|---:|---:|---:|---:|---:|---:|"]
for r in report_data:
    for kind in ("SNP", "INDEL"):
        x = r["rows"][kind]
        lines.append(f"| {r['title']} | {kind} | {x['TRUTH.TP']} / {x['TRUTH.TOTAL']} | "
                     f"{x['TRUTH.FN']} | {x['QUERY.FP']} | {x['QUERY.UNK']} / {x['QUERY.TOTAL']} | "
                     f"{x['METRIC.Precision']:.3f} | {x['METRIC.Recall']:.3f} | {x['METRIC.F1_Score']:.3f} |")
lines += ["", "## Combined descriptive totals", "",
          "| Region | Truth units | Truth TP | Truth FN | Query FP | Query UNK | Query units |",
          "|---|---:|---:|---:|---:|---:|---:|"]
for r in report_data:
    t = r["totals"]
    lines.append(f"| {r['title']} | {t['TRUTH.TOTAL']} | {t['TRUTH.TP']} | "
                 f"{t['TRUTH.FN']} | {t['QUERY.FP']} | {t['QUERY.UNK']} | {t['QUERY.TOTAL']} |")
lines += [
    "", "## Interpretation and limitations", "",
    "- Precision, recall and F1 are **only** for the GIAB-evaluable portions of the specified intervals. They do not estimate genome-wide accuracy.",
    "- `QUERY.UNK` denotes unassessed benchmark units, not automatically true or false calls. Report their number explicitly.",
    "- hap.py uses variant-representation-aware comparison. Its SNP/indel query units need not equal the raw GATK VCF record counts from BCFtools.",
    "- Candidate GATK VCF calls were assessed in the `ALL` category. The raw caller's `FILTER=.` is not an independently validated quality filter.",
    "- Two regions from one biological sample are not biological replication. Region selection and uneven benchmark difficulty constrain generalizability.",
    "- Regional CRAM extraction and reads crossing interval boundaries can affect edge-site evaluation. Avoid overgeneralizing a perfect regional score.",
    "- Source reference FASTA provenance still requires a full provenance note even though chromosome 20 CRAM MD5 compatibility was verified.",
    "", "## Evidence and methods", "",
    "- Truth: GIAB HG001 v4.2.1, GRCh38. Query: GATK HaplotypeCaller on NA12878 30× CRAM regional reads.",
    "- High-confidence BED sets assessable regions; hap.py was run with `-f` for confident regions and `-R` for the target interval.",
    f"- Initial summary: `{first['csv_path']}`.",
    f"- Second summary: `{second['csv_path']}`.",
    "- Machine-readable table: `results/benchmark/benchmark_metrics_comparison.tsv`.",
    "- Diagram: `docs/figures/giab_truth_tp_by_region.svg`.",
    "", "## Remaining before publication", "",
    "1. Record exact hap.py and Conda environment versions, reference provenance and input checksums.",
    "2. Audit the final scripts, data-download commands, `.gitignore`, documentation and generated figures.",
    "3. Publish the code, lightweight summary results and reproducibility instructions only after review.",
]
(OUT / "FINAL_BENCHMARK_COMPARISON.md").write_text("\n".join(lines)+"\n", encoding="utf-8")

# Dependable, dependency-free SVG; all counts are taken from hap.py files.
bars = [(r["title"], typ, r["rows"][typ]["TRUTH.TP"]) for r in report_data for typ in ("SNP", "INDEL")]
maxval = max(v for _,_,v in bars) or 1
svg = [
    '<svg xmlns="http://www.w3.org/2000/svg" width="880" height="330" viewBox="0 0 880 330">',
    '<rect width="880" height="330" fill="#ffffff"/>',
    '<text x="35" y="35" font-family="Arial,sans-serif" font-size="22" font-weight="bold" fill="#15283d">GIAB truth variants recovered in each test interval</text>',
    '<text x="35" y="58" font-family="Arial,sans-serif" font-size="13" fill="#526577">Haplotype-aware comparison within confident regions · NA12878/HG001</text>',
]
for i, (reg, typ, count) in enumerate(bars):
    y = 86 + i*58
    width = int(560*count/maxval)
    label = escape(f"{reg} · {typ}")
    svg += [
        f'<text x="35" y="{y+22}" font-family="Arial,sans-serif" font-size="14" fill="#15283d">{label}</text>',
        f'<rect x="260" y="{y+5}" width="560" height="25" rx="4" fill="#eef2f6"/>',
        f'<rect x="260" y="{y+5}" width="{width}" height="25" rx="4" fill="{"#227f8e" if typ == "SNP" else "#d2a047"}"/>',
        f'<text x="{270+width}" y="{y+23}" font-family="Arial,sans-serif" font-size="14" fill="#15283d">{count}</text>',
    ]
svg += [
    '<text x="35" y="319" font-family="Arial,sans-serif" font-size="12" fill="#526577">Counts are truth-side TP units, not estimates of genome-wide performance.</text>',
    '</svg>',
]
(OUT / "figures" / "giab_truth_tp_by_region.svg").write_text("\n".join(svg)+"\n", encoding="utf-8")
print("Created docs/FINAL_BENCHMARK_COMPARISON.md")
print("Created results/benchmark/benchmark_metrics_comparison.tsv")
print("Created docs/figures/giab_truth_tp_by_region.svg")
for r in report_data:
    t = r["totals"]
    print(f"{r['title']}: truth TP {t['TRUTH.TP']}/{t['TRUTH.TOTAL']}; FP {t['QUERY.FP']}; FN {t['TRUTH.FN']}; query UNK {t['QUERY.UNK']}; confident bases {r['confident']:,}/{r['size']:,}")
