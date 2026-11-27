#!/usr/bin/env python3
"""Read-only checks on pipeline outputs.

usage: check_outputs.py ROOT OUTDIR [--tables DIR]

ROOT is searched recursively, so it can be one run directory, one batch or
the whole data root. Writes to OUTDIR:

  hard_clips.tsv      supplementary records with a hard clip, per BAM
  guide_lengths.tsv   guide hit counts by orientation and length
  strands.tsv         formatter records by strand
  samples.tsv         one row per sample: files found, validity, record counts
  references.tsv      reference used for each aligned sample
  batches.tsv         finished samples per batch table (only with --tables)
"""
import argparse
import csv
import json
from collections import Counter
from pathlib import Path

import pysam

parser = argparse.ArgumentParser()
parser.add_argument("root", type=Path)
parser.add_argument("outdir", type=Path)
parser.add_argument("--tables", type=Path, help="directory of batch_*.tsv sample tables")
args = parser.parse_args()
args.outdir.mkdir(parents=True, exist_ok=True)


def write(name, header, rows):
    with open(args.outdir / name, "w", newline="") as handle:
        out = csv.writer(handle, delimiter="\t")
        out.writerow(header)
        out.writerows(rows)
    print(f"{name}: {len(rows)} rows")


def sample_of(path):
    return path.name.split(".")[0].split("_is_records")[0].split("_circle")[0].split("_partial")[0]


samples = {}


def row(sample):
    return samples.setdefault(sample, {})


# 1. Hard clips in supplementary records
clip_rows = []
for kind, pattern in [
    ("alignment", "**/alignments/*.sorted.bam"),
    ("circle", "**/circle_output/*/*.circle.sorted.bam"),
    ("partial_circle", "**/partial_circle_output/*/*.partial.nsorted.bam"),
]:
    for bam_path in sorted(args.root.glob(pattern)):
        n_records = n_supp = n_hard = 0
        with pysam.AlignmentFile(str(bam_path), "rb", check_sq=False) as bam:
            for read in bam.fetch(until_eof=True):
                n_records += 1
                if read.is_supplementary:
                    n_supp += 1
                    if read.cigartuples and (read.cigartuples[0][0] == 5 or read.cigartuples[-1][0] == 5):
                        n_hard += 1
        clip_rows.append([kind, sample_of(bam_path), n_records, n_supp, n_hard, bam_path])
write("hard_clips.tsv", ["bam_kind", "sample", "records", "supplementary", "supplementary_hard_clipped", "path"], clip_rows)

# 2 and 3. Guide lengths by orientation, formatter strands; JSON validity
guide_counts = Counter()
strand_counts = Counter()
for path in sorted(args.root.glob("**/is_formatter_output/*/*_is_records.json")):
    sample = sample_of(path)
    try:
        records = json.load(open(path))
    except ValueError:
        row(sample)["formatter_json"] = "unparseable"
        continue
    row(sample)["formatter_json"] = "ok"
    row(sample)["formatter_records"] = len(records)
    row(sample)["records_with_sequence"] = sum(1 for r in records if r["is_element"]["sequence"])
    for record in records:
        strand_counts[(record.get("assembly_method"), record["is_element"].get("strand"))] += 1
for path in sorted(args.root.glob("**/is_formatter_output/*/*_is_records_annotated.json")):
    try:
        json.load(open(path))
        row(sample_of(path))["annotated_json"] = "ok"
    except ValueError:
        row(sample_of(path))["annotated_json"] = "unparseable"
for path in sorted(args.root.glob("**/is_formatter_output/*/*_is_records_guide.json")):
    sample = sample_of(path)
    try:
        records = json.load(open(path))
    except ValueError:
        row(sample)["guide_json"] = "unparseable"
        continue
    row(sample)["guide_json"] = "ok"
    hits = [hit for record in records for hit in record.get("guide_hits") or []]
    row(sample)["guide_hits"] = len(hits)
    for hit in hits:
        guide_counts[(hit["orientation"], hit["length"], hit["mismatches"])] += 1
write("guide_lengths.tsv", ["orientation", "length", "mismatches", "hits"],
      [[*key, n] for key, n in sorted(guide_counts.items())])
write("strands.tsv", ["assembly_method", "strand", "records"],
      [[*key, n] for key, n in sorted(strand_counts.items(), key=str)])

# 4. Completeness of each sample's files
for path in sorted(args.root.glob("**/sra_downloads/download_status.tsv")):
    for line in csv.DictReader(open(path), delimiter="\t"):
        row(line["srr_accession"])["download"] = line["download_status"]
reference_rows = []
for path in sorted(args.root.glob("**/alignments/alignment_status.tsv")):
    for line in csv.DictReader(open(path), delimiter="\t"):
        row(line["sample_id"])["alignment"] = line["status"]
        reference = Path(line["reference"]) if line["reference"] else None
        reference_rows.append([line["sample_id"], reference.parent.name if reference else "", line["reference"]])
write("references.tsv", ["sample", "accession", "fasta"], reference_rows)

for path in sorted(args.root.glob("**/sniffles_output/*/*.vcf")):
    sample = sample_of(path)
    text = open(path).read()
    lines = text.split("\n")
    body = [line for line in lines if line and not line.startswith("#")]
    header_ok = any(line.startswith("#CHROM") for line in lines)
    last_ok = not body or len(body[-1].split("\t")) >= 8
    complete = header_ok and text.endswith("\n") and last_ok
    row(sample)["vcf"] = "ok" if complete else "incomplete"
    row(sample)["vcf_records"] = len(body)
    row(sample)["vcf_symbolic_ins"] = sum(1 for line in body if line.split("\t")[4:5] == ["<INS>"])
    row(sample)["vcf_bnd"] = sum(1 for line in body if "SVTYPE=BND" in line)
for path in sorted(args.root.glob("**/sniffles_output/*/*.table.txt")):
    row(sample_of(path))["sniffles_table_rows"] = sum(1 for _ in open(path)) - 1
for path in sorted(args.root.glob("**/circle_output/*/*_circle_summary.tsv")):
    row(sample_of(path))["circle_rows"] = sum(1 for _ in open(path)) - 1
for path in sorted(args.root.glob("**/partial_circle_output/*/*_partial_circle_summary.json")):
    try:
        row(sample_of(path))["partial_circle_calls"] = len(json.load(open(path)))
    except ValueError:
        row(sample_of(path))["partial_circle_calls"] = "unparseable"

steps = ["download", "alignment", "vcf", "sniffles_table_rows", "circle_rows", "formatter_json",
         "annotated_json", "guide_json", "partial_circle_calls"]
columns = ["download", "alignment", "vcf", "vcf_records", "vcf_symbolic_ins", "vcf_bnd", "sniffles_table_rows",
           "circle_rows", "formatter_json", "formatter_records", "records_with_sequence", "annotated_json",
           "guide_json", "guide_hits", "partial_circle_calls"]
sample_rows = []
for sample in sorted(samples):
    found = samples[sample]
    reached = "none"
    for step in steps:
        if found.get(step) in (None, "", "unparseable", "incomplete") or (step in ("download", "alignment") and found[step] != "ok"):
            break
        reached = step
    sample_rows.append([sample, reached] + [found.get(column, "") for column in columns])
write("samples.tsv", ["sample", "last_step_with_output"] + columns, sample_rows)

if args.tables:
    batch_rows = []
    for table in sorted(args.tables.glob("batch_*.tsv")):
        wanted = [line["srr_accession"] for line in csv.DictReader(open(table), delimiter="\t")]
        batch_rows.append([
            table.stem, len(wanted),
            sum(1 for s in wanted if samples.get(s, {}).get("vcf") == "ok"),
            sum(1 for s in wanted if samples.get(s, {}).get("vcf") == "incomplete"),
            sum(1 for s in wanted if samples.get(s, {}).get("guide_json") == "ok"),
            sum(1 for s in wanted if samples.get(s, {}).get("guide_json") == "unparseable"),
        ])
    write("batches.tsv", ["batch", "samples", "vcf_ok", "vcf_incomplete", "guide_json_ok", "guide_json_unparseable"], batch_rows)
