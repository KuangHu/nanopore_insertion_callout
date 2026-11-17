"""Visualize long IS621 candidates (>2500 bp) with full annotation.

Steps:
1. Get IS621 IDs with strong hits (≥95% id, ≥80% cov)
2. Filter by IS length > 2500 bp
3. Pull full formatter records (with flanking, ORFs, guide hits, circle evidence)
4. Cross-reference with partial circle data
5. Generate PNG + GBK for each
"""
import csv
import json
import glob
import logging
import os
import sys
from collections import defaultdict

csv.field_size_limit(10_000_000)
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))  # repository root
from cycle.site_config import data_root  # noqa: E402
from cycle.is110_filter import IS110Filter
from cycle.visualizer.visualizer import ISElementVisualizer
from cycle.visualizer.genbank import ISElementGenBank

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(name)s: %(message)s')
logger = logging.getLogger(__name__)

BASE = str(data_root())
OUT = f"{BASE}/is621_search/long_visualizations"
os.makedirs(OUT, exist_ok=True)

# Step 1: Get IS621 IDs with strong hits
hits_by_is = defaultdict(list)
with open(f"{BASE}/is621_search/is621_hits.tsv") as f:
    for line in f:
        cols = line.strip().split("\t")
        pident = float(cols[2])
        length = int(cols[3])
        qlen = int(cols[12])
        is_id = cols[1].split("__")[0]
        hits_by_is[is_id].append({"pident": pident, "coverage": length / qlen})

strong_is621 = {is_id for is_id, hits in hits_by_is.items()
                if any(h["pident"] >= 95 and h["coverage"] >= 0.8 for h in hits)}
logger.info("Strong IS621 candidates: %d", len(strong_is621))

# Step 2: Get IS lengths from Sniffles tables, filter > 2500
long_is621 = set()
for b in range(6, 79):
    snif_dir = f"{BASE}/batch_{b:03d}/sniffles_output"
    scan_dirs = []
    if os.path.isdir(snif_dir):
        scan_dirs = [(snif_dir, sd) for sd in os.listdir(snif_dir)]
    else:
        for cdir in glob.glob(f"{BASE}/batch_{b:03d}/chunks/chunk_*/sniffles_output"):
            scan_dirs.extend([(cdir, sd) for sd in os.listdir(cdir)])
    for parent, sample_dir in scan_dirs:
        tf = os.path.join(parent, sample_dir, f"{sample_dir}.table.txt")
        if not os.path.exists(tf):
            continue
        with open(tf) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                uuid = row["UUID"]
                if uuid in strong_is621 and len(row.get("Consensus", "")) > 2500:
                    long_is621.add(uuid)

logger.info("Long IS621 (>2500 bp): %d", len(long_is621))

# Step 3: Collect formatter records (with flanking, ORFs, circle evidence)
formatter_records = {}
for b in range(6, 79):
    fmt_dir = f"{BASE}/batch_{b:03d}/is_formatter_output"
    if not os.path.isdir(fmt_dir):
        continue
    for jf in sorted(glob.glob(f"{fmt_dir}/*/*_is_records_guide.json")):
        with open(jf) as f:
            records = json.load(f)
        for rec in records:
            is_id = rec.get("is_id", "")
            if is_id in long_is621:
                rec["_batch"] = f"batch_{b:03d}"
                formatter_records[is_id] = rec

logger.info("Found formatter records for: %d/%d", len(formatter_records), len(long_is621))
# (records in formatter mean they pass require_th_reads=True → all have full circle)

# Step 4: Cross-reference partial circle data
pc_by_is = defaultdict(list)
for b in range(6, 79):
    pc_dir = f"{BASE}/batch_{b:03d}/partial_circle_output"
    for pcf in glob.glob(f"{pc_dir}/**/*_partial_circle_summary.json", recursive=True):
        with open(pcf) as f:
            data = json.load(f)
        if isinstance(data, list):
            for call in data:
                is_id = call.get("is_id", "")
                if is_id in long_is621:
                    pc_by_is[is_id].append(call)

# Annotate records with partial circle status + domains
filt = IS110Filter()
protein_hits = {}
for name in ("DEDD", "Tnp20"):
    tbl = f"{BASE}/is110_all_006_070/{name}_hits.tbl"
    if os.path.exists(tbl):
        protein_hits[name] = filt.parse_tblout(tbl)
_, trans_domains = filt._identify_is110(protein_hits) if protein_hits else (set(), {})

with_full = 0
with_partial = 0
with_both = 0

for is_id, rec in formatter_records.items():
    has_full = rec.get("circle_evidence", {}).get("n_tail_head_reads", 0) > 0
    has_partial = is_id in pc_by_is
    rec["has_full_circle"] = has_full
    rec["has_partial_circle"] = has_partial
    if has_partial:
        rec["partial_circles"] = pc_by_is[is_id]
    if has_full:
        with_full += 1
    if has_partial:
        with_partial += 1
    if has_full and has_partial:
        with_both += 1
    if trans_domains:
        filt._annotate_orf_domains(rec, trans_domains)

logger.info("Long IS621 with full cycle: %d", with_full)
logger.info("Long IS621 with partial cycle: %d", with_partial)
logger.info("Long IS621 with BOTH: %d", with_both)

# IS621 in long systems NOT in formatter (no full cycle, may or may not have partial)
in_formatter = set(formatter_records.keys())
not_in_formatter = long_is621 - in_formatter
n_with_partial_only = sum(1 for uid in not_in_formatter if uid in pc_by_is)
logger.info("Long IS621 NOT in formatter (no full circle): %d", len(not_in_formatter))
logger.info("  ... but have partial circle: %d", n_with_partial_only)

# Step 5: Generate PNG + GBK for all in formatter (these have flanking data)
vis = ISElementVisualizer()
gb = ISElementGenBank()
n_done = 0
n_err = 0
for is_id, rec in formatter_records.items():
    sample_id = rec["sample_id"]
    tag = f"{sample_id}__{is_id[:8]}"
    png_path = os.path.join(OUT, f"{tag}.png")
    gbk_path = os.path.join(OUT, f"{tag}.gbk")
    alignments = rec.get("guide_hits", [])
    circle_info = rec.get("circle_evidence")
    partial_circles = pc_by_is.get(is_id)
    try:
        vis.save_element_png(rec, alignments, png_path, dpi=150,
                             circle_info=circle_info, partial_circles=partial_circles)
        gb.save_genbank(rec, alignments, gbk_path, circle_info=circle_info)
        n_done += 1
    except Exception as e:
        logger.error("Failed %s: %s", tag, e)
        n_err += 1
logger.info("Visualized: %d PNG+GBK, %d errors", n_done, n_err)

# Save records as JSON for reference
with open(f"{OUT}/long_is621_records.json", "w") as f:
    json.dump(list(formatter_records.values()), f, indent=2)
logger.info("Done. Output: %s", OUT)
