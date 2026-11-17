"""Visualize ALL long IS621 candidates (>2500 bp), regardless of circle evidence.

For records in formatter: use full records (with flanking)
For records only in Sniffles: build from consensus + Prodigal ORFs (no flanking)
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

# Step 1: Strong IS621 IDs
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

# Step 2: Get long ones (>2500 bp) + their Sniffles consensus
long_is621_data = {}  # uuid -> {sample_id, consensus, chrom, start, end, batch}
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
                if uuid not in strong_is621 or uuid in long_is621_data:
                    continue
                consensus = row.get("Consensus", "")
                if len(consensus) > 2500:
                    long_is621_data[uuid] = {
                        "sample_id": sample_dir,
                        "consensus": consensus,
                        "chrom": row["Chrom"],
                        "start": int(row["Start"]),
                        "end": int(row["End"]),
                        "batch": f"batch_{b:03d}",
                    }
logger.info("Long IS621 (>2500 bp): %d", len(long_is621_data))

# Step 3: Try to get formatter records (have flanking + assembled IS)
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
            if is_id in long_is621_data:
                rec["_batch"] = f"batch_{b:03d}"
                formatter_records[is_id] = rec

# Step 4: Get Prodigal ORFs for those not in formatter
need_orfs = set(long_is621_data.keys()) - set(formatter_records.keys())
uuid_orfs = defaultdict(list)
with open(f"{BASE}/is110_all_006_070/all_proteins.faa") as f:
    header = None
    seq_lines = []
    for line in f:
        if line.startswith(">"):
            if header and header["uuid"] in need_orfs and seq_lines:
                seq = "".join(seq_lines).rstrip("*")
                uuid_orfs[header["uuid"]].append({
                    "start": header["start"], "end": header["end"],
                    "strand": header["strand"], "protein_sequence": seq,
                    "length_nt": header["end"] - header["start"] + 1,
                })
            parts = line[1:].strip().split(" # ")
            if len(parts) < 4:
                header = None; seq_lines = []; continue
            prot_id = parts[0]
            uuid = "_".join(prot_id.rsplit("_", 1)[:-1])
            if uuid not in need_orfs:
                header = None; seq_lines = []; continue
            strand = "+" if parts[3] == "1" else "-"
            header = {"uuid": uuid, "start": int(parts[1]), "end": int(parts[2]), "strand": strand}
            seq_lines = []
        else:
            if header:
                seq_lines.append(line.strip())
    if header and header["uuid"] in need_orfs and seq_lines:
        seq = "".join(seq_lines).rstrip("*")
        uuid_orfs[header["uuid"]].append({
            "start": header["start"], "end": header["end"],
            "strand": header["strand"], "protein_sequence": seq,
            "length_nt": header["end"] - header["start"] + 1,
        })

# Step 5: Build records for those not in formatter
for is_id in need_orfs:
    data = long_is621_data[is_id]
    orfs = uuid_orfs.get(is_id, [])
    formatter_records[is_id] = {
        "is_id": is_id,
        "sample_id": data["sample_id"],
        "organism": "",
        "is_element": {
            "sequence": data["consensus"],
            "length": len(data["consensus"]),
            "chrom": data["chrom"],
            "start": data["start"],
            "end": data["end"],
        },
        "flanking_upstream": {"sequence": ""},
        "flanking_downstream": {"sequence": ""},
        "insertion_site": {},
        "circle_evidence": {
            "n_tail_head_reads": 0,
            "n_genome_head_reads": 0,
            "n_tail_genome_reads": 0,
            "n_total_mapped": 0,
        },
        "assembly_method": "sniffles_consensus",
        "orf_annotation": {"num_orfs": len(orfs), "orfs": orfs},
        "guide_hits": [],
        "guide_summary": {"n_hits": 0, "best_length": 0, "has_revcomp_hit": False},
        "_batch": data["batch"],
        "_source": "sniffles",
    }

# Step 6: Load partial circle data
pc_by_is = defaultdict(list)
for b in range(6, 79):
    pc_dir = f"{BASE}/batch_{b:03d}/partial_circle_output"
    for pcf in glob.glob(f"{pc_dir}/**/*_partial_circle_summary.json", recursive=True):
        with open(pcf) as f:
            data = json.load(f)
        if isinstance(data, list):
            for call in data:
                is_id = call.get("is_id", "")
                if is_id in long_is621_data:
                    pc_by_is[is_id].append(call)

# Step 7: Visualize all
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

with open(f"{OUT}/long_is621_records.json", "w") as f:
    json.dump(list(formatter_records.values()), f, indent=2)
logger.info("Done. Output: %s", OUT)
