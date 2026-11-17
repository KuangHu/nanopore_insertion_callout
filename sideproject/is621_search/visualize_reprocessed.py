"""Find IS621 records in the reprocessed samples by protein match (not UUID).

For each reprocessed sample, scan its guide JSON for any IS record whose ORF
matches the IS621 transposase by BLASTp (≥95% id, ≥80% cov). Visualize those.
"""
import csv
import glob
import json
import logging
import os
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))  # repository root
from cycle.site_config import data_root  # noqa: E402
from cycle.visualizer.visualizer import ISElementVisualizer
from cycle.visualizer.genbank import ISElementGenBank

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(name)s: %(message)s')
logger = logging.getLogger(__name__)

BASE = str(data_root())
WORK = f"{BASE}/is621_search/reprocess_work"
OUT = f"{BASE}/is621_search/long_visualizations_flanked"
QUERY = os.path.join(os.path.dirname(os.path.abspath(__file__)), "is621_query.faa")
os.makedirs(OUT, exist_ok=True)

# Step 1: Collect all proteins from reprocessed samples
proteins_fasta = f"{WORK}/all_proteins.faa"
record_by_protein = {}  # protein_id -> (sample_id, is_id, rec, orf_idx)

with open(proteins_fasta, "w") as fout:
    for guide_json in glob.glob(f"{WORK}/*/is_formatter_output/*/*_is_records_guide.json"):
        sample_id = os.path.basename(guide_json).replace("_is_records_guide.json", "")
        with open(guide_json) as f:
            records = json.load(f)
        for rec in records:
            is_id = rec.get("is_id", "")
            orf_ann = rec.get("orf_annotation") or {}
            for i, orf in enumerate(orf_ann.get("orfs", [])):
                seq = orf.get("protein_sequence", "")
                if not seq or len(seq) < 100:
                    continue
                prot_id = f"{sample_id}__{is_id[:8]}__{i}"
                fout.write(f">{prot_id}\n{seq}\n")
                record_by_protein[prot_id] = (sample_id, is_id, rec, i)

logger.info("Total proteins from reprocessed: %d", len(record_by_protein))
if not record_by_protein:
    logger.error("No proteins found, exiting")
    sys.exit(1)

# Step 2: BLASTp IS621 against these
db = f"{WORK}/blast_db"
subprocess.run(["makeblastdb", "-in", proteins_fasta, "-dbtype", "prot", "-out", db],
               check=True, capture_output=True)
blast_out = f"{WORK}/is621_hits.tsv"
subprocess.run([
    "blastp", "-query", QUERY, "-db", db, "-out", blast_out,
    "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen",
    "-evalue", "1e-5", "-num_threads", "8", "-max_target_seqs", "10000",
], check=True, capture_output=True)

# Filter: ≥95% id, ≥80% cov
matching = set()
with open(blast_out) as f:
    for line in f:
        cols = line.strip().split("\t")
        pident = float(cols[2])
        length = int(cols[3])
        qlen = int(cols[12])
        if pident >= 95 and length / qlen >= 0.8:
            matching.add(cols[1])

logger.info("Proteins with strong IS621 match: %d", len(matching))

# Find unique records
unique_records = {}
for prot_id in matching:
    sample_id, is_id, rec, orf_idx = record_by_protein[prot_id]
    # Only keep IS elements > 2500 bp
    if rec.get("is_element", {}).get("length", 0) < 2500:
        continue
    key = f"{sample_id}__{is_id[:8]}"
    if key not in unique_records:
        unique_records[key] = rec

logger.info("Unique long IS621 records: %d", len(unique_records))

# Step 3: Load partial circle data
pc_by_is = {}
for b in range(6, 79):
    pc_dir = f"{BASE}/batch_{b:03d}/partial_circle_output"
    for pcf in glob.glob(f"{pc_dir}/**/*_partial_circle_summary.json", recursive=True):
        with open(pcf) as f:
            data = json.load(f)
        if isinstance(data, list):
            for call in data:
                is_id = call.get("is_id", "")
                if is_id:
                    pc_by_is.setdefault(is_id, []).append(call)

# Step 4: Visualize
vis = ISElementVisualizer()
gb = ISElementGenBank()
n_done = 0
n_err = 0
final_records = []
for tag, rec in unique_records.items():
    png_path = os.path.join(OUT, f"{tag}.png")
    gbk_path = os.path.join(OUT, f"{tag}.gbk")
    alignments = rec.get("guide_hits", [])
    circle_info = rec.get("circle_evidence")
    partial_circles = pc_by_is.get(rec["is_id"])
    rec["_source"] = "reprocessed"
    final_records.append(rec)
    try:
        vis.save_element_png(rec, alignments, png_path, dpi=150,
                             circle_info=circle_info, partial_circles=partial_circles)
        gb.save_genbank(rec, alignments, gbk_path, circle_info=circle_info)
        n_done += 1
        logger.info("Done %s: flanking=%dx%d, guide_hits=%d, n_th=%d",
                    tag,
                    len(rec.get("flanking_upstream", {}).get("sequence", "")),
                    len(rec.get("flanking_downstream", {}).get("sequence", "")),
                    len(alignments),
                    circle_info.get("n_tail_head_reads", 0) if circle_info else 0)
    except Exception as e:
        logger.error("Failed %s: %s", tag, e)
        n_err += 1

logger.info("Visualized: %d PNG+GBK, %d errors", n_done, n_err)
with open(f"{OUT}/long_is621_records_flanked.json", "w") as f:
    json.dump(final_records, f, indent=2)
