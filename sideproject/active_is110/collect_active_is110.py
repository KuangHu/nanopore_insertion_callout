"""Collect IS110 candidates with circle activity (full or partial cycle evidence).

Pipeline:
1. Run hmmsearch on Prodigal proteins from all batches → IS110 IDs (DEDD + Tnp20)
2. Cross-reference with circle summary TSVs → tail-head reads (full cycle)
3. Cross-reference with partial circle output → partial cycle calls
4. Collect IS110 records that have TH>0 OR partial circle
5. Build full standard records (flanking, ORFs, guide hits, etc.) per record
"""
import csv
import json
import glob
import logging
import os
import sys
from collections import defaultdict

csv.field_size_limit(10_000_000)
sys.path.insert(0, '/home/kuangh/tools/Cycle')
from cycle.is110_filter import IS110Filter

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(name)s: %(message)s')
logger = logging.getLogger(__name__)

BASE = "/groups/rubin/projects/kuang/out/IS_cycle"
OUT = f"{BASE}/active_is110"
os.makedirs(OUT, exist_ok=True)

# Range: 006-070
BATCH_RANGE = range(6, 71)


def get_batch_dirs(batch_num):
    """Return list of (sniffles_dir, parent_label) tuples for a batch."""
    dirs = []
    bdir = f"{BASE}/batch_{batch_num:03d}"
    snif = f"{bdir}/sniffles_output"
    if os.path.isdir(snif):
        dirs.append((snif, "main"))
    for cdir in sorted(glob.glob(f"{bdir}/chunks/chunk_*/sniffles_output")):
        dirs.append((cdir, "chunk"))
    return dirs


# Step 1: Build protein FASTA from all annotated JSONs (Prodigal proteins)
logger.info("Step 1: Building protein FASTA from formatter output...")
fasta_path = f"{OUT}/all_proteins.faa"
n_proteins = 0
n_records_seen = 0
seen_ids = set()
with open(fasta_path, "w") as fout:
    for b in BATCH_RANGE:
        fmt_dir = f"{BASE}/batch_{b:03d}/is_formatter_output"
        if not os.path.isdir(fmt_dir):
            continue
        for jf in sorted(glob.glob(f"{fmt_dir}/*/*_is_records_annotated.json")):
            with open(jf) as f:
                records = json.load(f)
            for rec in records:
                n_records_seen += 1
                is_id = rec.get("is_id", "")
                if not is_id or is_id in seen_ids:
                    continue
                seen_ids.add(is_id)
                orf_ann = rec.get("orf_annotation")
                if not orf_ann:
                    continue
                for orf in orf_ann.get("orfs", []):
                    seq = orf.get("protein_sequence", "")
                    if not seq:
                        continue
                    header = f"{is_id}__{orf['start']}_{orf['end']}_{orf['strand']}"
                    fout.write(f">{header}\n{seq}\n")
                    n_proteins += 1

logger.info("Wrote %d proteins from %d unique IS records (%d total scanned)",
            n_proteins, len(seen_ids), n_records_seen)

# Step 2: Run hmmsearch
logger.info("Step 2: Running hmmsearch...")
filt = IS110Filter()
is110_ids, trans_domains = filt.run_hmmsearch(fasta_path, OUT, cpus=8)
with open(f"{OUT}/is110_ids.json", "w") as f:
    json.dump(list(is110_ids), f)
logger.info("IS110 IDs: %d", len(is110_ids))

# Step 3: Load partial circle data, keyed by is_id
logger.info("Step 3: Loading partial circle data...")
pc_by_is = {}
for b in BATCH_RANGE:
    pc_dir = f"{BASE}/batch_{b:03d}/partial_circle_output"
    for pcf in glob.glob(f"{pc_dir}/**/*_partial_circle_summary.json", recursive=True):
        with open(pcf) as f:
            data = json.load(f)
        if isinstance(data, list):
            for call in data:
                is_id = call.get("is_id", "")
                if is_id:
                    pc_by_is.setdefault(is_id, []).append(call)
logger.info("Partial circle data: %d IS elements", len(pc_by_is))

# Step 4: Collect IS110 records with TH>0 or partial circle
logger.info("Step 4: Collecting IS110 records with circle activity...")
active_records = []
seen_active = set()
n_active_full = 0
n_active_partial_only = 0

for b in BATCH_RANGE:
    fmt_dir = f"{BASE}/batch_{b:03d}/is_formatter_output"
    if not os.path.isdir(fmt_dir):
        continue
    for jf in sorted(glob.glob(f"{fmt_dir}/*/*_is_records_guide.json")):
        with open(jf) as f:
            records = json.load(f)
        for rec in records:
            is_id = rec.get("is_id", "")
            if is_id not in is110_ids or is_id in seen_active:
                continue

            ce = rec.get("circle_evidence", {})
            n_th = ce.get("n_tail_head_reads", 0)
            has_full = n_th > 0
            has_partial = is_id in pc_by_is

            if not (has_full or has_partial):
                continue

            seen_active.add(is_id)
            rec["_batch"] = f"batch_{b:03d}"
            rec["_circle_type"] = "full" if has_full else "partial"
            rec["has_full_circle"] = has_full
            rec["has_partial_circle"] = has_partial
            if has_partial:
                rec["partial_circles"] = pc_by_is[is_id]
            filt._annotate_orf_domains(rec, trans_domains)
            active_records.append(rec)

            if has_full:
                n_active_full += 1
            else:
                n_active_partial_only += 1

logger.info("Active IS110: %d total (%d with full, %d partial-only)",
            len(active_records), n_active_full, n_active_partial_only)

# Step 5: Export
out_json = f"{OUT}/active_is110_records.json"
with open(out_json, "w") as f:
    json.dump(active_records, f, indent=2)
logger.info("Wrote %s", out_json)

# Summary TSV
fields = ["is_id", "sample_id", "organism", "batch", "is_length", "n_orfs",
          "circle_type", "has_full_circle", "has_partial_circle",
          "n_tail_head_reads", "n_partial_circles",
          "n_guide_hits", "best_guide_length"]
out_tsv = f"{OUT}/active_is110_summary.tsv"
with open(out_tsv, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    for rec in active_records:
        ce = rec.get("circle_evidence", {})
        gs = rec.get("guide_summary", {})
        writer.writerow({
            "is_id": rec.get("is_id", ""),
            "sample_id": rec.get("sample_id", ""),
            "organism": rec.get("organism", ""),
            "batch": rec.get("_batch", ""),
            "is_length": rec.get("is_element", {}).get("length", 0),
            "n_orfs": (rec.get("orf_annotation") or {}).get("num_orfs", 0),
            "circle_type": rec.get("_circle_type", ""),
            "has_full_circle": rec.get("has_full_circle", False),
            "has_partial_circle": rec.get("has_partial_circle", False),
            "n_tail_head_reads": ce.get("n_tail_head_reads", 0),
            "n_partial_circles": len(rec.get("partial_circles", [])),
            "n_guide_hits": gs.get("n_hits", 0),
            "best_guide_length": gs.get("best_length", 0),
        })
logger.info("Wrote %s", out_tsv)

logger.info("Done.")
