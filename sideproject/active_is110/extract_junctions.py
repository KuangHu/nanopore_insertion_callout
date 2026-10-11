"""Extract junction sequences for full and partial cycle IS110 elements.

For each active IS110:
- Full cycle junction: [last 20bp of IS] + [first 20bp of IS]
- Partial cycle junction(s): [seq[end-20:end]] + [seq[start:start+20]]
  for each partial circle [start, end] breakpoint

Output: TSV with junction sequences and FASTA file.
"""
import csv
import json
import os
import sys

BASE = "/groups/rubin/projects/kuang/out/IS_cycle/active_is110"

with open(f"{BASE}/active_is110_records.json") as f:
    records = json.load(f)

WINDOW = 20

junctions = []

for rec in records:
    is_id = rec["is_id"]
    sample_id = rec["sample_id"]
    organism = rec.get("organism", "")
    seq = rec.get("is_element", {}).get("sequence", "")
    is_len = len(seq)

    if not seq or is_len < 2 * WINDOW:
        continue

    # Full circle junction
    if rec.get("has_full_circle"):
        upstream = seq[-WINDOW:]   # last 20bp of IS (= just before junction)
        downstream = seq[:WINDOW]  # first 20bp of IS (= just after junction)
        junction_seq = upstream + downstream
        junctions.append({
            "is_id": is_id,
            "sample_id": sample_id,
            "organism": organism,
            "junction_type": "full_circle",
            "circle_start": 0,
            "circle_end": is_len,
            "circle_size": is_len,
            "upstream_20bp": upstream,
            "downstream_20bp": downstream,
            "junction_seq": junction_seq,
        })

    # Partial circle junctions
    for pc in rec.get("partial_circles", []):
        cs = pc.get("circle_start", 0)
        ce = pc.get("circle_end", 0)
        if ce - cs < 2 * WINDOW:
            continue
        if ce > is_len:
            ce = is_len
        if cs < 0:
            cs = 0
        # Upstream = seq right before the back-jump (end of partial circle)
        upstream = seq[max(0, ce - WINDOW):ce]
        # Downstream = seq at the start of the partial circle
        downstream = seq[cs:cs + WINDOW]
        if len(upstream) < WINDOW or len(downstream) < WINDOW:
            continue
        junction_seq = upstream + downstream
        junctions.append({
            "is_id": is_id,
            "sample_id": sample_id,
            "organism": organism,
            "junction_type": "partial_circle",
            "circle_start": cs,
            "circle_end": ce,
            "circle_size": ce - cs,
            "n_supporting_reads": pc.get("n_supporting_reads", 0),
            "upstream_20bp": upstream,
            "downstream_20bp": downstream,
            "junction_seq": junction_seq,
        })

print(f"Total junctions: {len(junctions)}")
print(f"  full circle: {sum(1 for j in junctions if j['junction_type'] == 'full_circle')}")
print(f"  partial circle: {sum(1 for j in junctions if j['junction_type'] == 'partial_circle')}")

# Write TSV
tsv_path = f"{BASE}/junction_sites.tsv"
fields = ["is_id", "sample_id", "organism", "junction_type",
          "circle_start", "circle_end", "circle_size",
          "n_supporting_reads", "upstream_20bp", "downstream_20bp", "junction_seq"]
with open(tsv_path, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
    writer.writeheader()
    for j in junctions:
        writer.writerow(j)
print(f"Wrote {tsv_path}")

# Write FASTA
fa_path = f"{BASE}/junction_sites.fna"
with open(fa_path, "w") as f:
    for j in junctions:
        header = f">{j['is_id'][:8]}_{j['sample_id']}_{j['junction_type']}_{j['circle_start']}-{j['circle_end']}"
        f.write(f"{header}\n{j['junction_seq']}\n")
print(f"Wrote {fa_path}")
