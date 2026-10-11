"""Within each IS110 system, compare the 5bp+5bp breakpoint motif across junctions.

For each IS110 with multiple junctions:
- Extract 5bp upstream + 5bp downstream of each junction (10bp total)
- Compute pairwise identity within the system
- Report systems with very similar vs very different breakpoint motifs
"""
import csv
import json
from collections import defaultdict
from itertools import combinations

BASE = "/groups/rubin/projects/kuang/out/IS_cycle/active_is110"

# Load junctions
junctions_by_is = defaultdict(list)
with open(f"{BASE}/junction_sites.tsv") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        is_id = row["is_id"]
        upstream = row["upstream_20bp"]
        downstream = row["downstream_20bp"]
        # 5bp+5bp = last 5bp upstream + first 5bp downstream
        motif_5 = upstream[-5:] + downstream[:5]
        junctions_by_is[is_id].append({
            "type": row["junction_type"],
            "start": int(row["circle_start"]),
            "end": int(row["circle_end"]),
            "size": int(row["circle_size"]),
            "motif_5": motif_5,
            "motif_full": upstream + downstream,
        })


def hamming_identity(s1: str, s2: str) -> float:
    if len(s1) != len(s2):
        return 0.0
    matches = sum(1 for a, b in zip(s1, s2) if a == b)
    return matches / len(s1)


# For each IS110 with 2+ junctions, compute pairwise motif identity
results = []
for is_id, juncs in junctions_by_is.items():
    if len(juncs) < 2:
        continue
    motifs_5 = [j["motif_5"] for j in juncs]
    motifs_full = [j["motif_full"] for j in juncs]
    types = [j["type"] for j in juncs]

    # Pairwise identity (5bp+5bp)
    ids_5 = [hamming_identity(a, b) for a, b in combinations(motifs_5, 2)]
    ids_full = [hamming_identity(a, b) for a, b in combinations(motifs_full, 2)]

    # Number of unique 5bp+5bp motifs
    n_unique_5 = len(set(motifs_5))

    results.append({
        "is_id": is_id,
        "n_junctions": len(juncs),
        "n_full": types.count("full_circle"),
        "n_partial": types.count("partial_circle"),
        "n_unique_motif_5": n_unique_5,
        "min_motif_5_identity": min(ids_5) if ids_5 else 1.0,
        "mean_motif_5_identity": sum(ids_5) / len(ids_5) if ids_5 else 1.0,
        "max_motif_5_identity": max(ids_5) if ids_5 else 1.0,
        "min_motif_full_identity": min(ids_full) if ids_full else 1.0,
        "mean_motif_full_identity": sum(ids_full) / len(ids_full) if ids_full else 1.0,
    })

print(f"IS110 with 2+ junctions: {len(results)}")
print()

# Distribution of motif identity within systems
buckets = {"identical (100%)": 0, "very similar (>=80%)": 0,
           "similar (>=60%)": 0, "different (<60%)": 0}
for r in results:
    mean_id = r["mean_motif_5_identity"]
    if mean_id == 1.0:
        buckets["identical (100%)"] += 1
    elif mean_id >= 0.8:
        buckets["very similar (>=80%)"] += 1
    elif mean_id >= 0.6:
        buckets["similar (>=60%)"] += 1
    else:
        buckets["different (<60%)"] += 1

print("Distribution of mean 5bp+5bp motif identity within IS110 systems:")
for k, v in buckets.items():
    pct = 100 * v / len(results) if results else 0
    print(f"  {k:30s}: {v:5d} ({pct:.1f}%)")

# Save results
out_tsv = f"{BASE}/junction_motif_comparison.tsv"
fields = ["is_id", "n_junctions", "n_full", "n_partial", "n_unique_motif_5",
          "min_motif_5_identity", "mean_motif_5_identity", "max_motif_5_identity",
          "min_motif_full_identity", "mean_motif_full_identity"]
with open(out_tsv, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    for r in sorted(results, key=lambda x: -x["n_junctions"]):
        writer.writerow({k: f"{v:.3f}" if isinstance(v, float) else v for k, v in r.items()})
print(f"\nWrote {out_tsv}")

# Show top 10 systems with most junctions
print("\nTop 10 IS110 systems by number of junctions:")
print(f"{'IS ID':<16} {'#junc':>6} {'#full':>6} {'#part':>6} {'#unique':>8} {'mean ident':>10}")
for r in sorted(results, key=lambda x: -x["n_junctions"])[:10]:
    print(f"{r['is_id'][:16]:<16} {r['n_junctions']:>6} {r['n_full']:>6} "
          f"{r['n_partial']:>6} {r['n_unique_motif_5']:>8} {r['mean_motif_5_identity']:>10.3f}")

# Show 10 examples with all-identical motifs
print("\nExamples of IS110 with identical 5bp+5bp motifs across all junctions:")
identical = [r for r in results if r["mean_motif_5_identity"] == 1.0]
for r in identical[:5]:
    juncs = junctions_by_is[r["is_id"]]
    print(f"  {r['is_id'][:16]} ({r['n_junctions']} junctions, motif: {juncs[0]['motif_5']})")

# Show 10 examples with very different motifs
print("\nExamples of IS110 with very different motifs (mean identity < 0.5):")
diff = sorted([r for r in results if r["mean_motif_5_identity"] < 0.5],
              key=lambda x: x["mean_motif_5_identity"])
for r in diff[:5]:
    juncs = junctions_by_is[r["is_id"]]
    motifs = [j["motif_5"] for j in juncs]
    print(f"  {r['is_id'][:16]} ({r['n_junctions']} junctions, mean ident: {r['mean_motif_5_identity']:.2f})")
    for m, j in zip(motifs[:5], juncs[:5]):
        print(f"    {m}  type={j['type']}  start={j['start']} end={j['end']}")
