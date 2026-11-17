"""Build comprehensive IS110 JSON with all annotations.

For IS110 with circle evidence: use full records from formatter guide JSONs.
For IS110 without: build records from Sniffles consensus + Prodigal ORFs.
Annotate all with DEDD/Tnp20 domains from HMM results.
"""
import csv
import json
import glob
import logging
import os
import sys
from collections import defaultdict

csv.field_size_limit(10_000_000)
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))  # repository root
from cycle.site_config import data_root  # noqa: E402
from cycle.is110_filter import IS110Filter

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(name)s: %(message)s')
logger = logging.getLogger(__name__)

BASE = str(data_root())
HMM_DIR = f"{BASE}/is110_all_006_052"
OUT = f"{BASE}/is110_all_006_052"

# Load IS110 IDs
with open(f"{HMM_DIR}/is110_ids.json") as f:
    is110_ids = set(json.load(f))
logger.info("IS110 IDs: %d", len(is110_ids))

# Load domain annotations
filt = IS110Filter()
protein_hits = {}
for name in ("DEDD", "Tnp20"):
    protein_hits[name] = filt.parse_tblout(f"{HMM_DIR}/{name}_hits.tbl")
_, trans_domains = filt._identify_is110(protein_hits)

# Step 1: Collect full records from formatter guide JSONs
logger.info("Step 1: Collecting from formatter guide JSONs...")
formatter_records = {}
for b in range(6, 53):
    fmt_dir = f"{BASE}/batch_{b:03d}/is_formatter_output"
    if not os.path.isdir(fmt_dir):
        continue
    for jf in sorted(glob.glob(f"{fmt_dir}/*/*_is_records_guide.json")):
        with open(jf) as f:
            records = json.load(f)
        for rec in records:
            is_id = rec.get("is_id", "")
            if is_id in is110_ids and is_id not in formatter_records:
                rec["_batch"] = f"batch_{b:03d}"
                rec["_source"] = "formatter"
                filt._annotate_orf_domains(rec, trans_domains)
                formatter_records[is_id] = rec
logger.info("From formatter: %d", len(formatter_records))

# Step 2: Load Prodigal ORFs for remaining IS110s
missing_ids = is110_ids - set(formatter_records.keys())
logger.info("Step 2: Need Sniffles data for %d IS110s", len(missing_ids))

uuid_orfs = defaultdict(list)
with open(f"{HMM_DIR}/all_proteins.faa") as f:
    header = None
    seq_lines = []
    for line in f:
        if line.startswith(">"):
            if header and header["uuid"] in missing_ids and seq_lines:
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
            if uuid not in missing_ids:
                header = None; seq_lines = []; continue
            strand = "+" if parts[3] == "1" else "-"
            header = {"uuid": uuid, "start": int(parts[1]), "end": int(parts[2]), "strand": strand}
            seq_lines = []
        else:
            if header:
                seq_lines.append(line.strip())
    if header and header["uuid"] in missing_ids and seq_lines:
        seq = "".join(seq_lines).rstrip("*")
        uuid_orfs[header["uuid"]].append({
            "start": header["start"], "end": header["end"],
            "strand": header["strand"], "protein_sequence": seq,
            "length_nt": header["end"] - header["start"] + 1,
        })
logger.info("Loaded ORFs for %d IS110s", len(uuid_orfs))

# Step 3: Load Sniffles consensus + sample info for missing IS110s
logger.info("Step 3: Loading Sniffles data for missing IS110s...")
# Load organism metadata
sample_to_org = {}
for b in range(6, 53):
    for cf in glob.glob(f"{BASE}/batch_{b:03d}/chunks/chunk_*_all.tsv"):
        with open(cf) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                sample_to_org[row["srr_accession"]] = row.get("organism", "")

sniffles_records = {}
for b in range(6, 53):
    snif_dir = f"{BASE}/batch_{b:03d}/sniffles_output"
    dirs_to_scan = []
    if os.path.isdir(snif_dir):
        dirs_to_scan = [(snif_dir, sd) for sd in os.listdir(snif_dir)]
    else:
        for cdir in glob.glob(f"{BASE}/batch_{b:03d}/chunks/chunk_*/sniffles_output"):
            dirs_to_scan.extend([(cdir, sd) for sd in os.listdir(cdir)])

    for parent, sample_dir in dirs_to_scan:
        table_file = os.path.join(parent, sample_dir, f"{sample_dir}.table.txt")
        if not os.path.exists(table_file):
            continue
        with open(table_file) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                uuid = row["UUID"]
                if uuid not in missing_ids or uuid in sniffles_records:
                    continue
                consensus = row.get("Consensus", "")
                orfs = uuid_orfs.get(uuid, [])
                organism = sample_to_org.get(sample_dir, "")

                # Annotate domains on ORFs
                domains_info = trans_domains.get(uuid, {})
                for orf in orfs:
                    orf_domains = []
                    for domain_name in ("DEDD", "Tnp20"):
                        for _, s, e, _ in domains_info.get(domain_name, []):
                            if s == orf["start"] and e == orf["end"]:
                                orf_domains.append(domain_name)
                    if orf_domains:
                        orf["domains"] = sorted(set(orf_domains))

                sniffles_records[uuid] = {
                    "is_id": uuid,
                    "sample_id": sample_dir,
                    "organism": organism,
                    "is_element": {
                        "sequence": consensus,
                        "length": len(consensus),
                        "chrom": row["Chrom"],
                        "start": int(row["Start"]),
                        "end": int(row["End"]),
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
                    "orf_annotation": {
                        "num_orfs": len(orfs),
                        "orfs": orfs,
                    },
                    "guide_hits": [],
                    "guide_summary": {"n_hits": 0, "best_length": 0, "has_revcomp_hit": False},
                    "_batch": f"batch_{b:03d}",
                    "_source": "sniffles",
                }

logger.info("From Sniffles: %d", len(sniffles_records))

# Step 4: Load circle evidence from circle_summary TSVs for sniffles records
logger.info("Step 4: Loading circle evidence...")
for b in range(6, 53):
    for cdir in glob.glob(f"{BASE}/batch_{b:03d}/chunks/chunk_*/circle_output"):
        if not os.path.isdir(cdir):
            continue
        for sample_dir in os.listdir(cdir):
            summary = os.path.join(cdir, sample_dir, f"{sample_dir}_circle_summary.tsv")
            if not os.path.exists(summary):
                continue
            with open(summary) as f:
                reader = csv.DictReader(f, delimiter="\t")
                for row in reader:
                    uuid = row.get("is_uuid", "")
                    if uuid in sniffles_records:
                        sniffles_records[uuid]["circle_evidence"] = {
                            "n_tail_head_reads": int(row.get("n_tail_head_reads", 0)),
                            "n_genome_head_reads": int(row.get("n_genome_head_reads", 0)),
                            "n_tail_genome_reads": int(row.get("n_tail_genome_reads", 0)),
                            "n_total_mapped": int(row.get("n_total_mapped", 0)),
                        }

# Combine and annotate circle_type
all_records = list(formatter_records.values()) + list(sniffles_records.values())
n_with = 0
n_without = 0
for rec in all_records:
    n_th = rec.get("circle_evidence", {}).get("n_tail_head_reads", 0)
    if n_th > 0:
        rec["circle_type"] = "with_circle_evidence"
        n_with += 1
    else:
        rec["circle_type"] = "without_circle_evidence"
        n_without += 1

logger.info("Total: %d records (%d with circle, %d without)",
            len(all_records), n_with, n_without)

# Export single JSON + TSV
with open(f"{OUT}/is110_all_records.json", "w") as f:
    json.dump(all_records, f, indent=2)
logger.info("Wrote %s/is110_all_records.json", OUT)

fields = ["is_id", "sample_id", "organism", "batch", "source", "is_length", "n_orfs",
          "circle_type", "n_tail_head_reads", "n_genome_head_reads", "n_tail_genome_reads",
          "n_guide_hits", "best_guide_length"]
with open(f"{OUT}/is110_all_summary.tsv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    for rec in all_records:
        ce = rec.get("circle_evidence", {})
        gs = rec.get("guide_summary", {})
        writer.writerow({
            "is_id": rec.get("is_id", ""),
            "sample_id": rec.get("sample_id", ""),
            "organism": rec.get("organism", ""),
            "batch": rec.get("_batch", ""),
            "source": rec.get("_source", ""),
            "is_length": rec.get("is_element", {}).get("length", 0),
            "n_orfs": (rec.get("orf_annotation") or {}).get("num_orfs", 0),
            "circle_type": rec.get("circle_type", ""),
            "n_tail_head_reads": ce.get("n_tail_head_reads", 0),
            "n_genome_head_reads": ce.get("n_genome_head_reads", 0),
            "n_tail_genome_reads": ce.get("n_tail_genome_reads", 0),
            "n_guide_hits": gs.get("n_hits", 0),
            "best_guide_length": gs.get("best_length", 0),
        })
logger.info("Wrote %s/is110_all_summary.tsv", OUT)

logger.info("Done.")
