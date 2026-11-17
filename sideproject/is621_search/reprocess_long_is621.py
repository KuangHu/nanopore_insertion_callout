"""Re-process 43 long IS621 samples to get flanking sequences.

For each sample:
1. Find original batch and reference
2. Download FASTQ, align, run Sniffles (if not cached), format with require_th_reads=False
3. Extract record for the IS621 UUID
4. Clean up FASTQ/BAM
5. Visualize
"""
import csv
import json
import glob
import logging
import os
import sys
import shutil
import subprocess
from pathlib import Path
from collections import defaultdict

csv.field_size_limit(10_000_000)
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))  # repository root
from cycle.site_config import data_root  # noqa: E402
from cycle.download_manager.sra_download import SRADownloader
from cycle.preprocess.aligner import Aligner
from cycle.preprocess.sniffles_runner import SnifflesRunner
from cycle.circle_detect.circle_finder import CircleFinder
from cycle.is_formatter.formatter import ISFormatter
from cycle.orf_annotator.annotator import ORFAnnotator
from cycle.guide_finder.finder import GuideFinder
from cycle.visualizer.visualizer import ISElementVisualizer
from cycle.visualizer.genbank import ISElementGenBank
from cycle.utils import find_fastq

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(name)s: %(message)s')
logger = logging.getLogger(__name__)

BASE = str(data_root())
OUT = f"{BASE}/is621_search/long_visualizations_flanked"
WORK = f"{BASE}/is621_search/reprocess_work"
os.makedirs(OUT, exist_ok=True)
os.makedirs(WORK, exist_ok=True)

# Step 1: Load the 43 long IS621 candidates (those without formatter records)
with open(f"{BASE}/is621_search/long_visualizations/long_is621_records.json") as f:
    existing_records = json.load(f)
need_reprocess = [r for r in existing_records if r.get("_source") == "sniffles"]
logger.info("Need to reprocess: %d samples", len(need_reprocess))

# Step 2: Build sample → batch + reference map from chunk TSVs
sample_meta = {}
for b in range(6, 79):
    for cf in glob.glob(f"{BASE}/batch_{b:03d}/chunks/chunk_*_all.tsv"):
        with open(cf) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                sample_meta[row["srr_accession"]] = {
                    "batch": f"batch_{b:03d}",
                    "organism": row.get("organism", ""),
                }

# Also try metadata_for_sniffles.tsv for reference paths
sample_ref = {}
for b in range(6, 79):
    meta_file = f"{BASE}/batch_{b:03d}/metadata_for_sniffles.tsv"
    if os.path.exists(meta_file):
        with open(meta_file) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                sample_ref[row["srr_accession"]] = row["reference"]

# For samples without metadata_for_sniffles, find reference from existing reference_genomes
def find_reference(sample_id, batch):
    """Try to find the reference for this sample."""
    if sample_id in sample_ref:
        ref = sample_ref[sample_id]
        if os.path.exists(ref):
            return ref
    # Try to use the existing reference_genomes folder
    ref_base = f"{BASE}/{batch}/reference_genomes"
    if not os.path.isdir(ref_base):
        return None
    # Get all references and use the first one (might be wrong - need better mapping)
    refs = sorted(glob.glob(f"{ref_base}/*/*.fna"))
    return refs[0] if refs else None

# Step 3: Process each sample
vis = ISElementVisualizer()
gb = ISElementGenBank()
filt_records = []
n_done = 0
n_err = 0

for rec in need_reprocess:
    is_id = rec["is_id"]
    sample_id = rec["sample_id"]
    meta = sample_meta.get(sample_id, {})
    batch = meta.get("batch", "")
    organism = meta.get("organism", "")

    if not batch:
        logger.warning("Cannot find batch for %s, skipping", sample_id)
        continue

    # Find reference for this sample (use specific mapping if available)
    ref_fasta = find_reference(sample_id, batch)
    if not ref_fasta:
        logger.warning("No reference for %s, skipping", sample_id)
        continue
    ref_fasta = Path(ref_fasta)

    work_dir = Path(WORK) / sample_id
    work_dir.mkdir(parents=True, exist_ok=True)
    sra_dir = work_dir / "sra_downloads"
    bam_dir = work_dir / "alignments"
    snif_dir = work_dir / "sniffles_output"
    circle_dir = work_dir / "circle_output"
    fmt_dir = work_dir / "is_formatter_output"

    try:
        # Download FASTQ
        sra_dir.mkdir(exist_ok=True)
        downloader = SRADownloader(str(sra_dir))
        downloader.download_one(sample_id)
        fq = find_fastq(sra_dir, sample_id)
        if not fq:
            logger.error("No FASTQ for %s", sample_id)
            n_err += 1
            continue

        # Align
        aligner = Aligner(output_dir=str(bam_dir), threads=8)
        bam = aligner.align(fq, ref_fasta, sample_id=sample_id)
        if not bam:
            logger.error("No BAM for %s", sample_id)
            n_err += 1
            continue

        # Sniffles
        sniff = SnifflesRunner(output_dir=str(snif_dir), alignment_dir=str(bam_dir))
        table = sniff.run_sample(bam=Path(bam), ref_fasta=ref_fasta, sample_id=sample_id)
        if not table:
            logger.error("No Sniffles for %s", sample_id)
            n_err += 1
            continue

        # Circle detection
        cf = CircleFinder(output_dir=str(circle_dir), threads=8)
        cf.run_sample(tldr_table=table, sample_id=sample_id, fastq_path=fq)

        # Formatter with require_th_reads=False
        formatter = ISFormatter(output_dir=str(fmt_dir), require_th_reads=False)
        formatter.run_sample(circle_dir=circle_dir / sample_id, sample_id=sample_id)

        # ORF annotation
        raw_json = fmt_dir / sample_id / f"{sample_id}_is_records.json"
        ann_json = fmt_dir / sample_id / f"{sample_id}_is_records_annotated.json"
        if raw_json.exists():
            ann = ORFAnnotator()
            ann.annotate_sample(str(raw_json))

        # Guide finder
        guide_json = fmt_dir / sample_id / f"{sample_id}_is_records_guide.json"
        if ann_json.exists():
            gf = GuideFinder(min_length=9, max_mismatches=1, min_length_for_mismatch=12)
            gf.find_guides_sample(str(ann_json))

        # Extract the IS621 record
        if guide_json.exists():
            with open(guide_json) as f:
                records = json.load(f)
            target = next((r for r in records if r["is_id"] == is_id), None)
            if target:
                target["organism"] = organism
                target["_batch"] = batch
                target["_source"] = "reprocessed"
                filt_records.append(target)

                # Visualize
                tag = f"{sample_id}__{is_id[:8]}"
                png_path = os.path.join(OUT, f"{tag}.png")
                gbk_path = os.path.join(OUT, f"{tag}.gbk")

                # Load partial circles for this is_id
                pc = []
                for pcf in glob.glob(f"{BASE}/batch_*/partial_circle_output/**/*_partial_circle_summary.json", recursive=True):
                    with open(pcf) as f:
                        data = json.load(f)
                    if isinstance(data, list):
                        for call in data:
                            if call.get("is_id") == is_id:
                                pc.append(call)

                vis.save_element_png(target, target.get("guide_hits", []), png_path,
                                     dpi=150, circle_info=target.get("circle_evidence"),
                                     partial_circles=pc if pc else None)
                gb.save_genbank(target, target.get("guide_hits", []), gbk_path,
                                circle_info=target.get("circle_evidence"))
                n_done += 1
                logger.info("Done %s: flanking_up=%d, flanking_down=%d, guide_hits=%d",
                            tag,
                            len(target.get("flanking_upstream", {}).get("sequence", "")),
                            len(target.get("flanking_downstream", {}).get("sequence", "")),
                            len(target.get("guide_hits", [])))

    except Exception as e:
        logger.error("Failed %s: %s", sample_id, e)
        n_err += 1
    finally:
        # Cleanup FASTQ + BAM to save space
        if sra_dir.exists():
            shutil.rmtree(sra_dir, ignore_errors=True)
        if bam_dir.exists():
            shutil.rmtree(bam_dir, ignore_errors=True)

logger.info("Reprocessed: %d done, %d errors", n_done, n_err)
with open(f"{OUT}/long_is621_records_flanked.json", "w") as f:
    json.dump(filt_records, f, indent=2)
logger.info("Done. Output: %s", OUT)
