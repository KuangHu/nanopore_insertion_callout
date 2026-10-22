#!/usr/bin/env python3
"""Run partial-circle jobs, retrying incomplete samples and reporting failures."""

import argparse
import csv
import logging
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from cycle.utils import find_fastq

logger = logging.getLogger(__name__)


def run_batch(input_dir, fastq_dir, output_dir, sample_list=None, threads=8):
    input_dir, fastq_dir, output_dir = map(Path, (input_dir, fastq_dir, output_dir))
    if sample_list:
        with open(sample_list, newline="") as fh:
            samples = [row["srr_accession"] for row in csv.DictReader(fh, delimiter="\t")]
    else:
        samples = [p.parent.name for p in sorted(input_dir.glob("*/*_is_records_guide.json"))]
    failures = []
    runner = Path(__file__).with_name("run_partial_circle.py")
    for sid in dict.fromkeys(samples):
        records = input_dir / sid / f"{sid}_is_records_guide.json"
        if not records.is_file():
            logger.info("No guide records for %s, skipping", sid)
            continue
        try:
            fastq = find_fastq(fastq_dir, sid)
            if fastq is None:
                raise FileNotFoundError(f"No FASTQ for {sid}")
            # The single-sample runner verifies completion against the inputs
            # and writes output_dir/sample_id itself (avoid double nesting).
            subprocess.run([
                sys.executable, str(runner), "--is-records", str(records),
                "--fastq", str(fastq), "--sample-id", sid,
                "--output-dir", str(output_dir), "--threads", str(threads),
            ], check=True)
        except (OSError, ValueError, subprocess.CalledProcessError) as exc:
            logger.error("Partial circle failed for %s: %s", sid, exc)
            failures.append(sid)
    return failures


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--fastq-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--sample-list")
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    failures = run_batch(**vars(args))
    if failures:
        logger.error("Partial circle incomplete: %d failed samples", len(failures))
        sys.exit(1)


if __name__ == "__main__":
    main()
