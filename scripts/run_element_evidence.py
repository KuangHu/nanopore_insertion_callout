#!/usr/bin/env python3
"""Add ISfinder hits, transposase domains and orientation to annotated records.

Runs after run_orf_annotator.py and before run_guide_finder.py.

Usage:
    python scripts/run_element_evidence.py --input-dir is_formatter_output \
        --isfinder-fasta IS.fna --transposase-hmm Pfam_transposase.hmm \
        --db-dir isfinder_db [--sample-list samples.tsv] [--parallel 8]
"""

from __future__ import annotations

import argparse
import csv
import logging
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from cycle.element_evidence import annotate_sample, build_isfinder_db

logger = logging.getLogger(__name__)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--isfinder-fasta", required=True)
    parser.add_argument("--transposase-hmm", required=True)
    parser.add_argument("--db-dir", required=True, help="where the ISfinder BLAST database is kept")
    parser.add_argument("--sample-list")
    parser.add_argument("--parallel", type=int, default=1)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")

    db = build_isfinder_db(Path(args.isfinder_fasta), Path(args.db_dir))
    paths = sorted(Path(args.input_dir).glob("*/*_is_records_annotated.json"))
    if args.sample_list:
        with open(args.sample_list, newline="") as handle:
            wanted = {row["srr_accession"] for row in csv.DictReader(handle, delimiter="\t")}
        paths = [p for p in paths if p.parent.name in wanted]

    failed = []
    with ProcessPoolExecutor(max_workers=args.parallel) as pool:
        futures = {pool.submit(annotate_sample, p, db, Path(args.transposase_hmm)): p for p in paths}
        for future in as_completed(futures):
            try:
                future.result()
            except Exception as exc:
                logger.error(f"Element evidence failed for {futures[future].parent.name}: {exc}")
                failed.append(futures[future].parent.name)
    logger.info(f"Done: {len(paths) - len(failed)}/{len(paths)} samples")
    if paths and len(failed) == len(paths):
        sys.exit(1)


if __name__ == "__main__":
    main()
