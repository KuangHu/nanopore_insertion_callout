#!/usr/bin/env python3
"""Pick the fixed real-data panel used for before/after comparisons.

Reads data/batches/batch_*.tsv and writes data/panels/panel_v1.tsv.
The draw is seeded, so rerunning gives the same table.
"""
import glob
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SEED = 20261126
MIN_BASES = 2e8       # enough depth for a bacterial genome
MAX_BASES = 1.5e9     # keeps download and alignment cheap
MAX_PER_ORGANISM = 3  # per instrument, so no species dominates
QUOTA = {
    ("OXFORD_NANOPORE", "MinION"): 60,
    ("OXFORD_NANOPORE", "GridION"): 40,
    ("OXFORD_NANOPORE", "PromethION"): 20,
    ("PACBIO_SMRT", "PacBio RS"): 15,
    ("PACBIO_SMRT", "PacBio RS II"): 20,
    ("PACBIO_SMRT", "Sequel"): 15,
    ("PACBIO_SMRT", "Sequel II"): 15,
    ("PACBIO_SMRT", "Sequel IIe"): 8,
    ("PACBIO_SMRT", "Revio"): 7,
}

tables = []
for path in sorted(glob.glob(str(ROOT / "data/batches/batch_*.tsv"))):
    table = pd.read_csv(path, sep="\t")
    table["source_batch"] = Path(path).stem
    tables.append(table)
samples = pd.concat(tables)
samples = samples[samples.total_bases.between(MIN_BASES, MAX_BASES)]

picked = []
for (platform, instrument), n in QUOTA.items():
    group = samples[(samples.platform == platform) & (samples.instrument == instrument)]
    group = group.sample(frac=1, random_state=SEED)
    group = group[group.groupby("organism").cumcount() < MAX_PER_ORGANISM]
    picked.append(group.head(n))
panel = pd.concat(picked).sort_values("srr_accession")

out = ROOT / "data/panels/panel_v1.tsv"
panel.to_csv(out, sep="\t", index=False)
print(f"{len(panel)} samples -> {out}")
