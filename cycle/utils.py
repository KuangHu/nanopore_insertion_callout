"""Shared utility functions used across Cycle subpackages."""

import re
from pathlib import Path


def slugify(name: str) -> str:
    """Convert organism name to a filesystem-safe slug."""
    slug = name.lower().strip()
    slug = re.sub(r"[^a-z0-9]+", "_", slug)
    return slug.strip("_")


def find_fastq(fastq_dir: str | Path, sample_id: str) -> Path | None:
    """Find FASTQ file for a sample, trying common naming patterns."""
    fastq_dir = Path(fastq_dir)
    for pattern in (
        f"{sample_id}.fastq.gz",
        f"{sample_id}.fastq",
        f"{sample_id}_1.fastq.gz",
        f"{sample_id}_pass.fastq.gz",
    ):
        path = fastq_dir / pattern
        if path.is_file():
            return path

    # Require a delimiter after the complete sample ID: SRR123 must never
    # select SRR1234. Exclude unfinished downloads and ambiguous shard sets.
    name_pattern = re.compile(
        re.escape(sample_id) + r"(?:[_.-][A-Za-z0-9_.-]+)?\.(?:fastq|fq)(?:\.gz)?$"
    )
    matches = sorted(
        p for p in fastq_dir.iterdir()
        if p.is_file() and name_pattern.fullmatch(p.name)
    ) if fastq_dir.is_dir() else []
    if len(matches) > 1:
        raise ValueError(f"Ambiguous FASTQ files for {sample_id}: {matches}")
    if len(matches) == 1:
        return matches[0]

    return None
