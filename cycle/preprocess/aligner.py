"""Align long reads to reference genomes with minimap2 → sorted BAM."""

import gzip
import logging
import math
import os
import statistics
import shutil
import subprocess
from pathlib import Path
from typing import Optional

import pysam

from ..processes import run_checked_pipe
from .config import (
    DEFAULT_ALIGNMENT_DIR,
    DEFAULT_MINIMAP2_PRESET,
    DEFAULT_SORT_MEMORY,
    DEFAULT_THREADS,
)

logger = logging.getLogger(__name__)

# PacBio reads whose median per-read accuracy reaches Q20 are HiFi; below that, or without base
# qualities (one value repeated over the whole read), they are CLR.
HIFI_MIN_READ_Q = 20
PRESET_READS = 20000


def pick_preset(fastq: Path, platform: str) -> str:
    """minimap2 preset for one run: map-ont for nanopore, map-hifi or map-pb for PacBio by read quality."""
    if platform != "PACBIO_SMRT":
        return DEFAULT_MINIMAP2_PRESET
    accuracies, flat = [], 0
    opener = gzip.open if str(fastq).endswith(".gz") else open
    with opener(fastq, "rt") as handle:
        for n, line in enumerate(handle):
            if n % 4 == 3:
                quality = line.rstrip("\n")
                if len(set(quality)) <= 1:
                    flat += 1
                error = sum(10 ** (-(ord(c) - 33) / 10) for c in quality) / max(len(quality), 1)
                accuracies.append(-10 * math.log10(max(error, 1e-6)))
                if len(accuracies) >= PRESET_READS:
                    break
    if not accuracies or flat > len(accuracies) / 2 or statistics.median(accuracies) < HIFI_MIN_READ_Q:
        return "map-pb"
    return "map-hifi"


def zmw_key(header: str) -> str:
    """The ZMW a read came from: ``movie/zmw`` of a PacBio name, else the SRA spot (first word).

    SRA keeps the subreads of one ZMW as the reads of one spot, all named ``<run>.<spot>``.
    """
    words = header[1:].split()
    for word in words:
        fields = word.split("/")
        if len(fields) >= 3 and fields[0].startswith("m") and fields[1].isdigit():
            return f"{fields[0]}/{fields[1]}"
    return words[0] if words else ""


def longest_per_zmw(fastq: Path, out: Path) -> tuple[int, int]:
    """Write the longest read of each ZMW to ``out`` (plain FASTQ); returns (reads in, reads out).

    CLR runs hold several subreads per ZMW, passes over one molecule, which pile up as
    separate reads at one site and split one insertion into several calls. The subreads
    of a ZMW are consecutive in the file.
    """
    opener = gzip.open if str(fastq).endswith(".gz") else open
    n_in = n_out = 0
    best, best_key = None, None
    pending = out.with_name(f".{out.name}.pending")
    with opener(fastq, "rt") as handle, open(pending, "w") as sink:
        while True:
            record = [handle.readline() for _ in range(4)]
            key = zmw_key(record[0]) if record[0] else None
            if key != best_key or not record[0]:
                if best is not None:
                    sink.writelines(best)
                    n_out += 1
                best, best_key = None, key
            if not record[0]:
                break
            n_in += 1
            if best is None or len(record[1]) > len(best[1]):
                best = record
    pending.replace(out)
    return n_in, n_out


class Aligner:
    """Minimap2 index + align → sorted BAM + .bai.

    One ``.mmi`` index is built per reference genome and preset (shared across
    samples of the same organism), since an index keeps the k-mer settings of
    the preset it was built with. With no preset given, each sample gets one
    from its platform and read quality (``pick_preset``). CLR runs (map-pb) are
    aligned with the longest subread of each ZMW only (``longest_per_zmw``). Alignment pipes
    minimap2 directly into ``samtools sort`` with no intermediate SAM on disk.
    """

    # Default per-sample alignment timeout in seconds (2 hours)
    DEFAULT_TIMEOUT = 7200

    def __init__(
        self,
        output_dir: str = DEFAULT_ALIGNMENT_DIR,
        preset: Optional[str] = None,
        threads: int = DEFAULT_THREADS,
        sort_memory: str = DEFAULT_SORT_MEMORY,
        timeout: int = DEFAULT_TIMEOUT,
    ):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.preset = preset
        self.threads = threads
        self.sort_memory = sort_memory
        self.timeout = timeout

        for tool in ("minimap2", "samtools"):
            if not shutil.which(tool):
                raise RuntimeError(f"{tool} not found in PATH")

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def index(self, fasta: Path, preset: str) -> Optional[Path]:
        """Build the minimap2 index of a reference FASTA for one preset, ``<fasta stem>.<preset>.mmi``.

        Skips if the index already exists. The index is written under a temporary name and
        renamed, so a half-written index is never used. Returns its path, or None on failure.
        """
        fasta = Path(fasta)
        mmi = fasta.with_suffix(f".{preset}.mmi")

        if mmi.exists():
            logger.info(f"Index exists: {mmi}")
            return mmi

        pending = mmi.with_name(f".{mmi.name}.{os.getpid()}.pending")
        cmd = [
            "minimap2",
            "-x", preset,
            "-t", str(self.threads),
            "-d", str(pending),
            str(fasta),
        ]
        logger.info(f"Indexing {fasta.name} for {preset}")
        ret = subprocess.run(cmd, capture_output=True, text=True)
        if ret.returncode != 0:
            logger.error(f"minimap2 index failed: {ret.stderr.strip()}")
            pending.unlink(missing_ok=True)
            return None
        if not pending.exists():
            logger.error(f"minimap2 left no index for {fasta.name}")
            return None
        pending.replace(mmi)
        logger.info(f"  → {mmi}")
        return mmi

    def align(
        self,
        fastq: Path,
        reference: Path,
        sample_id: Optional[str] = None,
        preset: Optional[str] = None,
    ) -> Optional[Path]:
        """Align reads to a reference and produce a sorted BAM + index.

        Args:
            fastq: Path to input FASTQ (or .fastq.gz).
            reference: Path to reference FASTA.
            sample_id: Used for output filename. Defaults to FASTQ stem.
            preset: minimap2 preset for this sample; defaults to the aligner's, then map-ont.
                The index ``<fasta stem>.<preset>.mmi`` is used when it exists.

        Returns:
            Path to the sorted BAM, or None on failure.
        """
        fastq = Path(fastq)
        reference = Path(reference)
        if sample_id is None:
            sample_id = fastq.name.split(".")[0]

        bam = self.output_dir / f"{sample_id}.sorted.bam"
        pending_bam = bam.with_name(f".{bam.stem}.pending.bam")

        preset = preset or self.preset or DEFAULT_MINIMAP2_PRESET
        mmi = reference.with_suffix(f".{preset}.mmi")
        ref = mmi if mmi.exists() else reference

        # Reuse a BAM only when its minimap2 command line names this reference, preset and read file
        if bam.exists() and Path(str(bam) + ".bai").exists():
            try:
                with pysam.AlignmentFile(str(bam), "rb", check_sq=False) as handle:
                    commands = [pg.get("CL", "") for pg in handle.header.to_dict().get("PG", []) if pg.get("PN") == "minimap2"]
            except (OSError, ValueError):
                commands = []
            if any(str(ref) in cl.split() and f"-x {preset}" in cl and fastq.name in cl for cl in commands):
                logger.info(f"BAM exists: {bam}")
                return bam
            logger.info(f"BAM {bam} was made with another reference or preset; realigning")

        # minimap2 -a | samtools sort → BAM (pipe-based, no intermediate SAM)
        # --MD flag is required for Sniffles2 to parse alignments correctly
        mm2_cmd = [
            "minimap2",
            "-a",
            "--MD",
            "-x", preset,
            "-t", str(self.threads),
            str(ref),
            str(fastq),
        ]
        sort_cmd = [
            "samtools", "sort",
            "-@", str(self.threads),
            "-m", self.sort_memory,
            "-o", str(pending_bam),
            "-",
        ]

        logger.info(f"Aligning {fastq.name} → {bam.name} with {preset}")
        try:
            run_checked_pipe(mm2_cmd, sort_cmd, timeout=self.timeout)
            subprocess.run(
                ["samtools", "index", str(pending_bam)],
                capture_output=True, text=True, check=True, timeout=self.timeout,
            )
            # Publish only a complete BAM and index. Interrupted pending files
            # are never eligible for the existing-output shortcut above.
            pending_bam.replace(bam)
            Path(str(pending_bam) + ".bai").replace(Path(str(bam) + ".bai"))
        except Exception as exc:
            logger.error(f"Alignment pipeline failed for {sample_id}: {exc}; {getattr(exc, 'stderr', '')}")
            return None
        finally:
            pending_bam.unlink(missing_ok=True)
            Path(str(pending_bam) + ".bai").unlink(missing_ok=True)

        logger.info(f"  → {bam} + .bai")
        return bam

    def align_batch(
        self,
        sample_map: list[dict],
    ) -> list[dict]:
        """Align a batch of samples.

        Args:
            sample_map: List of dicts with keys
                {sample_id, fastq, reference_fasta} and optionally platform.

        Returns:
            List of dicts with added ``bam`` key (Path or None).
        """
        results = []
        for i, entry in enumerate(sample_map, 1):
            sid = entry["sample_id"]
            fastq = Path(entry["fastq"])
            ref = Path(entry["reference_fasta"])

            logger.info(
                f"[{i}/{len(sample_map)}] {sid}"
            )

            preset = self.preset or pick_preset(fastq, entry.get("platform", ""))
            self.index(ref, preset)

            # CLR: one read per ZMW, the longest subread
            reads = fastq
            if preset == "map-pb":
                reads = self.output_dir / f"{sid}.longest_subread.fastq"
                n_in, n_out = longest_per_zmw(fastq, reads)
                logger.info(f"{sid}: kept the longest subread of each ZMW, {n_out}/{n_in} reads")
            try:
                bam = self.align(reads, ref, sample_id=sid, preset=preset)
            finally:
                if reads != fastq:
                    reads.unlink(missing_ok=True)
            results.append({**entry, "bam": bam})

        ok = sum(1 for r in results if r["bam"])
        logger.info(f"Alignment complete: {ok}/{len(results)} succeeded")
        return results
