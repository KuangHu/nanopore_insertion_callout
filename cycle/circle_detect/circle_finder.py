"""Detect IS circular intermediates and insertion boundaries via bait mapping.

For each IS consensus from tldr, a concatemer bait ``[IS][IS]`` is built.
Reads from circular intermediates map across the tail-head junction at the
center seam.  Genome-head and tail-genome junctions are detected via
soft-clip analysis on reads aligned near IS copy boundaries within the
same bait — no reference genome flanking sequence is needed.
"""

import csv
import logging
import math
import shutil
import subprocess
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import pysam

from ..processes import run_checked_pipe

from ..utils import find_fastq
from .config import (
    DEFAULT_BOUNDARY_TOLERANCE,
    DEFAULT_MAX_TARGET_DUPLICATION,
    DEFAULT_MIN_DUPLICATION_READS,
    DEFAULT_CIRCLE_OUTPUT_DIR,
    DEFAULT_JUNCTION_WINDOW,
    DEFAULT_MIN_CONSENSUS_ENTROPY,
    DEFAULT_MIN_CONSENSUS_LENGTH,
    DEFAULT_MIN_JUNCTION_OVERLAP,
    JUNCTION_FLANK,
    MAX_RECORDED_READS,
)

logger = logging.getLogger(__name__)

READ_COLUMNS = [
    "read_id", "sample_id", "is_uuid", "junction_type", "chrom",
    "start", "end", "family", "subfamily", "read_length",
    "alignment_start", "alignment_end", "junction_pos",
    "overlap_left", "overlap_right", "mapping_quality",
    "strand", "junction_read_pos", "junction_seq", "junction_deletion",
]
SUMMARY_COLUMNS = [
    "is_uuid", "chrom", "start", "end", "family", "subfamily",
    "consensus_length", "consensus", "target_duplication", "n_tail_head_reads",
    "n_genome_head_reads", "n_tail_genome_reads", "n_total_mapped",
    "example_th_read", "example_th_read_id", "junction_seq", "tail_head_read_ids",
]


@dataclass
class ISEntry:
    """A single IS element insertion from a tldr table."""

    uuid: str
    chrom: str
    start: int
    end: int
    family: str
    subfamily: str
    consensus: str


class CircleFinder:
    """Detect IS circular intermediates and insertion boundaries via bait mapping.

    For each sample, concatemer bait ``[IS][IS]`` sequences are built for all
    IS elements from tldr/sniffles.  Tail-head junctions are detected by reads
    spanning the center seam.  Genome-head and tail-genome junctions are
    detected via soft-clip analysis on reads aligned near IS copy boundaries.
    """

    def __init__(
        self,
        output_dir: str = DEFAULT_CIRCLE_OUTPUT_DIR,
        min_overlap: int = DEFAULT_MIN_JUNCTION_OVERLAP,
        min_consensus_length: int = DEFAULT_MIN_CONSENSUS_LENGTH,
        boundary_tolerance: int = DEFAULT_BOUNDARY_TOLERANCE,
        min_consensus_entropy: float = DEFAULT_MIN_CONSENSUS_ENTROPY,
        threads: int = 8,
        sort_memory: str = "4G",
        flank_length: int | None = None,  # deprecated, ignored
    ):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.min_overlap = min_overlap
        self.min_consensus_length = min_consensus_length
        self.boundary_tolerance = boundary_tolerance
        self.min_consensus_entropy = min_consensus_entropy
        self.threads = threads
        self.sort_memory = sort_memory

        for tool in ("minimap2", "samtools"):
            if not shutil.which(tool):
                raise RuntimeError(f"{tool} not found in PATH")

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _seq_entropy(seq: str) -> float:
        """Shannon entropy in bits per base (0 = homopolymer, 2 = random)."""
        if not seq:
            return 0.0
        seq = seq.upper()
        n = len(seq)
        counts = Counter(seq)
        return -sum((c / n) * math.log2(c / n) for c in counts.values())

    def _parse_tldr_table(self, table_path: Path) -> list[ISEntry]:
        """Parse a tldr .table.txt and return IS entries with valid consensus.

        Filters out entries whose consensus is empty, 'NA', shorter than
        ``min_consensus_length``, or low-complexity (entropy below
        ``min_consensus_entropy``).
        """
        table_path = Path(table_path)
        entries: list[ISEntry] = []
        skipped_low_complexity = 0

        with open(table_path) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                consensus = row.get("Consensus", "").strip()
                if not consensus or consensus == "NA":
                    continue
                if len(consensus) < self.min_consensus_length:
                    continue
                if self._seq_entropy(consensus) < self.min_consensus_entropy:
                    skipped_low_complexity += 1
                    continue

                entries.append(ISEntry(
                    uuid=row["UUID"],
                    chrom=row["Chrom"],
                    start=int(row["Start"]),
                    end=int(row["End"]),
                    family=row.get("Family", "NA"),
                    subfamily=row.get("Subfamily", "NA"),
                    consensus=consensus,
                ))

        logger.info(
            f"Parsed {len(entries)} IS entries with consensus >= "
            f"{self.min_consensus_length} bp from {table_path.name}"
            f" (skipped {skipped_low_complexity} low-complexity)"
        )
        return entries

    def _build_bait_fasta(
        self,
        entries: list[ISEntry],
        output_path: Path,
    ) -> tuple[Path, dict[str, list[dict]]]:
        """Write concatemer bait FASTA ``[IS][IS]`` for each IS element.

        Header format::

            >{uuid}__th__j{N}

        Returns:
            ``(fasta_path, junction_map)`` where *junction_map* maps
            ref_name to a list of ``{type, position, uuid}`` dicts.
        """
        output_path = Path(output_path)
        junction_map: dict[str, list[dict]] = {}
        lines = []
        for entry in entries:
            n = len(entry.consensus)
            th_name = f"{entry.uuid}__th__j{n}"
            lines.append(f">{th_name}\n{entry.consensus}{entry.consensus}\n")
            junction_map[th_name] = [
                {"type": "tail_head", "position": n, "uuid": entry.uuid},
            ]
        text = "".join(lines)

        if output_path.exists() and output_path.read_text() == text:
            logger.info(f"Bait FASTA exists: {output_path}")
            return output_path, junction_map

        # New or changed candidates: reads mapped to an older bait no longer apply.
        for old_bam in output_path.parent.glob("*.circle.sorted.bam*"):
            old_bam.unlink()
        output_path.write_text(text)

        logger.info(
            f"Wrote bait FASTA with {len(junction_map)} tail-head "
            f"sequences to {output_path}"
        )
        return output_path, junction_map

    def _parse_bait_headers(self, fasta_path: Path) -> dict[str, list[dict]]:
        """Reconstruct junction_map from existing bait FASTA headers.

        Handles current ``__th__`` format as well as the legacy ``__len{N}``
        format for backward compatibility.  Old ``__ic__`` headers are ignored.
        """
        junction_map: dict[str, list[dict]] = {}

        with open(fasta_path) as fh:
            for line in fh:
                if not line.startswith(">"):
                    continue
                ref_name = line[1:].strip().split()[0]

                if "__th__j" in ref_name:
                    parts = ref_name.split("__th__j")
                    uuid = parts[0]
                    junction_pos = int(parts[1])
                    junction_map[ref_name] = [
                        {"type": "tail_head", "position": junction_pos, "uuid": uuid},
                    ]

                elif "__len" in ref_name:
                    parts = ref_name.rsplit("__len", 1)
                    uuid = parts[0]
                    junction_pos = int(parts[1])
                    junction_map[ref_name] = [
                        {"type": "tail_head", "position": junction_pos, "uuid": uuid},
                    ]

        return junction_map

    def _map_reads(
        self, fastq: Path, bait_fasta: Path, output_bam: Path,
    ) -> Optional[Path]:
        """Map reads to the bait reference with minimap2.

        Uses ``--secondary=no`` so each read maps to its best bait target only.
        Pipes minimap2 directly into samtools sort (no intermediate SAM).
        """
        fastq = Path(fastq)
        output_bam = Path(output_bam)
        pending_bam = output_bam.with_name(f".{output_bam.stem}.pending.bam")

        if output_bam.exists() and Path(str(output_bam) + ".bai").exists():
            logger.info(f"BAM exists: {output_bam}")
            return output_bam

        mm2_cmd = [
            "minimap2",
            "-a",
            "--secondary=no",
            "-x", "map-ont",
            "-t", str(self.threads),
            str(bait_fasta),
            str(fastq),
        ]
        sort_cmd = [
            "samtools", "sort",
            "-@", str(self.threads),
            "-m", self.sort_memory,
            "-o", str(pending_bam),
            "-",
        ]

        logger.info(f"Mapping {fastq.name} -> {output_bam.name}")
        try:
            run_checked_pipe(mm2_cmd, sort_cmd)
            subprocess.run(
                ["samtools", "index", str(pending_bam)],
                capture_output=True, text=True, check=True,
            )
            pending_bam.replace(output_bam)
            Path(str(pending_bam) + ".bai").replace(Path(str(output_bam) + ".bai"))
        except Exception as exc:
            logger.error(f"Mapping pipeline failed: {exc}; {getattr(exc, 'stderr', '')}")
            return None
        finally:
            pending_bam.unlink(missing_ok=True)
            Path(str(pending_bam) + ".bai").unlink(missing_ok=True)

        logger.info(f"  -> {output_bam} + .bai")
        return output_bam

    def _find_junction_reads(
        self,
        bam_path: Path,
        entries_by_uuid: dict[str, ISEntry],
        junction_map: dict[str, list[dict]],
        sample_id: str,
        trim_duplication: bool = False,
    ) -> tuple[list[dict], dict[str, dict]]:
        """Scan BAM for reads spanning bait junctions.

        **Tail-head**: read alignment spans the junction at position N
        (center of ``[IS][IS]`` bait) by at least ``min_overlap`` on each side.

        **Genome-head**: read has a large left soft-clip (genomic) and its
        alignment starts near an IS copy boundary (0 or N on the TH bait).

        **Tail-genome**: read has a large right soft-clip (genomic) and its
        alignment ends near an IS copy boundary (N or 2N on the TH bait).

        Returns:
            ``(junction_reads, summary_by_uuid)`` where *junction_reads* is a
            list of per-read dicts, and *summary_by_uuid* maps uuid to counts
            split by junction type.
        """
        junction_reads: list[dict] = []
        total_mapped: dict[str, int] = {}
        counts: dict[str, dict[str, int]] = {}  # uuid -> {type -> count}
        # uuid -> (mapq, length, sequence, read id, junction sequence) for best tail-head example read
        th_example: dict[str, tuple[int, int, str, str, str]] = {}
        th_read_ids: dict[str, list[str]] = {}
        junction_gaps: dict[str, list[tuple[int, int]]] = {}  # uuid -> (start, length) per tail-head read

        tol = self.boundary_tolerance
        min_overlap = self.min_overlap

        bam = pysam.AlignmentFile(str(bam_path), "rb")
        for read in bam.fetch():
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue

            ref_name = read.reference_name
            if not ref_name or ref_name not in junction_map:
                continue

            junctions = junction_map[ref_name]
            uuid = junctions[0]["uuid"]

            entry = entries_by_uuid.get(uuid)
            if entry is None:
                continue

            total_mapped[uuid] = total_mapped.get(uuid, 0) + 1

            ref_start = read.reference_start   # 0-based
            ref_end = read.reference_end       # 0-based, exclusive
            aligned_len = ref_end - ref_start

            # Soft-clip lengths from CIGAR
            cigar = read.cigartuples
            left_clip = cigar[0][1] if cigar and cigar[0][0] == 4 else 0
            right_clip = cigar[-1][1] if cigar and cigar[-1][0] == 4 else 0

            # The junction_map for a TH bait has one entry with position = N
            junction_pos = junctions[0]["position"]
            n = junction_pos  # IS consensus length

            detected_types: list[str] = []

            # Count query-aligned anchors and require the junction to be crossed,
            # either gap-free or with one short deletion (a target site duplication).
            gap = self._junction_gap(cigar, ref_start, junction_pos, min_overlap)
            if gap is not None:
                detected_types.append("tail_head")
                junction_gaps.setdefault(uuid, []).append(gap)

            # Genome-head: left soft-clip is genomic, alignment starts
            # near an IS copy boundary (position 0 or N)
            if (left_clip >= min_overlap
                    and aligned_len >= min_overlap
                    and (ref_start < tol or abs(ref_start - n) < tol)):
                detected_types.append("genome_head")

            # Tail-genome: right soft-clip is genomic, alignment ends
            # near an IS copy boundary (position N or 2N)
            if (right_clip >= min_overlap
                    and aligned_len >= min_overlap
                    and (abs(ref_end - n) < tol or abs(ref_end - 2 * n) < tol)):
                detected_types.append("tail_genome")

            # Primary records are never hard clipped, so query positions are positions in the read.
            sequence = read.query_sequence or ""
            for junction_type in detected_types:
                if uuid not in counts:
                    counts[uuid] = {}
                counts[uuid][junction_type] = (
                    counts[uuid].get(junction_type, 0) + 1
                )

                # Where the junction falls in the read as aligned, and the bases around it.
                if junction_type == "tail_head":
                    query_pos = self._query_pos(cigar, ref_start, junction_pos)
                elif junction_type == "genome_head":
                    query_pos = left_clip
                else:
                    query_pos = read.query_length - right_clip
                junction_seq = sequence[max(0, query_pos - JUNCTION_FLANK):query_pos + JUNCTION_FLANK]

                if junction_type == "tail_head":
                    th_read_ids.setdefault(uuid, []).append(read.query_name)
                    if sequence:
                        candidate = (
                            read.mapping_quality,
                            len(sequence),
                            sequence,
                            read.query_name,
                            junction_seq,
                        )
                        if uuid not in th_example or candidate[:2] > th_example[uuid][:2]:
                            th_example[uuid] = candidate

                junction_reads.append({
                    "read_id": read.query_name,
                    "sample_id": sample_id,
                    "is_uuid": uuid,
                    "junction_type": junction_type,
                    "chrom": entry.chrom,
                    "start": entry.start,
                    "end": entry.end,
                    "family": entry.family,
                    "subfamily": entry.subfamily,
                    "read_length": read.query_length,
                    "strand": "-" if read.is_reverse else "+",
                    # Position in the read as sequenced, so it can be found again in the FASTQ.
                    "junction_read_pos": read.query_length - query_pos if read.is_reverse else query_pos,
                    "junction_seq": junction_seq,
                    "junction_deletion": gap[1] if junction_type == "tail_head" else "",
                    "alignment_start": ref_start,
                    "alignment_end": ref_end,
                    "junction_pos": junction_pos,
                    "overlap_left": junction_pos - ref_start,
                    "overlap_right": ref_end - junction_pos,
                    "mapping_quality": read.mapping_quality,
                })

        bam.close()

        # Build summary dict
        summary_by_uuid: dict[str, dict] = {}
        for uuid, entry in entries_by_uuid.items():
            uuid_counts = counts.get(uuid, {})
            # The insertion caller cannot tell the element from its target site duplication: it
            # reports the two together, so the bait junction holds the duplicated bases and reads
            # of the excised circle show them as a deletion. When most tail-head reads agree on
            # one deletion, record its length; with ``trim_duplication`` also cut the consensus
            # back to the element and move the position with it.
            consensus = entry.consensus
            shift = 0
            duplication = ""
            gaps = junction_gaps.get(uuid, [])
            if gaps:
                duplication, votes = Counter(length for _, length in gaps).most_common(1)[0]
                if votes < DEFAULT_MIN_DUPLICATION_READS or votes * 2 <= len(gaps):
                    duplication = ""
                elif duplication and trim_duplication:
                    starts = Counter(start for start, length in gaps if length == duplication and start is not None)
                    if starts:
                        gap_start = starts.most_common(1)[0][0]
                        n = len(entry.consensus)
                        shift = gap_start + duplication - n
                        consensus = (entry.consensus * 3)[gap_start + duplication:gap_start + n]
            summary_by_uuid[uuid] = {
                "is_uuid": uuid,
                "chrom": entry.chrom,
                "start": entry.start + shift,
                # A Sniffles insertion is a point (start == end) and stays one; a table that
                # gives the element as a span has its end follow the trimmed length.
                "end": entry.end + shift + (len(consensus) - len(entry.consensus) if entry.end != entry.start else 0),
                "family": entry.family,
                "subfamily": entry.subfamily,
                "consensus_length": len(consensus),
                "consensus": consensus,
                "target_duplication": duplication,
                "n_tail_head_reads": uuid_counts.get("tail_head", 0),
                "n_genome_head_reads": uuid_counts.get("genome_head", 0),
                "n_tail_genome_reads": uuid_counts.get("tail_genome", 0),
                "n_total_mapped": total_mapped.get(uuid, 0),
                "example_th_read": th_example[uuid][2] if uuid in th_example else "",
                "example_th_read_id": th_example[uuid][3] if uuid in th_example else "",
                "junction_seq": th_example[uuid][4] if uuid in th_example else "",
                "tail_head_read_ids": ",".join(th_read_ids.get(uuid, [])[:MAX_RECORDED_READS]),
            }

        return junction_reads, summary_by_uuid

    @staticmethod
    def _query_pos(cigar, ref_start, ref_pos):
        """Query position aligned to ``ref_pos``, or the next one when it falls in a deletion."""
        query = 0
        ref = ref_start
        for op, length in cigar or []:
            if op in (0, 7, 8):
                if ref + length > ref_pos:
                    return query + ref_pos - ref
                query += length
                ref += length
            elif op in (1, 4):
                query += length
            elif op in (2, 3):
                if ref + length > ref_pos:
                    return query
                ref += length
        return query

    @staticmethod
    def _junction_gap(cigar, ref_start, junction_pos, min_overlap):
        """Return the deletion a read has at the bait junction, or None if it is not junction evidence.

        The result is ``(start, length)``. ``length`` is the number of deleted bases within
        DEFAULT_JUNCTION_WINDOW of the junction, at most DEFAULT_MAX_TARGET_DUPLICATION;
        ``start`` is where they begin when they form one deletion, else None. A read aligned
        straight across gives ``(junction_pos, 0)``. Both sides need ``min_overlap`` aligned
        bases; indels away from the junction may occur in noisy reads but add nothing to
        either anchor. The aligner may place the deletion a base or two off the junction where
        neighbouring bases are equal, which is why a window is used.
        """
        low = junction_pos - DEFAULT_JUNCTION_WINDOW
        high = junction_pos + DEFAULT_JUNCTION_WINDOW
        ref_pos = ref_start
        block_start = None
        left_bases = right_bases = 0
        bridged = False
        deletions = []
        for op, length in cigar or []:
            if op in (0, 7, 8):  # consumes both reference and query
                if block_start is None:
                    block_start = ref_pos
                end = ref_pos + length
                left_bases += max(0, min(end, junction_pos) - ref_pos)
                right_bases += max(0, end - max(ref_pos, junction_pos))
                bridged |= block_start < junction_pos < end
                ref_pos = end
            else:
                if op in (2, 3):  # deletion or reference skip
                    if ref_pos <= high and ref_pos + length >= low:
                        if op == 3:
                            return None
                        deletions.append((ref_pos, length))
                    ref_pos += length
                block_start = None
        deleted = sum(length for _, length in deletions)
        if left_bases < min_overlap or right_bases < min_overlap or deleted > DEFAULT_MAX_TARGET_DUPLICATION:
            return None
        if deletions:
            return (deletions[0][0] if len(deletions) == 1 else None, deleted)
        return (junction_pos, 0) if bridged else None

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def run_sample(
        self,
        tldr_table: Path,
        sample_id: str,
        fastq_path: Path,
        ref_fasta: Path | None = None,
        trim_duplication: bool = False,
    ) -> Optional[Path]:
        """Run circle detection for one sample.

        Args:
            tldr_table: Path to the tldr/sniffles .table.txt for this sample.
            sample_id: Sample accession.
            fastq_path: Path to the sample's FASTQ file.
            ref_fasta: Accepted for backward compatibility but ignored.
            trim_duplication: Cut a target site duplication seen in the junction
                reads out of the consensus. Only safe when the consensus is accurate.

        Returns:
            Path to the sample output directory, or ``None`` on failure.
        """
        sample_dir = self.output_dir / sample_id
        sample_dir.mkdir(parents=True, exist_ok=True)

        # 1. Parse tldr table
        entries = self._parse_tldr_table(tldr_table)
        if not entries:
            # A sample without insertions is finished, not failed: leave empty tables.
            logger.warning(f"No valid IS entries for {sample_id}, nothing to detect")
            _write_tsv(sample_dir / f"{sample_id}_circle_reads.tsv", [], READ_COLUMNS)
            _write_tsv(sample_dir / f"{sample_id}_circle_summary.tsv", [], SUMMARY_COLUMNS)
            return sample_dir

        entries_by_uuid = {e.uuid: e for e in entries}

        # 2. Build bait FASTA (concatemer)
        bait_path = sample_dir / f"{sample_id}_bait.fa"
        bait_path, junction_map = self._build_bait_fasta(entries, bait_path)

        # 3. Map FASTQ and collect junction reads
        bam_path = sample_dir / f"{sample_id}.circle.sorted.bam"
        bam = self._map_reads(fastq_path, bait_path, bam_path)
        if bam is None:
            logger.warning(f"Mapping failed for {sample_id}, skipping")
            return None

        junction_reads, summary_by_uuid = self._find_junction_reads(
            bam, entries_by_uuid, junction_map, sample_id, trim_duplication,
        )

        logger.info(
            f"  {sample_id}: {len(junction_reads)} junction reads found"
        )

        # 4. Write output TSVs
        reads_tsv = sample_dir / f"{sample_id}_circle_reads.tsv"
        _write_tsv(reads_tsv, junction_reads, READ_COLUMNS)

        summary_tsv = sample_dir / f"{sample_id}_circle_summary.tsv"
        summary_rows = list(summary_by_uuid.values())
        _write_tsv(summary_tsv, summary_rows, SUMMARY_COLUMNS)

        total_junctions = sum(
            s["n_tail_head_reads"] + s["n_genome_head_reads"]
            + s["n_tail_genome_reads"]
            for s in summary_by_uuid.values()
        )
        logger.info(
            f"Circle detection for {sample_id}: {total_junctions} junction "
            f"reads across {len(entries)} IS elements"
        )
        return sample_dir

    def run_batch(
        self,
        tldr_results: dict[str, Optional[Path]],
        metadata,
        fastq_dir: str | Path,
        accession_col: str = "srr_accession",
        ref_map: dict[str, Path] | None = None,
        parallel: int = 1,
    ) -> dict[str, Optional[Path]]:
        """Run circle detection for all samples.

        Args:
            tldr_results: sample_id -> tldr table path mapping.
            metadata: DataFrame with accession column.
            fastq_dir: Directory containing FASTQ files.
            accession_col: Column name for SRR accession.
            ref_map: Accepted for backward compatibility but ignored.
            parallel: Number of samples to process in parallel.

        Returns:
            Dict mapping sample_id -> circle output directory (or None).
        """
        fastq_dir = Path(fastq_dir)
        results: dict[str, Optional[Path]] = {}

        # Collect tasks
        tasks: list[tuple[str, Path, Path, Path | None, bool]] = []
        for _, row in metadata.iterrows():
            sid = row[accession_col]
            table_path = tldr_results.get(sid)
            if not table_path:
                logger.warning(
                    f"No insertion table for {sid}, skipping circle detection"
                )
                results[sid] = None
                continue

            fq = find_fastq(fastq_dir, sid)
            if not fq:
                logger.warning(f"No FASTQ found for {sid}")
                results[sid] = None
                continue

            ref_fasta = None
            if ref_map:
                ref_fasta = ref_map.get(sid)

            # Nanopore consensus sequences are too rough to cut by a few bases.
            trim = row.get("platform") == "PACBIO_SMRT"
            tasks.append((sid, fq, table_path, ref_fasta, trim))

        logger.info(f"Running circle detection for {len(tasks)} samples")

        if parallel <= 1:
            for sid, fq, table_path, ref_fasta, trim in tasks:
                try:
                    out = self.run_sample(
                        tldr_table=table_path,
                        sample_id=sid,
                        fastq_path=fq,
                        ref_fasta=ref_fasta,
                        trim_duplication=trim,
                    )
                    results[sid] = out
                except Exception as exc:
                    logger.error(
                        f"Circle detection failed for {sid}, skipping: {exc}"
                    )
                    results[sid] = None
        else:
            logger.info(f"Running {len(tasks)} samples with {parallel} in parallel")
            with ProcessPoolExecutor(max_workers=parallel) as pool:
                futures = {}
                for sid, fq, table_path, ref_fasta, trim in tasks:
                    fut = pool.submit(
                        _run_circle_worker,
                        sample_id=sid,
                        fastq_path=fq,
                        table_path=table_path,
                        ref_fasta=ref_fasta,
                        trim_duplication=trim,
                        output_dir=self.output_dir,
                        threads=self.threads,
                        sort_memory=self.sort_memory,
                        min_overlap=self.min_overlap,
                        min_consensus_length=self.min_consensus_length,
                        boundary_tolerance=self.boundary_tolerance,
                        min_consensus_entropy=self.min_consensus_entropy,
                    )
                    futures[fut] = sid

                for fut in as_completed(futures):
                    sid = futures[fut]
                    try:
                        out = fut.result()
                        results[sid] = out
                        if out:
                            logger.info(f"  -> {out}")
                    except Exception as exc:
                        logger.error(
                            f"Circle detection worker failed for {sid}: {exc}"
                        )
                        results[sid] = None

        ok = sum(1 for v in results.values() if v)
        logger.info(
            f"Circle detection complete: {ok}/{len(results)} samples processed"
        )
        return results


def _run_circle_worker(
    sample_id: str,
    fastq_path: Path,
    table_path: Path,
    ref_fasta: Path | None,
    output_dir: Path,
    threads: int,
    sort_memory: str,
    min_overlap: int,
    min_consensus_length: int,
    boundary_tolerance: int = DEFAULT_BOUNDARY_TOLERANCE,
    min_consensus_entropy: float = DEFAULT_MIN_CONSENSUS_ENTROPY,
    flank_length: int | None = None,  # deprecated, ignored
    trim_duplication: bool = False,
) -> Optional[Path]:
    """Standalone worker function for parallel circle detection."""
    finder = CircleFinder(
        output_dir=str(output_dir),
        threads=threads,
        sort_memory=sort_memory,
        min_overlap=min_overlap,
        min_consensus_length=min_consensus_length,
        boundary_tolerance=boundary_tolerance,
        min_consensus_entropy=min_consensus_entropy,
    )
    return finder.run_sample(
        tldr_table=table_path,
        sample_id=sample_id,
        fastq_path=fastq_path,
        ref_fasta=ref_fasta,
        trim_duplication=trim_duplication,
    )


def _write_tsv(path: Path, rows: list[dict], columns: list[str]) -> None:
    """Write a list of dicts to a TSV file with the given column order."""
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=columns, delimiter="\t", extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)
    logger.info(f"Wrote {len(rows)} rows to {path}")
