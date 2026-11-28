"""Run Sniffles2 per sample to detect insertions as IS candidates."""

import logging
import math
import re
import statistics
import shutil
import subprocess
import uuid
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Optional

import pandas as pd
import pysam

from .config import DEFAULT_ALIGNMENT_DIR

logger = logging.getLogger(__name__)

MERGE_DISTANCE = 10  # bp; insertion calls closer than this are one site


class SnifflesRunner:
    """Run Sniffles2 once per sample to detect insertions.

    Sniffles2 is a structural variant caller that can detect insertions
    without requiring a transposon reference library. We use it as a
    faster alternative to tldr for finding IS insertion candidates.
    """

    def __init__(
        self,
        output_dir: str,
        alignment_dir: str = DEFAULT_ALIGNMENT_DIR,
    ):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.alignment_dir = Path(alignment_dir)

        if not shutil.which("sniffles"):
            raise RuntimeError("sniffles not found in PATH")

    @staticmethod
    def _seq_entropy(seq: str) -> float:
        """Shannon entropy in bits per base (0 = homopolymer, 2 = random)."""
        if not seq:
            return 0.0
        seq = seq.upper()
        n = len(seq)
        counts = Counter(seq)
        return -sum((c / n) * math.log2(c / n) for c in counts.values())

    def run_sample(
        self,
        bam: Path,
        ref_fasta: Path,
        sample_id: str,
        min_size: int = 500,
        max_size: int = 20000,
        min_support: int = 3,
        disable_qc: bool = False,
    ) -> Optional[Path]:
        """Run Sniffles2 on a single BAM for one sample.

        Args:
            bam: Sorted BAM file.
            ref_fasta: Reference genome FASTA used for alignment.
            sample_id: Sample accession (used for output directory naming).
            min_size: Minimum insertion size in bp.
            max_size: Maximum insertion size in bp.
            min_support: Minimum number of supporting reads.
            disable_qc: Disable Sniffles2 quality control filters for max sensitivity.

        Returns:
            Path to the Sniffles2-derived insertion table, or None on failure.
        """
        sample_dir = self.output_dir / sample_id
        sample_dir.mkdir(parents=True, exist_ok=True)

        output_prefix = sample_dir / sample_id
        table_path = Path(f"{output_prefix}.table.txt")
        vcf_path = Path(f"{output_prefix}.vcf")

        if table_path.exists():
            with open(table_path) as f:
                line_count = sum(1 for _ in f)
            if line_count > 1:
                logger.info(f"Sniffles2 output exists for {sample_id} ({line_count - 1} insertions): {table_path}")
                return table_path
            else:
                logger.info(f"Reprocessing {sample_id} — table exists but has no insertions")

        logger.info(f"Running Sniffles2 for {sample_id}: {bam.name}")

        # Skip if already processed with actual variants
        skip_existing = False
        if vcf_path.exists():
            with open(vcf_path) as f:
                variant_count = sum(1 for line in f if not line.startswith('#'))
            if variant_count > 0:
                logger.info(f"  {sample_id} - using existing VCF ({variant_count} variants)")
                skip_existing = True
            else:
                logger.info(f"  {sample_id} - reprocessing empty VCF")

        if not skip_existing:
            # Sniffles writes to a pending file, renamed only after a clean exit, so a VCF
            # under the final name always comes from a finished run.
            pending_vcf = sample_dir / f".{sample_id}.pending.vcf"
            cmd = [
                "sniffles",
                "--input", str(bam),
                "--vcf", str(pending_vcf),
                "--reference", str(ref_fasta),
                "--threads", "1",
                "--minsvlen", str(min_size),
                "--minsupport", str(min_support),
                "--output-rnames",
                "--allow-overwrite",
                # Sniffles 2.7 skips contigs shorter than 1 Mb (plasmids, small chromosomes) otherwise
                "--all-contigs",
            ]
            if disable_qc:
                cmd.append("--no-qc")

            try:
                ret = subprocess.run(
                    cmd, capture_output=True, text=True, timeout=600,
                    cwd=str(sample_dir),
                )
                if ret.returncode != 0:
                    err_msg = ret.stderr.strip() or ret.stdout.strip() or "(no output)"
                    logger.warning(f"Sniffles2 failed for {sample_id} (exit={ret.returncode}): {err_msg[:300]}")
                    pending_vcf.unlink(missing_ok=True)
                    return None
            except subprocess.TimeoutExpired:
                logger.warning(f"Sniffles2 timed out for {sample_id}")
                pending_vcf.unlink(missing_ok=True)
                return None
            pending_vcf.replace(vcf_path)

        # Parse VCF and collect insertions
        insertions = []
        if vcf_path.exists():
            try:
                insertions = self._parse_vcf(vcf_path, min_size, max_size, sample_id)
            except Exception as e:
                logger.warning(f"Failed to parse VCF for {sample_id}: {e}")
                return None

        # A symbolic ALT gives no bases: take them from the supporting reads
        symbolic = [ins for ins in insertions if not ins["Sequence"]]
        n_found = 0
        for ins in symbolic:
            # A duplication may be reported at either end of the copied stretch
            sites = [ins["Start"], ins["Start"] + ins["Length"]] if ins["Family"].endswith("DUP") else [ins["Start"]]
            for site in sites:
                inserted = self._inserted_from_reads(bam, ins["Chrom"], site, ins["Length"], ins["Rnames"])
                if inserted and self._seq_entropy(inserted) >= 1.5 and min_size <= len(inserted) <= max_size:
                    ins["Sequence"] = inserted
                    n_found += 1
                    break
        if symbolic:
            logger.info(f"  {sample_id}: inserted bases found in reads for {n_found}/{len(symbolic)} symbolic ALTs")

        # Remove entries that still have no usable sequence
        insertions = [ins for ins in insertions if ins["Sequence"]]

        # Several calls at one site (reads disagreeing on the inserted length) are one insertion
        n_called = len(insertions)
        insertions = self._merge_nearby(insertions)
        if len(insertions) < n_called:
            logger.info(f"  {sample_id}: {n_called - len(insertions)} calls merged into a call within {MERGE_DISTANCE} bp")

        # Insertions of a sequence found elsewhere in the reference, which Sniffles does not call
        recovered = []
        if Path(bam).exists():
            recovered = self._split_read_insertions(bam, insertions, min_size, max_size, min_support, sample_id)
        else:
            logger.warning(f"  {sample_id}: no BAM at {bam}, split reads not searched")
        if recovered:
            logger.info(f"  {sample_id}: {len(recovered)} insertions from split reads with the inserted part elsewhere")
        insertions += recovered

        # The supporting reads stay in the table: they tie reads to this locus downstream
        for ins in insertions:
            ins["Rnames"] = ",".join(ins["Rnames"])
            ins["Length"] = len(ins["Sequence"])

        if insertions:
            self._write_table(insertions, table_path)
            logger.info(f"  -> {table_path} ({len(insertions)} insertions)")
            return table_path
        else:
            logger.warning(f"No insertions found for {sample_id}")
            pd.DataFrame(columns=["UUID", "Chrom", "Start", "End", "Length", "Family", "Subfamily", "Consensus", "Support", "Rnames"]).to_csv(
                table_path, sep="\t", index=False
            )
            return table_path

    def _parse_vcf(
        self, vcf_path: Path, min_size: int, max_size: int, sample_id: str,
        min_entropy: float = 1.5,
    ) -> list[dict]:
        """Parse Sniffles2 VCF and extract insertions.

        An insertion is a point on the reference: ``Start`` and ``End`` are both
        the VCF position (the inserted bases follow that base) and the inserted
        length is in ``Length``. Symbolic ALTs (``<INS>``, ``<DUP>``) get
        ``Sequence=""`` and the supporting read names in ``Rnames``. Breakends
        are not insertions and are skipped. Resolved sequences with Shannon
        entropy below *min_entropy* are discarded.
        """
        insertions = []
        n_symbolic = 0
        n_low_entropy = 0
        n_total = 0

        with open(vcf_path) as f:
            for line in f:
                if line.startswith("#"):
                    continue

                fields = line.strip().split("\t")
                if len(fields) < 8:
                    continue

                chrom = fields[0]
                pos = int(fields[1])
                ref = fields[3]
                alt = fields[4]
                info = fields[7]

                # Parse INFO field
                info_dict = {}
                for item in info.split(";"):
                    if "=" in item:
                        key, val = item.split("=", 1)
                        info_dict[key] = val

                svtype = info_dict.get("SVTYPE", "")
                if svtype not in ("INS", "DUP"):
                    continue

                svlen = info_dict.get("SVLEN")
                if svlen:
                    svlen = abs(int(svlen))
                    if svlen < min_size or svlen > max_size:
                        continue
                else:
                    if alt.startswith("<"):
                        continue
                    svlen = len(alt) - len(ref)
                    if svlen < min_size or svlen > max_size:
                        continue

                n_total += 1
                support = int(info_dict.get("SUPPORT", "0"))

                # Extract read names (from --output-rnames)
                rnames_str = info_dict.get("RNAMES", "")
                rnames = [r for r in rnames_str.split(",") if r] if rnames_str else []

                if alt.startswith("<"):
                    # Symbolic ALT — mark for assembly instead of filling N's
                    n_symbolic += 1
                    sequence = ""
                else:
                    sequence = alt[len(ref):]
                    if self._seq_entropy(sequence) < min_entropy:
                        n_low_entropy += 1
                        continue

                insertions.append({
                    "UUID": str(uuid.uuid4()),
                    "Chrom": chrom,
                    "Start": pos,
                    "End": pos,
                    "Length": svlen,
                    "Sequence": sequence,
                    "Support": support,
                    "Family": f"Sniffles2_{svtype}",
                    "Subfamily": f"{sample_id}_{chrom}_{pos}",
                    "Rnames": rnames,
                })

        logger.info(
            f"  {sample_id} VCF: {n_total} calls, {n_symbolic} symbolic, "
            f"{n_low_entropy} low-entropy filtered, {len(insertions)} kept"
        )
        return insertions

    @staticmethod
    def _merge_nearby(insertions: list[dict]) -> list[dict]:
        """Keep one call per site: calls on a contig whose starts chain within MERGE_DISTANCE bp
        are one group, and the best-supported call of the group is kept (the first in VCF order
        on a tie)."""
        order = sorted(range(len(insertions)), key=lambda i: (insertions[i]["Chrom"], insertions[i]["Start"], i))
        kept, group = [], []
        for i in order + [None]:
            ins = insertions[i] if i is not None else None
            if group and (ins is None or ins["Chrom"] != insertions[group[-1]]["Chrom"]
                          or ins["Start"] - insertions[group[-1]]["Start"] > MERGE_DISTANCE):
                kept.append(max(group, key=lambda j: (insertions[j]["Support"], -j)))
                group = []
            if i is not None:
                group.append(i)
        return [insertions[i] for i in sorted(kept)]

    @staticmethod
    def _inserted_from_reads(
        bam_path: Path, chrom: str, pos: int, svlen: int, rnames: list[str], slack: int = 100,
    ) -> Optional[str]:
        """Inserted bases of a symbolic call, taken from one supporting read.

        A read shows the insertion either as one insertion inside an alignment
        or, when the aligner split the read, as the read bases between the piece
        that ends at the site and the piece that starts there. Of all supporting
        reads the one whose inserted length is closest to the called length is
        used. The bases are given on the reference strand. Returns None when no
        read shows the insertion.
        """
        wanted = set(rnames)
        by_read: dict[str, list] = {}
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            for read in bam.fetch(chrom, max(0, pos - slack), pos + slack):
                if read.query_name in wanted and not read.is_secondary and read.cigartuples:
                    by_read.setdefault(read.query_name, []).append(read)

        candidates = []
        for records in by_read.values():
            for read in records:
                if not read.query_sequence:
                    continue
                stored = 0  # position in the stored sequence, which leaves out hard clips
                ref = read.reference_start
                for op, length in read.cigartuples:
                    if op == 1:
                        if abs(ref - pos) <= slack and length * 2 >= svlen:
                            candidates.append(read.query_sequence[stored:stored + length])
                        stored += length
                    elif op in (0, 7, 8):
                        stored += length
                        ref += length
                    elif op == 4:
                        stored += length
                    elif op in (2, 3):
                        ref += length

            # Split read: only the primary record stores the whole read
            primary = next((r for r in records if not r.is_supplementary and r.query_sequence), None)
            if primary is None:
                continue
            spans = []  # (read start, read end, record), read positions counting every clip
            for read in records:
                if read.is_reverse != primary.is_reverse:
                    continue
                lead = 0
                for op, length in read.cigartuples:
                    if op not in (4, 5):
                        break
                    lead += length
                spans.append((lead, lead + read.query_alignment_length, read))
            for left_start, left_end, left in spans:
                for right_start, right_end, right in spans:
                    if (left is not right and abs(left.reference_end - pos) <= slack
                            and abs(right.reference_start - pos) <= slack and left_end < right_start):
                        candidates.append(primary.query_sequence[left_end:right_start])

        if not candidates:
            return None
        return min(candidates, key=lambda seq: abs(len(seq) - svlen))

    @staticmethod
    def _split_read_insertions(
        bam_path: Path, called: list[dict], min_size: int, max_size: int, min_reads: int, sample_id: str,
        join: int = 50, away: int = 2000, cluster: int = 30, near: int = 200, min_entropy: float = 1.5,
    ) -> list[dict]:
        """Insertions whose reads align the inserted part to another copy in the reference.

        A read across such an insertion aligns in three parts, in read order: the left flank,
        the inserted part on its other location (one or more pieces, at least *away* bp from the
        site), and the right flank, which meets the left flank on the reference within *join* bp
        on the same strand. Sniffles reports these reads as breakends at most, so an element
        inserted beside a resident copy of itself is missed. The read bases between the flanks
        are the insertion; at least half of them must sit on the other location. Sites within
        *cluster* bp are one insertion; one with at least *min_reads* reads and no call within
        *near* bp is returned in the format of ``_parse_vcf``, with the bases of the read whose
        inserted length is the median.
        """
        complement = str.maketrans("ACGTNacgtn", "TGCANtgcan")

        def piece(contig, start, reverse, cigar):
            lead = cigar[0][1] if cigar[0][0] in (4, 5) else 0
            aligned = sum(n for op, n in cigar if op in (0, 1, 7, 8))
            total = sum(n for op, n in cigar if op in (0, 1, 4, 5, 7, 8))
            ref_len = sum(n for op, n in cigar if op in (0, 2, 3, 7, 8))
            q0 = total - lead - aligned if reverse else lead
            return contig, reverse, q0, q0 + aligned, start, start + ref_len

        def cigar_ops(text):
            ops = re.findall(r"(\d+)([MIDNSHP=X])", text)
            return [("MIDNSHP=X".index(op), int(n)) for n, op in ops]

        sites: dict[str, list] = {}
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            for read in bam.fetch():
                if read.is_secondary or read.is_supplementary or not read.query_sequence or not read.has_tag("SA"):
                    continue
                pieces = [piece(read.reference_name, read.reference_start, read.is_reverse, read.cigartuples)]
                for entry in read.get_tag("SA").rstrip(";").split(";"):
                    contig, pos, strand, cigar = entry.split(",")[:4]
                    pieces.append(piece(contig, int(pos) - 1, strand == "-", cigar_ops(cigar)))
                pieces.sort(key=lambda x: x[2])
                sequence = read.query_sequence.translate(complement)[::-1] if read.is_reverse else read.query_sequence
                for i, left in enumerate(pieces):
                    for j in range(i + 2, len(pieces)):
                        right = pieces[j]
                        if right[0] != left[0] or right[1] != left[1]:
                            continue
                        # Along the read a reverse-strand alignment runs towards the reference start
                        site, other = (left[4], right[5]) if left[1] else (left[5], right[4])
                        if abs(other - site) > join:
                            continue
                        gap = right[2] - left[3]
                        elsewhere = [x for x in pieces[i + 1:j]
                                     if x[0] != left[0] or min(abs(x[4] - site), abs(x[5] - site)) > away]
                        if min_size <= gap <= max_size and sum(x[3] - x[2] for x in elsewhere) * 2 >= gap:
                            inserted = sequence[left[3]:right[2]]
                            if left[1]:
                                inserted = inserted.translate(complement)[::-1]
                            sites.setdefault(left[0], []).append((site, read.query_name, inserted))
                        break

        recovered = []
        for contig, items in sites.items():
            items.sort()
            groups = [[items[0]]]
            for item in items[1:]:
                if item[0] - groups[-1][-1][0] <= cluster:
                    groups[-1].append(item)
                else:
                    groups.append([item])
            for group in groups:
                reads = sorted({x[1] for x in group})
                pos = int(statistics.median(x[0] for x in group))
                if len(reads) < min_reads or any(c["Chrom"] == contig and abs(c["Start"] - pos) <= near for c in called):
                    continue
                length = statistics.median(len(x[2]) for x in group)
                inserted = min(group, key=lambda x: abs(len(x[2]) - length))[2]
                if SnifflesRunner._seq_entropy(inserted) < min_entropy:
                    continue
                recovered.append({
                    "UUID": str(uuid.uuid4()), "Chrom": contig, "Start": pos, "End": pos, "Length": len(inserted),
                    "Sequence": inserted, "Support": len(reads), "Family": "split_reads",
                    "Subfamily": f"{sample_id}_{contig}_{pos}", "Rnames": reads,
                })
        return recovered

    def _write_table(self, insertions: list[dict], output_path: Path):
        """Write insertion table in tldr-compatible format."""
        df = pd.DataFrame(insertions)
        df = df.rename(columns={"Sequence": "Consensus"})
        cols = ["UUID", "Chrom", "Start", "End", "Length", "Family", "Subfamily", "Consensus", "Support", "Rnames"]
        df = df[cols]
        df.to_csv(output_path, sep="\t", index=False)

    def run_batch(
        self,
        metadata: pd.DataFrame,
        ref_map: dict[str, Path],
        parallel: int = 1,
        min_size: int = 500,
        max_size: int = 20000,
        disable_qc: bool = False,
        accession_col: str = "srr_accession",
    ) -> dict[str, Optional[Path]]:
        """Run Sniffles2 for each sample.

        Args:
            metadata: DataFrame with accession column.
            ref_map: sample_id -> ref_fasta path.
            parallel: Number of samples to run in parallel.
            min_size: Minimum insertion size in bp.
            max_size: Maximum insertion size in bp.
            disable_qc: Disable Sniffles2 quality control filters.
            accession_col: Column name for SRR accession.

        Returns:
            Dict mapping sample_id -> insertion table path (or None).
        """
        results: dict[str, Optional[Path]] = {}

        # Build task list
        tasks: list[tuple[str, Path, Path]] = []
        for _, row in metadata.iterrows():
            sid = row[accession_col]
            ref_fasta = ref_map.get(sid)
            if not ref_fasta:
                logger.warning(f"No reference for {sid}, skipping Sniffles2")
                results[sid] = None
                continue

            bam = self.alignment_dir / f"{sid}.sorted.bam"
            if not bam.exists():
                logger.warning(f"BAM not found for {sid}: {bam}")
                results[sid] = None
                continue

            tasks.append((sid, bam, ref_fasta))

        logger.info(f"Running Sniffles2 for {len(tasks)} samples")

        if parallel <= 1:
            for sid, bam, ref_fasta in tasks:
                table = self.run_sample(
                    bam=bam, ref_fasta=ref_fasta,
                    sample_id=sid, min_size=min_size, max_size=max_size,
                    disable_qc=disable_qc,
                )
                results[sid] = table
        else:
            logger.info(f"Running {len(tasks)} samples with {parallel} in parallel")
            with ProcessPoolExecutor(max_workers=parallel) as pool:
                futures = {}
                for sid, bam, ref_fasta in tasks:
                    fut = pool.submit(
                        _run_sniffles_worker,
                        bam=bam,
                        ref_fasta=ref_fasta,
                        sample_id=sid,
                        min_size=min_size,
                        max_size=max_size,
                        disable_qc=disable_qc,
                        output_dir=self.output_dir,
                        alignment_dir=self.alignment_dir,
                    )
                    futures[fut] = sid

                for fut in as_completed(futures):
                    sid = futures[fut]
                    try:
                        table = fut.result()
                        results[sid] = table
                        if table:
                            logger.info(f"  -> {table}")
                    except Exception as e:
                        logger.error(f"Sniffles2 worker failed for {sid}: {e}")
                        results[sid] = None

        ok = sum(1 for v in results.values() if v)
        logger.info(f"Sniffles2 complete: {ok}/{len(results)} samples succeeded")
        return results


def _run_sniffles_worker(
    bam: Path,
    ref_fasta: Path,
    sample_id: str,
    min_size: int,
    max_size: int,
    disable_qc: bool,
    output_dir: Path,
    alignment_dir: Path,
) -> Optional[Path]:
    """Standalone worker function for parallel Sniffles2 execution."""
    runner = SnifflesRunner(
        output_dir=str(output_dir),
        alignment_dir=str(alignment_dir),
    )
    return runner.run_sample(
        bam=bam,
        ref_fasta=ref_fasta,
        sample_id=sample_id,
        min_size=min_size,
        max_size=max_size,
        disable_qc=disable_qc,
    )
