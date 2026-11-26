"""Round-1 review tests: partial circles (C5-C8), guide alignment (C9-C11), circle caches (C15).

Each test fails on the reviewed code because of one defect and should pass once
that defect is fixed. Inputs are synthetic; minimap2/samtools are replaced by
fakes that write small BAMs with pysam.
"""

import csv
import json
import random
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pysam

from cycle.circle_detect.circle_finder import CircleFinder
from cycle.guide_finder.alignment import ShortAlignmentFinder, reverse_complement
from cycle.partial_circle.detector import PartialCircleDetector

IS_ID = 'IS1'
SAMPLE = 'S1'
REF = ''.join(random.Random(1).choices('ACGT', k=1000))
CONSENSUS = ''.join(random.Random(2).choices('ACGT', k=400))


def segment(name, ref_start, cigar, seq, supplementary=False):
    """One forward-strand alignment record on the IS reference (tid 0)."""
    read = pysam.AlignedSegment()
    read.query_name = name
    read.flag = 2048 if supplementary else 0
    read.reference_id = 0
    read.reference_start = ref_start
    read.mapping_quality = 60
    read.cigarstring = cigar
    read.query_sequence = seq
    return read


def write_bam(path, records, ref_name=IS_ID, ref_length=len(REF), sort_order='queryname'):
    header = {'HD': {'VN': '1.6', 'SO': sort_order}, 'SQ': [{'SN': ref_name, 'LN': ref_length}]}
    with pysam.AlignmentFile(str(path), 'wb', header=header) as bam:
        for record in records:
            bam.write(record)


def back_jump_read(name, circle_start, circle_end, anchor=150):
    """A read that runs up to circle_end and continues at circle_start (soft clips only)."""
    seq = REF[circle_end - anchor:circle_end] + REF[circle_start:circle_start + anchor]
    return [
        segment(name, circle_end - anchor, f'{anchor}M{anchor}S', seq),
        segment(name, circle_start, f'{anchor}S{anchor}M', seq, supplementary=True),
    ]


class PartialCircleTest(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        self.records_json = self.root / 'is_records.json'
        self.records_json.write_text(json.dumps([{'is_id': IS_ID, 'is_element': {'sequence': REF}}]))
        self.fastq = self.root / f'{SAMPLE}.fastq'
        self.fastq.write_text('@r\nA\n+\nI\n')

    def run_detector(self, fake_mapper, **kwargs):
        """Run the public run_sample with the minimap2|samtools pipe replaced by fake_mapper.

        fake_mapper(mm2_cmd, output_bam) writes the name-sorted BAM the pipe would produce.
        Returns (per-read back-jump rows, summary calls).
        """
        def pipe(mm2_cmd, sort_cmd, *args, **kw):
            fake_mapper(mm2_cmd, sort_cmd[sort_cmd.index('-o') + 1])

        with patch('shutil.which', return_value='fake'):
            detector = PartialCircleDetector(output_dir=str(self.root / 'out'), **kwargs)
        with patch('cycle.partial_circle.detector.run_checked_pipe', side_effect=pipe):
            sample_dir = detector.run_sample(self.records_json, SAMPLE, self.fastq)
        self.assertIsNotNone(sample_dir, 'run_sample failed before any detection')
        with open(sample_dir / f'{SAMPLE}_partial_circle_reads.tsv') as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        calls = json.loads((sample_dir / f'{SAMPLE}_partial_circle_summary.json').read_text())
        return rows, calls

    def run_on_records(self, records, **kwargs):
        return self.run_detector(lambda mm2_cmd, out: write_bam(out, records), **kwargs)

    def test_c5_hard_clipped_supplementary_forward_read_is_not_a_circle(self):
        """A read whose pieces map in forward order (100, 400, 700) gives no back-jump."""
        read = REF[100:200] + REF[400:500] + REF[700:800]

        # Design assumption: the fix either puts hard-clipped records in the full-read
        # frame or asks minimap2 for soft-clipped supplementaries (-Y). The fake mapper
        # follows minimap2: hard-clipped supplementary records unless -Y is given.
        # Record order primary, 700, 400 is valid for a name-sorted BAM (one read name).
        def minimap2(mm2_cmd, out):
            if '-Y' in mm2_cmd:
                supplementary = [
                    segment('fwd', 700, '200S100M', read, supplementary=True),
                    segment('fwd', 400, '100S100M100S', read, supplementary=True),
                ]
            else:
                supplementary = [
                    segment('fwd', 700, '200H100M', read[200:300], supplementary=True),
                    segment('fwd', 400, '100H100M100H', read[100:200], supplementary=True),
                ]
            write_bam(out, [segment('fwd', 100, '100M200S', read)] + supplementary)

        rows, _ = self.run_detector(minimap2)
        self.assertEqual(
            rows, [],
            'forward read 100->400->700 with hard-clipped supplementary records was reported as a back-jump',
        )

    def test_c6_overlapping_forward_split_is_not_a_back_jump(self):
        """Two overlapping alignments on one forward diagonal give no back-jump."""
        # The read is the contiguous reference piece [100, 750). Segment A covers
        # query [0, 400) = ref [100, 500); segment B covers query [250, 650) = ref [350, 750).
        read = REF[100:750]
        records = [
            segment('collinear', 100, '400M250S', read),
            segment('collinear', 350, '250S400M', read, supplementary=True),
        ]
        rows, _ = self.run_on_records(records)
        self.assertEqual(
            rows, [],
            'collinear overlapping alignments of an unrearranged read were reported as a back-jump',
        )

    def test_c7_single_read_does_not_meet_two_read_threshold(self):
        """One read with two equal junctions is not a call at the default two-read threshold."""
        # Rolling-circle read: three consecutive copies of ref [100, 400).
        # Design assumption: min_supporting_reads counts distinct reads (its config
        # comment says "min reads per partial circle call").
        read = REF[100:400] * 3
        records = [
            segment('one_read', 100, '300M600S', read),
            segment('one_read', 100, '300S300M300S', read, supplementary=True),
            segment('one_read', 100, '600S300M', read, supplementary=True),
        ]
        rows, calls = self.run_on_records(records)
        self.assertTrue(rows, 'precondition: the rolling-circle read should still give back-jumps')
        self.assertEqual(calls, [], 'a single read satisfied the two-read support threshold')

    def test_c8_breakpoints_joined_through_a_third_read_form_one_cluster(self):
        """Breakpoints linked through a shared neighbour within tolerance end in one call."""
        # Tolerance 5: (100,400)-(104,405) and (104,405)-(101,410) are neighbours,
        # (100,400)-(101,410) are not. Design assumption: single linkage, as the
        # _cluster_breakpoints docstring states.
        records = (
            back_jump_read('r1', 100, 400)
            + back_jump_read('r2', 101, 410)
            + back_jump_read('r3', 104, 405)
        )
        rows, calls = self.run_on_records(records, breakpoint_tolerance=5)
        self.assertEqual(len(rows), 3, 'precondition: three reads, one back-jump each')
        self.assertEqual(
            [call['n_supporting_reads'] for call in calls], [3],
            'three reads connected by single linkage were not reported as one three-read call',
        )


class GuideAlignmentTest(unittest.TestCase):
    def test_c9_exact_nine_base_hit_survives_a_mismatching_next_base(self):
        """An exact 9 bp match is reported even when the next base mismatches."""
        finder = ShortAlignmentFinder(check_revcomp=False)
        hits = finder.find_alignments_between('ACGTACGTAA', 'ACGTACGTAC')
        self.assertTrue(
            any(hit['mismatches'] == 0 and hit['length'] >= 9 for hit in hits),
            f'exact 9 bp match ACGTACGTA at (0, 0) was not reported; hits: {hits}',
        )

    def test_c10_twelve_base_one_mismatch_alignment_is_found(self):
        """A 12 bp alignment with one central mismatch (allowed by the tiered rule) is reported."""
        finder = ShortAlignmentFinder(check_revcomp=False)
        hits = finder.find_alignments_between('ACGGTCATTGCA', 'ACGGTGATTGCA')
        self.assertTrue(
            any(hit['length'] == 12 and hit['mismatches'] == 1 for hit in hits),
            f'12 bp / 1 mismatch alignment was not reported; hits: {hits}',
        )

    def test_c11_reverse_complement_hit_is_extended_to_full_length(self):
        """A sequence against its own reverse complement gives one full-length exact hit."""
        seq = 'ACGTTGCAAGTCCTAGGCAT'
        finder = ShortAlignmentFinder(check_forward=False)
        hits = finder.find_alignments_between(seq, reverse_complement(seq))
        # Positions are not checked: the coordinate convention for reverse-complement
        # hits is left to the fix.
        self.assertTrue(
            any(hit['length'] == len(seq) and hit['mismatches'] == 0 for hit in hits),
            f'no 20 bp exact reverse-complement hit; lengths found: {[hit["length"] for hit in hits]}',
        )


class CircleCacheTest(unittest.TestCase):
    def test_c15_changed_candidate_table_is_not_scored_with_old_bait_and_bam(self):
        """After the candidate table changes, the new candidate keeps its junction read support."""
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        root = Path(temp.name)
        fastq = root / f'{SAMPLE}.fastq'
        fastq.write_text('@r\nA\n+\nI\n')

        def table(uuid):
            path = root / f'{uuid}.table.txt'
            path.write_text(
                'UUID\tChrom\tStart\tEnd\tFamily\tSubfamily\tConsensus\n'
                f'{uuid}\tchr1\t5000\t5400\tIS110\tNA\t{CONSENSUS}\n'
            )
            return path

        # Fake minimap2|samtools sort: one read across the tail-head junction of
        # every bait in the FASTA it is given.
        def pipe(mm2_cmd, sort_cmd, *args, **kw):
            bait = next(arg for arg in mm2_cmd if str(arg).endswith('.fa'))
            name, seq = Path(bait).read_text().split()[:2]
            junction = len(seq) // 2
            read = segment('junction_read', junction - 150, '300M', seq[junction - 150:junction + 150])
            write_bam(
                sort_cmd[sort_cmd.index('-o') + 1], [read],
                ref_name=name[1:], ref_length=len(seq), sort_order='coordinate',
            )

        def samtools_index(cmd, *args, **kw):
            pysam.index(cmd[-1])

        def tail_head_support(uuid):
            with patch('shutil.which', return_value='fake'):
                finder = CircleFinder(output_dir=str(root / 'out'))
            with patch('cycle.circle_detect.circle_finder.run_checked_pipe', side_effect=pipe), patch(
                'cycle.circle_detect.circle_finder.subprocess.run', side_effect=samtools_index,
            ):
                sample_dir = finder.run_sample(table(uuid), SAMPLE, fastq)
            self.assertIsNotNone(sample_dir, 'run_sample failed')
            with open(sample_dir / f'{SAMPLE}_circle_summary.tsv') as handle:
                return {row['is_uuid']: int(row['n_tail_head_reads']) for row in csv.DictReader(handle, delimiter='\t')}

        self.assertEqual(tail_head_support('old'), {'old': 1}, 'precondition: first run finds the junction read')
        # Same element and reads, new UUID (rebuilt Sniffles tables get new UUIDs).
        # Design assumption: the fix rebuilds stale bait/BAM rather than refusing to run.
        self.assertEqual(
            tail_head_support('new'), {'new': 1},
            'second run reused the bait/BAM of the old candidate table and lost the junction support',
        )


if __name__ == '__main__':
    unittest.main()
