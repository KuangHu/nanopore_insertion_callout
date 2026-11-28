"""Insertions next to a resident copy: the reads align the inserted part to the copy."""

import random
import tempfile
import unittest
from pathlib import Path

import pysam

from cycle.preprocess.sniffles_runner import SnifflesRunner

RNG = random.Random(11)
LEFT = "".join(RNG.choice("ACGT") for _ in range(1000))
ELEMENT = "".join(RNG.choice("ACGT") for _ in range(1500))
RIGHT = "".join(RNG.choice("ACGT") for _ in range(1000))
READ = LEFT + ELEMENT + RIGHT


def write_bam(path, copy_at, reads):
    """Each read: left flank 9000-10000 (primary), the element at copy_at, right flank from 10000."""
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 100000}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        for name, reverse in reads:
            strand = "-" if reverse else "+"
            read = pysam.AlignedSegment()
            read.query_name = name
            read.query_sequence = READ
            read.flag = 16 if reverse else 0
            read.reference_id = 0
            read.reference_start = 9000
            read.mapping_quality = 60
            read.cigartuples = [(0, 1000), (4, 2500)]
            read.set_tag("SA", f"chr1,{copy_at + 1},{strand},1000S1500M1000S,60,0;chr1,10001,{strand},2500S1000M,60,0;")
            bam.write(read)
    pysam.index(str(path))


class TestSplitReadInsertions(unittest.TestCase):
    def call(self, copy_at, reads, called=()):
        with tempfile.TemporaryDirectory() as tmp:
            bam = Path(tmp) / "S.sorted.bam"
            write_bam(bam, copy_at, reads)
            return SnifflesRunner._split_read_insertions(bam, list(called), 500, 20000, 3, "S")

    def test_element_aligned_to_a_copy_elsewhere_is_an_insertion(self):
        found = self.call(60000, [("r1", False), ("r2", False), ("r3", True)])
        self.assertEqual(len(found), 1)
        self.assertEqual((found[0]["Chrom"], found[0]["Start"], found[0]["Support"]), ("chr1", 10000, 3))
        self.assertEqual(found[0]["Sequence"], ELEMENT)
        self.assertEqual(found[0]["Rnames"], ["r1", "r2", "r3"])

    def test_too_few_reads_or_a_copy_beside_the_site_or_an_existing_call(self):
        self.assertEqual(self.call(60000, [("r1", False), ("r2", False)]), [])
        self.assertEqual(self.call(11000, [("r1", False), ("r2", False), ("r3", False)]), [])
        called = [{"Chrom": "chr1", "Start": 10050}]
        self.assertEqual(self.call(60000, [("r1", False), ("r2", False), ("r3", False)], called), [])


if __name__ == "__main__":
    unittest.main()
