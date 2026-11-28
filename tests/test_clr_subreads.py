"""CLR runs: one read per ZMW before alignment, one call per site after Sniffles."""

import gzip
import tempfile
import unittest
from pathlib import Path

from cycle.preprocess.aligner import longest_per_zmw, zmw_key
from cycle.preprocess.sniffles_runner import SnifflesRunner


def fastq(records):
    return "".join(f"@{name}\n{seq}\n+\n{'!' * len(seq)}\n" for name, seq in records)


class TestLongestPerZmw(unittest.TestCase):
    def test_zmw_from_sra_spot_and_pacbio_name(self):
        self.assertEqual(zmw_key("@SRR7818938.198 /3\n"), "SRR7818938.198")
        self.assertEqual(zmw_key("@SRR1.5 m54006_170101_000000/4194/0_1200\n"), "m54006_170101_000000/4194")
        self.assertEqual(zmw_key("@m54006_170101_000000/4194/1300_2500\n"), "m54006_170101_000000/4194")

    def test_keeps_the_longest_subread_of_each_zmw(self):
        records = [("SRR1.1 /1", "A" * 50), ("SRR1.2 /1", "C" * 30), ("SRR1.2 /3", "G" * 80),
                   ("SRR1.2 /5", "T" * 60), ("SRR1.3 /1", "A" * 10)]
        with tempfile.TemporaryDirectory() as tmp:
            source, out = Path(tmp) / "r.fastq.gz", Path(tmp) / "out.fastq"
            with gzip.open(source, "wt") as handle:
                handle.write(fastq(records))
            self.assertEqual(longest_per_zmw(source, out), (5, 3))
            lines = out.read_text().splitlines()
        self.assertEqual([lines[i] for i in range(0, 12, 4)], ["@SRR1.1 /1", "@SRR1.2 /3", "@SRR1.3 /1"])
        self.assertEqual(lines[5], "G" * 80)


def call(chrom, start, support):
    return {"Chrom": chrom, "Start": start, "Support": support, "Sequence": "ACGT"}


class TestMergeNearby(unittest.TestCase):
    def test_best_supported_call_of_a_site_is_kept(self):
        calls = [call("c1", 1000, 3), call("c1", 1004, 9), call("c1", 1012, 2), call("c1", 1500, 4), call("c2", 1003, 1)]
        kept = SnifflesRunner._merge_nearby(calls)
        self.assertEqual([(c["Chrom"], c["Start"]) for c in kept], [("c1", 1004), ("c1", 1500), ("c2", 1003)])

    def test_calls_further_apart_stay(self):
        calls = [call("c1", 1000, 3), call("c1", 1011, 9)]
        self.assertEqual(len(SnifflesRunner._merge_nearby(calls)), 2)


if __name__ == "__main__":
    unittest.main()
