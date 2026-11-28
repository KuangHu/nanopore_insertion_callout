"""Records beside their own copy in the reference: junction reads of a tandem look like circle reads."""

import random
import tempfile
import unittest
from pathlib import Path

from cycle.is_formatter import ISFormatter

RNG = random.Random(5)


def seq(n):
    return "".join(RNG.choice("ACGT") for _ in range(n))


UNIT = seq(1000)
REFERENCE = seq(10000) + UNIT + seq(19000)  # one copy of UNIT at 10000-11000
OTHER = seq(1000)                           # in no reference


def record(name, sequence, start, group=None):
    return {"is_id": name, "element_group": group or name, "circle_evidence": {},
            "is_element": {"sequence": sequence, "contig": "chr1", "start": start}}


class TestTandemFlag(unittest.TestCase):
    def flag(self, records, with_reference=True):
        with tempfile.TemporaryDirectory() as tmp:
            reference = Path(tmp) / "ref.fa"
            reference.write_text(f">chr1\n{REFERENCE}\n")
            formatter = ISFormatter(output_dir=tmp, threads=1)
            formatter._flag_tandem(records, reference if with_reference else None, Path(tmp))
        return {r["is_id"]: r["circle_evidence"] for r in records}

    def test_copy_beside_the_site_is_flagged(self):
        flags = self.flag([record("dup", UNIT, 11000), record("new", OTHER, 11000)])
        self.assertIs(flags["dup"]["tandem_confounded"], True)
        self.assertGreaterEqual(flags["dup"]["local_copy_fraction"], 0.9)
        self.assertIs(flags["new"]["tandem_confounded"], False)
        self.assertEqual(flags["new"]["local_copy_fraction"], 0.0)

    def test_copy_far_from_the_site_is_not_flagged(self):
        flags = self.flag([record("far", UNIT, 29500 + 20000)])
        self.assertIs(flags["far"]["tandem_confounded"], False)

    def test_copy_a_few_kb_away_is_another_insertion(self):
        flags = self.flag([record("near", UNIT, 16000)])
        self.assertIs(flags["near"]["tandem_confounded"], False)

    def test_copies_sharing_counts_are_all_flagged(self):
        flags = self.flag([record("dup", UNIT, 11000, "g"), record("copy", UNIT[:900] + OTHER[:100], 60000, "g")])
        self.assertIs(flags["copy"]["tandem_confounded"], True)

    def test_no_reference_leaves_the_flag_unset(self):
        flags = self.flag([record("dup", UNIT, 11000)], with_reference=False)
        self.assertIsNone(flags["dup"]["tandem_confounded"])


if __name__ == "__main__":
    unittest.main()
