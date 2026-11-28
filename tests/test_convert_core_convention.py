"""Cross_reference_IS (core at the 3' end) to Cycle (core at the 5' end) and back."""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from convert_core_convention import core_shift  # noqa: E402
from cycle.element_evidence import _place_end_repeat  # noqa: E402

LEFT = "GATCCGTTAGCAAT" + "CT"        # the flank ends with the core
BODY = "GGCATTTACGGATCCAAGTTTGCA"
RIGHT = "AGGTCCATTGAC"
# Plus element, core CT at its 3' end (Cross_reference_IS): ... CT | BODY CT | RIGHT
PLUS = LEFT + BODY + "CT" + RIGHT
PLUS_START, PLUS_END = len(LEFT), len(LEFT) + len(BODY) + 2


class TestCoreShift(unittest.TestCase):
    def test_plus_element_moves_left_to_5prime(self):
        shift = core_shift(PLUS, PLUS_START, PLUS_END, "+", "5prime")
        self.assertEqual(shift, -2)
        self.assertEqual(PLUS[PLUS_START + shift:PLUS_END + shift], "CT" + BODY)

    def test_minus_element_moves_right_to_5prime(self):
        # Minus element, core at its 3' end = its reference-left end: LEFT | AG BODY | AG ...
        genome = "TTGACCAGTA" + "AG" + BODY + "AG" + "CCTTAGGA"
        start, end = 10, 10 + 2 + len(BODY)
        self.assertEqual(core_shift(genome, start, end, "-", "5prime"), 2)

    def test_round_trip(self):
        shift = core_shift(PLUS, PLUS_START, PLUS_END, "+", "5prime")
        back = core_shift(PLUS, PLUS_START + shift, PLUS_END + shift, "+", "3prime")
        self.assertEqual(back, -shift)

    def test_no_shared_bases_no_shift(self):
        genome = "AAAA" + "CCGGTTGCG" + "TTTT"
        self.assertEqual(core_shift(genome, 4, 13, "+", "5prime"), 0)

    def test_agrees_with_the_pipeline_placement(self):
        record = {"orientation": "+", "is_element": {"sequence": PLUS[PLUS_START:PLUS_END], "start": 0, "end": 0},
                  "flanking_upstream": {"sequence": LEFT}, "flanking_downstream": {"sequence": RIGHT}}
        self.assertEqual(_place_end_repeat(record), core_shift(PLUS, PLUS_START, PLUS_END, "+", "5prime"))


if __name__ == "__main__":
    unittest.main()
