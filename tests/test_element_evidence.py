"""Orientation and evidence from ISfinder hits and transposase domains."""
import json
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from cycle.element_evidence import annotate_sample


def record(sequence, orfs, evidence="none"):
    return {"is_element": {"sequence": sequence}, "evidence": evidence,
            "orf_annotation": {"orfs": orfs}}


class ElementEvidence(unittest.TestCase):
    def run_on(self, records, blast_rows, domain_rows):
        def fake(cmd, **kwargs):
            if cmd[0] == "blastn":
                return subprocess.CompletedProcess(cmd, 0, "".join(row + "\n" for row in blast_rows), "")
            table = Path(cmd[cmd.index("--tblout") + 1])
            table.write_text("".join(f"{name} - {domain} PF0 1e-30 {score} 0\n" for name, domain, score in domain_rows))
            return subprocess.CompletedProcess(cmd, 0, "", "")
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "S_is_records_annotated.json"
            path.write_text(json.dumps(records))
            with patch("cycle.element_evidence.subprocess.run", side_effect=fake):
                annotate_sample(path, Path(tmp) / "db", Path(tmp) / "hmm")
            return json.loads(path.read_text())

    def test_minus_strand_hit_and_orf_agree(self):
        orfs = [{"start": 1, "end": 300, "strand": "-", "protein_sequence": "MKV*"}]
        out = self.run_on([record("A" * 1000, orfs)],
                          ["0\tIS1_IS110_IS1111\t98\t990\t5\t995\t1000\t10\t0\t1800\t1000"],
                          [("0_0", "DEDD_Tnp_IS110", 80.0)])[0]
        self.assertEqual((out["orientation"], out["orientation_status"]), ("-", "agree"))
        self.assertEqual(out["evidence"], "is_reference_hit")
        self.assertEqual(out["isfinder_hit"]["family"], "IS110")

    def test_disagreement_follows_isfinder(self):
        orfs = [{"start": 1, "end": 300, "strand": "+", "protein_sequence": "MKV*"}]
        out = self.run_on([record("A" * 1000, orfs)],
                          ["0\tIS1_IS110_IS1111\t98\t990\t5\t995\t1000\t10\t0\t1800\t1000"],
                          [("0_0", "DEDD_Tnp_IS110", 80.0)])[0]
        self.assertEqual((out["orientation"], out["orientation_status"]), ("-", "disagree"))

    def test_domain_only_and_short_hit(self):
        orfs = [{"start": 1, "end": 300, "strand": "-", "protein_sequence": "MKV*"}]
        # A hit over 10% of the record is not a reference hit
        out = self.run_on([record("A" * 1000, orfs, evidence="tail_head_reads")],
                          ["0\tIS1_IS110_IS1111\t98\t100\t5\t104\t1\t100\t0\t180\t1000"],
                          [("0_0", "DEDD_Tnp_IS110", 80.0)])[0]
        self.assertIsNone(out["isfinder_hit"])
        self.assertEqual((out["orientation"], out["orientation_status"]), ("-", "orf_only"))
        self.assertEqual(out["evidence"], "tail_head_reads")
        out = self.run_on([record("A" * 1000, orfs)], [], [("0_0", "DEDD_Tnp_IS110", 80.0)])[0]
        self.assertEqual(out["evidence"], "transposase_domain")


if __name__ == "__main__":
    unittest.main()
