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
    def run_on(self, records, blast_rows, domain_rows, entry_strands=None):
        def fake(cmd, **kwargs):
            if cmd[0] == "blastn":
                return subprocess.CompletedProcess(cmd, 0, "".join(row + "\n" for row in blast_rows), "")
            table = Path(cmd[cmd.index("--tblout") + 1])
            table.write_text("".join(f"{name} - {domain} PF0 1e-30 {score} 0\n" for name, domain, score in domain_rows))
            return subprocess.CompletedProcess(cmd, 0, "", "")
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "S_is_records_annotated.json"
            path.write_text(json.dumps(records))
            if entry_strands:
                (Path(tmp) / "db.tnp_strand.tsv").write_text("".join(f"{e}\t{st}\n" for e, st in entry_strands.items()))
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

    def test_entry_stored_with_its_transposase_on_the_minus_strand(self):
        """A plus-strand hit to an entry whose own transposase is on its minus strand predicts a minus ORF."""
        orfs = [{"start": 1, "end": 300, "strand": "-", "protein_sequence": "MKV*"}]
        out = self.run_on([record("A" * 1000, orfs)],
                          ["0\tIS609_IS200/IS605_IS605\t98\t990\t5\t995\t10\t1000\t0\t1800\t1000"],
                          [("0_0", "Y1_Tnp", 130.0)], {"IS609_IS200/IS605_IS605": "-"})[0]
        self.assertEqual((out["orientation"], out["orientation_status"]), ("+", "agree"))

    def test_transposase_beside_the_hit(self):
        """The best transposase ORF lies outside the matched element: two elements, no single orientation check."""
        orfs = [{"start": 1600, "end": 2700, "strand": "-", "protein_sequence": "MKV*"}]
        out = self.run_on([record("A" * 2800, orfs)],
                          ["0\tISLpn2_IS4_IS50\t99\t1500\t9\t1506\t1\t1498\t0\t2700\t1498"],
                          [("0_0", "DDE_Tnp_ISL3", 216.0)], {"ISLpn2_IS4_IS50": "+"})[0]
        self.assertEqual((out["orientation"], out["orientation_status"]), ("+", "two_elements"))


class EndRepeatPlacement(ElementEvidence):
    """Bases shared by the element's end and a flank go to the element's 5' end, as in ISfinder."""

    def placed(self, element, up, down, hit_strand, orfs=None):
        rec = record(element, orfs or [{"start": 5, "end": 40, "strand": hit_strand, "protein_sequence": "MKV*"}])
        rec["is_element"].update(start=1000, end=1000)
        rec["insertion_site"] = {"pos_5p": 1000, "pos_3p": 1000}
        rec["flanking_upstream"] = {"sequence": up, "length": len(up)}
        rec["flanking_downstream"] = {"sequence": down, "length": len(down)}
        n = len(element)
        sstart, send = (1, n) if hit_strand == "+" else (n, 1)
        return self.run_on([rec], [f"0\tIS1_IS110_IS1111\t99\t{n}\t1\t{n}\t{sstart}\t{send}\t0\t90\t{n}"], [])[0]

    def test_plus_element_ending_like_its_upstream_flank_moves_back(self):
        # filled locus GGGTTACT | ACGTACGTCT | AAA; the CT can sit at the element's 5' end
        out = self.placed("ACGTACGTCT" * 5, "GGGTTACT", "AAACCC", "+")
        self.assertEqual(out["end_repeat_shift"], -2)
        self.assertEqual(out["is_element"]["sequence"], "CT" + ("ACGTACGTCT" * 5)[:-2])
        self.assertEqual(out["flanking_upstream"]["sequence"], "GGGTTA")
        self.assertEqual(out["flanking_downstream"]["sequence"], "CTAAACCC")
        self.assertEqual((out["is_element"]["start"], out["insertion_site"]["pos_5p"]), (998, 998))
        self.assertEqual((out["orf_annotation"]["orfs"][0]["start"], out["orf_annotation"]["orfs"][0]["end"]), (7, 42))
        self.assertEqual((out["isfinder_hit"]["query_start"], out["isfinder_hit"]["query_end"]), (1, 50))

    def test_minus_element_starting_like_its_downstream_flank_moves_on(self):
        out = self.placed("GAACGTACGT" * 5, "GGGTTT", "GACCC", "-")
        self.assertEqual(out["end_repeat_shift"], 2)
        self.assertEqual(out["is_element"]["sequence"], ("GAACGTACGT" * 5)[2:] + "GA")
        self.assertEqual(out["flanking_upstream"]["sequence"], "GGGTTTGA")
        self.assertEqual(out["flanking_downstream"]["sequence"], "CCC")
        self.assertEqual(out["is_element"]["start"], 1002)
        self.assertEqual((out["orf_annotation"]["orfs"][0]["start"], out["orf_annotation"]["orfs"][0]["end"]), (3, 38))

    def test_repeat_already_at_the_5_prime_end_stays(self):
        # plus element starting like its downstream flank: the shared bases are at its 5' end already
        out = self.placed("CTACGTACGA" * 5, "GGGTTT", "CTAAAC", "+")
        self.assertEqual(out["end_repeat_shift"], 0)
        self.assertEqual(out["is_element"]["sequence"], "CTACGTACGA" * 5)
        self.assertEqual(out["is_element"]["start"], 1000)


if __name__ == "__main__":
    unittest.main()
