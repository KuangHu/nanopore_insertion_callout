"""Every output directory says which code, environment and parameters made it."""

import json
import tempfile
import unittest
from pathlib import Path

from cycle.preprocess.sniffles_runner import SnifflesRunner
from cycle.provenance import write_provenance


class TestProvenance(unittest.TestCase):
    def test_written_to_each_directory(self):
        with tempfile.TemporaryDirectory() as tmp:
            dirs = [Path(tmp) / "a", Path(tmp) / "b"]
            write_provenance(dirs, "demo", {"threads": 4, "out": Path("/x")}, {"note": "n"})
            infos = [json.loads((d / "provenance_demo.json").read_text()) for d in dirs]
        self.assertEqual(infos[0], infos[1])
        self.assertRegex(infos[0]["code_commit"], r"^([0-9a-f]{40}|unknown)$")
        self.assertEqual(infos[0]["parameters"], {"threads": 4, "out": "/x"})
        self.assertEqual(infos[0]["note"], "n")
        self.assertIn("minimap2", infos[0]["tools"])


class TestSnifflesOptions(unittest.TestCase):
    def test_paths_do_not_count_but_options_do(self):
        a = SnifflesRunner.command("/a.bam", "/a.vcf", "/r1.fa", 500, 3, False)
        b = SnifflesRunner.command("/b.bam", "/b.vcf", "/r2.fa", 500, 3, False)
        self.assertEqual(SnifflesRunner._options(a), SnifflesRunner._options(b))
        old = [token for token in a if token != "--all-contigs"]
        self.assertNotEqual(SnifflesRunner._options(old), SnifflesRunner._options(a))
        self.assertIn("--all-contigs", SnifflesRunner._options(a))


if __name__ == "__main__":
    unittest.main()
