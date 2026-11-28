"""Sniffles calls insertions on contigs shorter than 1 Mb (plasmids), not only on the chromosome."""

import random
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

from cycle.preprocess.sniffles_runner import SnifflesRunner


@unittest.skipUnless(all(shutil.which(t) for t in ("sniffles", "minimap2", "samtools")), "sniffles, minimap2, samtools needed")
class TestSmallContig(unittest.TestCase):
    def test_insertion_on_a_30_kb_contig_is_called(self):
        rng = random.Random(7)
        plasmid = "".join(rng.choice("ACGT") for _ in range(30000))
        element = "".join(rng.choice("ACGT") for _ in range(1500))
        filled = plasmid[:15000] + element + plasmid[15000:]
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            ref = root / "ref.fna"
            ref.write_text(f">plasmid\n{plasmid}\n")
            reads = root / "reads.fastq"
            with open(reads, "w") as handle:
                for i in range(20):
                    start = rng.randrange(9000, 12000)
                    seq = filled[start:start + 8000]
                    handle.write(f"@r{i}\n{seq}\n+\n{'I' * len(seq)}\n")
            bam = root / "S.sorted.bam"
            subprocess.run(f"minimap2 -a --MD -x map-ont {ref} {reads} | samtools sort -o {bam} -",
                           shell=True, check=True, capture_output=True)
            subprocess.run(["samtools", "index", str(bam)], check=True)
            runner = SnifflesRunner(output_dir=str(root / "sniffles"), alignment_dir=str(root))
            table = runner.run_sample(bam, ref, "S")
            self.assertIsNotNone(table)
            rows = Path(table).read_text().splitlines()[1:]
            starts = [int(row.split("\t")[2]) for row in rows if row.split("\t")[1] == "plasmid"]
            self.assertTrue(any(abs(s - 15000) <= 50 for s in starts), f"no call near 15000 on the plasmid: {rows}")


if __name__ == "__main__":
    unittest.main()
