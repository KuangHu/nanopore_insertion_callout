"""Per-sample minimap2 preset and one index per reference and preset."""

import gzip
import random
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

from cycle.preprocess.aligner import Aligner, pick_preset


def write_fastq(path, quality_char, n=50, length=200, gz=False):
    random.seed(1)
    text = "".join(f"@r{i}\n{''.join(random.choice('ACGT') for _ in range(length))}\n+\n{quality_char * length}\n"
                   for i in range(n))
    if gz:
        path.write_bytes(gzip.compress(text.encode()))
    else:
        path.write_text(text)
    return path


class TestPickPreset(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def test_nanopore_gets_map_ont_whatever_the_quality(self):
        fastq = write_fastq(self.root / "ont.fastq", "I")
        self.assertEqual(pick_preset(fastq, "OXFORD_NANOPORE"), "map-ont")

    def test_pacbio_q30_reads_are_hifi(self):
        fastq = self.root / "hifi.fastq.gz"
        random.seed(2)
        text = "".join(f"@r{i}\nACGTACGTAC\n+\n{''.join(random.choice('?@ABCDEFGHI') for _ in range(10))}\n"
                       for i in range(50))  # Q30 to Q40, not all one value
        fastq.write_bytes(gzip.compress(text.encode()))
        self.assertEqual(pick_preset(fastq, "PACBIO_SMRT"), "map-hifi")

    def test_pacbio_reads_without_qualities_are_clr(self):
        fastq = write_fastq(self.root / "flat.fastq", "!")
        self.assertEqual(pick_preset(fastq, "PACBIO_SMRT"), "map-pb")
        fastq = write_fastq(self.root / "flat_high.fastq", "I")  # one value repeated carries no quality
        self.assertEqual(pick_preset(fastq, "PACBIO_SMRT"), "map-pb")

    def test_pacbio_q10_reads_are_clr(self):
        fastq = self.root / "clr.fastq"
        fastq.write_text("".join(f"@r{i}\nACGTACGTAC\n+\n++,,++,,++\n" for i in range(50)))  # Q10 to Q11
        self.assertEqual(pick_preset(fastq, "PACBIO_SMRT"), "map-pb")


@unittest.skipUnless(shutil.which("minimap2") and shutil.which("samtools"), "minimap2 and samtools needed")
class TestIndexPerPreset(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        random.seed(3)
        genome = "".join(random.choice("ACGT") for _ in range(20000))
        self.fasta = self.root / "ref" / "GCF_1_genomic.fna"
        self.fasta.parent.mkdir()
        self.fasta.write_text(f">chr\n{genome}\n")
        self.reads = self.root / "reads.fastq"
        self.reads.write_text("".join(f"@r{i}\n{genome[i * 1000:i * 1000 + 3000]}\n+\n{'I' * 3000}\n" for i in range(15)))

    def tearDown(self):
        self.tmp.cleanup()

    def test_each_preset_has_its_own_index_and_alignment_uses_it(self):
        aligner = Aligner(output_dir=str(self.root / "aln"), threads=1)
        ont = aligner.index(self.fasta, "map-ont")
        pb = aligner.index(self.fasta, "map-pb")
        self.assertEqual(ont.name, "GCF_1_genomic.map-ont.mmi")
        self.assertEqual(pb.name, "GCF_1_genomic.map-pb.mmi")
        self.assertNotEqual(ont.read_bytes(), pb.read_bytes())
        self.assertEqual(list(self.fasta.parent.glob(".*pending")), [])

        results = aligner.align_batch([{"sample_id": "S1", "fastq": self.reads, "reference_fasta": self.fasta,
                                        "platform": "PACBIO_SMRT"}])
        bam = results[0]["bam"]
        header = subprocess.run(["samtools", "view", "-H", str(bam)], capture_output=True, text=True, check=True).stdout
        command = next(line for line in header.splitlines() if line.startswith("@PG") and "ID:minimap2" in line)
        self.assertIn("-x map-pb", command)  # one quality value repeated: CLR
        self.assertIn(str(pb), command.split())


if __name__ == "__main__":
    unittest.main()
