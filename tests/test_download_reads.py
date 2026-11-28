"""Runs labelled long-read that hold short reads are refused; long reads split over _1/_2 are joined."""

import gzip
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from cycle.download_manager.sra_download import SRADownloader
from cycle.utils import find_fastq


def fastq(length, n=5, tag="r"):
    return gzip.compress("".join(f"@{tag}{i}\n{'A' * length}\n+\n{'I' * length}\n" for i in range(n)).encode())


class TestDownloadedReads(unittest.TestCase):
    def download(self, files):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        out = Path(tmp.name)

        def kingfisher(cmd, **kwargs):
            for name, data in files.items():
                (out / name).write_bytes(data)
            return subprocess.CompletedProcess(cmd, 0, "", "")

        with patch("shutil.which", return_value="/fake/kingfisher"):
            downloader = SRADownloader(output_dir=str(out))
        with patch("cycle.download_manager.sra_download.subprocess.run", side_effect=kingfisher):
            table = downloader.download_batch(pd.DataFrame({"srr_accession": ["SRR1"]}))
        return out, table.iloc[0]

    def test_short_read_mates_are_refused(self):
        out, row = self.download({"SRR1_1.fastq.gz": fastq(301), "SRR1_2.fastq.gz": fastq(301)})
        self.assertEqual((row["download_status"], row["fastq_path"]), ("short_reads", ""))

    def test_short_reads_in_one_file_are_refused(self):
        out, row = self.download({"SRR1_subreads.fastq.gz": fastq(151)})
        self.assertEqual(row["download_status"], "short_reads")

    def test_long_read_mates_are_joined_as_single_end_reads(self):
        out, row = self.download({"SRR1_1.fastq.gz": fastq(5000, 3, "a"), "SRR1_2.fastq.gz": fastq(4000, 2, "b")})
        self.assertEqual((row["download_status"], row["fastq_path"]), ("ok", str(out / "SRR1.fastq.gz")))
        self.assertEqual(sorted(p.name for p in out.iterdir()), ["SRR1.fastq.gz"])
        names = gzip.decompress((out / "SRR1.fastq.gz").read_bytes()).decode().splitlines()[::4]
        self.assertEqual(names, ["@a0", "@a1", "@a2", "@b0", "@b1"])
        self.assertEqual(find_fastq(out, "SRR1"), out / "SRR1.fastq.gz")


if __name__ == "__main__":
    unittest.main()
