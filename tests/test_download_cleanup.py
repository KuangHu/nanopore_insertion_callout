"""The download step leaves no aria2 control files beside the FASTQ."""

import gzip
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from cycle.download_manager.sra_download import SRADownloader


class TestAria2Cleanup(unittest.TestCase):
    def test_control_file_is_removed_and_fastq_kept(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)

            def kingfisher(cmd, **kwargs):
                (out / "SRR1.fastq.gz").write_bytes(gzip.compress(f"@r\n{'A' * 1500}\n+\n{'I' * 1500}\n".encode()))
                (out / "SRR1.fastq.gz.aria2").write_bytes(b"control")
                return subprocess.CompletedProcess(cmd, 0, "", "")

            with patch("shutil.which", return_value="/fake/kingfisher"):
                downloader = SRADownloader(output_dir=str(out))
            with patch("cycle.download_manager.sra_download.subprocess.run", side_effect=kingfisher):
                result = downloader.download_one("SRR1")
            self.assertEqual(result, out / "SRR1.fastq.gz")
            self.assertEqual(sorted(p.name for p in out.iterdir()), ["SRR1.fastq.gz"])


if __name__ == "__main__":
    unittest.main()
