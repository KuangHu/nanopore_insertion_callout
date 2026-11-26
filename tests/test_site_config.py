"""config/load_site.sh: unfilled template values stop the job; filled ones load."""

import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
BASH = os.environ.get('BASH') or shutil.which('bash')


def load(site_text):
    with tempfile.TemporaryDirectory() as tmp:
        site = Path(tmp) / 'site.env'
        site.write_text(site_text, encoding='utf-8')
        env = {k: v for k, v in os.environ.items() if not k.startswith('CYCLE_')}
        env.update(CYCLE_ROOT=ROOT.as_posix(), CYCLE_SITE_ENV=site.as_posix())
        return subprocess.run(
            [BASH, '-c', 'source "$CYCLE_ROOT/config/load_site.sh" && echo "loaded ${CYCLE_SBATCH_SITE[*]}"'],
            env=env, capture_output=True, text=True)


class TestSitePlaceholders(unittest.TestCase):
    def setUp(self):
        self.template = (ROOT / 'config' / 'lawrencium.env.example').read_text(encoding='utf-8')

    def filled(self):
        text = self.template.replace(':-CHANGE_ME/', ':-/data/')
        return text.replace(':-CHANGE_ME}', ':-12}')

    def test_untouched_template_is_refused(self):
        result = load(self.template)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('CHANGE_ME', result.stderr)
        self.assertNotIn('loaded', result.stdout)

    def test_unfilled_environment_prefix_is_refused(self):
        text = self.filled().replace(':-/data/envs/opfi', ':-CHANGE_ME/envs/opfi')
        result = load(text)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('CYCLE_ENV_PREFIX', result.stderr)

    def test_unfilled_gtdb_paths_are_allowed(self):
        text = self.filled()
        for name in ('GTDB/gtdb_genomes_reps_r220', 'gtdb_blastdb/gtdb_reps_r220'):
            text = text.replace(':-/data/' + name, ':-CHANGE_ME/' + name)
        result = load(text)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('loaded --partition=12 --qos=12 --account=12', result.stdout)

    def test_igi_settings_load(self):
        result = load((ROOT / 'config' / 'igi.env').read_text(encoding='utf-8'))
        self.assertEqual(result.returncode, 0, result.stderr)


if __name__ == '__main__':
    unittest.main()
