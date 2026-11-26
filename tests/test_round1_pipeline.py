"""Round 1 review findings for the batch pipeline, IS110 collection and Hub upload.

Each test states the correct behaviour and fails on the code as reviewed. External tools
(sbatch, conda, prodigal, hmmsearch, minimap2, samtools, kingfisher, the Hub upload) are replaced;
all inputs are synthetic and live in a temporary directory.
"""

import gzip
import importlib.util
import json
import logging
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import MagicMock, patch

import pandas as pd

from cycle.preprocess.aligner import Aligner
from cycle.utils import find_fastq

ROOT = Path(__file__).resolve().parents[1]
BASH = os.environ.get('BASH') or shutil.which('bash')
REAL_RUN = subprocess.run
REAL_WHICH = shutil.which
FAKED_TOOLS = ('minimap2', 'samtools', 'kingfisher')
CONSENSUS = 'ACGT' * 150  # 600 nt, above the 500 nt cut-off of the IS110 collectors


def load_script(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / 'scripts' / f'{name}.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def fake_which(name, *args, **kwargs):
    return f'/fake/{name}' if name in FAKED_TOOLS else REAL_WHICH(name, *args, **kwargs)


def fake_run(cmd, **kwargs):
    """Stand in for samtools index and minimap2 indexing; every other command really runs."""
    if list(cmd[:2]) == ['samtools', 'index']:
        Path(cmd[-1] + '.bai').write_bytes(b'index')
    elif cmd[0] not in FAKED_TOOLS:
        return REAL_RUN(cmd, **kwargs)
    return subprocess.CompletedProcess(cmd, 0, '', '')


def aligner_factory(aligner):
    """Aligner replacement that, like the real class, creates its output directory."""
    def make(output_dir, **kwargs):
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        return aligner
    return make


def write_tsv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, sep='\t', index=False)


def write_executable(path, text):
    path.write_text(text, encoding='utf-8')
    path.chmod(0o755)


def run_main(pipeline, argv):
    """Run run_pipeline.main() and return its exit code (0 when it returns normally)."""
    with patch.object(sys, 'argv', ['run_pipeline.py'] + [str(a) for a in argv]):
        try:
            pipeline.main()
        except SystemExit as exc:
            return exc.code or 0
    return 0


class TemporaryTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        logging.disable(logging.CRITICAL)
        self.addCleanup(logging.disable, logging.NOTSET)
        if not BASH:
            self.skipTest('bash is required')
        self.stubs = self.root / 'stubs'
        self.stubs.mkdir()

    def shell_env(self, **extra):
        env = dict(os.environ)
        env.update({
            'PATH': os.pathsep.join([str(self.stubs), str(Path(sys.executable).parent), env.get('PATH', '')]),
            'CYCLE_ROOT': ROOT.as_posix(),
            'CYCLE_SITE_ENV': (ROOT / 'config' / 'igi.env').as_posix(),
            'CYCLE_DATA_ROOT': (self.root / 'outputs').as_posix(),
            'CYCLE_LOG_DIR': (self.root / 'logs').as_posix(),
            'CYCLE_ENV_PREFIX': (self.root / 'env').as_posix(),
            'CYCLE_HMM_DIR': (self.root / 'hmm').as_posix(),
            'CYCLE_QUIET': '1',
        })
        env.update(extra)
        return env

    # ── submit_one_batch.sh with sbatch replaced on PATH ─────────────
    def submit(self, batch_arg, nodes=1):
        """Run the real submitter; return (process, jobs) where jobs are the captured sbatch calls."""
        jobs_dir = self.root / 'sbatch_calls'
        jobs_dir.mkdir()
        write_executable(self.stubs / 'sbatch', (
            '#!/bin/bash\n'
            f'dir={shlex.quote(jobs_dir.as_posix())}\n'
            'n=$(( $(cat "$dir/count" 2>/dev/null || echo 100) + 1 ))\n'
            'echo "$n" > "$dir/count"\n'
            'printf "%s\\0" "$@" > "$dir/job_$n.args"\n'
            'echo "$n"\n'))
        write_executable(self.stubs / 'python3', f'#!/bin/bash\nexec {shlex.quote(sys.executable)} "$@"\n')
        env = self.shell_env(CYCLE_BATCH_METADATA_DIR=(self.root / 'batch_metadata').as_posix())
        proc = REAL_RUN([BASH, str(ROOT / 'submit_one_batch.sh'), batch_arg, '--nodes', str(nodes)],
                        capture_output=True, text=True, env=env)
        jobs = []
        for path in sorted(jobs_dir.glob('job_*.args'), key=lambda p: int(p.stem.split('_')[1])):
            args = path.read_bytes().decode().split('\0')
            option = lambda key: next((a[len(key):] for a in args if a.startswith(key)), '')
            jobs.append({'id': path.stem.split('_')[1], 'name': option('--job-name='),
                         'dependency': option('--dependency='), 'wrap': option('--wrap=')})
        return proc, jobs

    # ── scripts/is110_run_all.sh with conda, prodigal and hmmsearch replaced on PATH ──
    def run_is110_all(self):
        write_executable(self.stubs / 'conda', (
            '#!/bin/bash\n'
            '[ "$1" = run ] && shift\n'
            'case "$1" in -p|-n) shift 2 ;; esac\n'
            'exec "$@"\n'))
        # Prodigal stand-in: one protein per input record; a record named failme* makes the worker exit 7.
        write_executable(self.stubs / 'prodigal', (
            f'#!{sys.executable}\n'
            'import sys\n'
            'args = sys.argv[1:]\n'
            'src, dst = args[args.index("-i") + 1], args[args.index("-a") + 1]\n'
            'names = [line[1:].strip() for line in open(src) if line.startswith(">")]\n'
            'open(src + ".called", "w").write("\\n".join(names))\n'
            'if any(name.startswith("failme") for name in names):\n'
            '    sys.exit(7)\n'
            'with open(dst, "w") as out:\n'
            '    for name in names:\n'
            '        out.write(f">{name}_1 # 1 # 300 # 1 # ID=1_1;partial=00\\nMKV*\\n")\n'))
        write_executable(self.stubs / 'hmmsearch', (
            f'#!{sys.executable}\n'
            'import sys\n'
            'open(sys.argv[sys.argv.index("--tblout") + 1], "w").write("# no hits\\n")\n'))
        return REAL_RUN([BASH, str(ROOT / 'scripts' / 'is110_run_all.sh')],
                        capture_output=True, text=True, env=self.shell_env())

    def sniffles_table(self, directory, sample, uuids):
        write_tsv(directory / sample / f'{sample}.table.txt', [
            {'UUID': uuid, 'Chrom': 'chr', 'Start': 100, 'End': 101, 'Consensus': CONSENSUS} for uuid in uuids])

    @property
    def is110_out(self):
        return self.root / 'outputs' / 'is110_all_006_052'

    # ── scripts/is110_build_json.py on synthetic HMM and Prodigal results for one element "u" ──
    def build_is110_json(self):
        out = self.is110_out
        out.mkdir(parents=True, exist_ok=True)
        (out / 'is110_ids.json').write_text('["u"]')
        for name in ('DEDD', 'Tnp20'):
            (out / f'{name}_hits.tbl').write_text('# no hits\n')
        (out / 'all_proteins.faa').write_text('>u_1 # 1 # 300 # 1 # ID=1_1;partial=00\nMKV*\n')
        proc = REAL_RUN([sys.executable, str(ROOT / 'scripts' / 'is110_build_json.py')],
                        capture_output=True, text=True, env=self.shell_env())
        self.assertEqual(proc.returncode, 0, proc.stderr)
        return {rec['is_id']: rec for rec in json.loads((out / 'is110_all_records.json').read_text())}

    # ── scripts/upload_batch_to_hf.py with the Hub replaced ──────────
    def staged_upload(self, batch_dir):
        """Run upload_batch and return {path in repo: bytes} of what would have been sent.

        Folder uploads are listed by the installed huggingface_hub itself, so symlink handling is
        the library's, not this test's. Assumes the fix keeps upload_batch(repo_id, batch_dir,
        batch_name) and sends files through HfApi (upload_folder, upload_large_folder, upload_file
        or create_commit).
        """
        import huggingface_hub
        real_api = huggingface_hub.HfApi
        sent = {}

        def content(source):
            if isinstance(source, bytes):
                return source
            if isinstance(source, (str, os.PathLike)):
                return Path(source).read_bytes()
            return source.read()

        class RecordingApi:
            def __init__(self, *args, **kwargs):
                pass

            def create_repo(self, *args, **kwargs):
                pass

            def upload_folder(self, *, folder_path, path_in_repo=None, allow_patterns=None,
                              ignore_patterns=None, **kwargs):
                for op in real_api._prepare_upload_folder_additions(
                        real_api(), folder_path, path_in_repo or '',
                        allow_patterns=allow_patterns, ignore_patterns=ignore_patterns):
                    sent[op.path_in_repo] = content(op.path_or_fileobj)

            upload_large_folder = upload_folder

            def upload_file(self, *, path_or_fileobj, path_in_repo, **kwargs):
                sent[path_in_repo] = content(path_or_fileobj)

            def create_commit(self, repo_id=None, operations=(), **kwargs):
                for op in operations:
                    if hasattr(op, 'path_or_fileobj'):
                        sent[op.path_in_repo] = content(op.path_or_fileobj)

        uploader = load_script('upload_batch_to_hf')
        (self.root / 'tmp').mkdir()
        with patch.dict(os.environ, {'TMPDIR': str(self.root / 'tmp')}), \
                patch.object(huggingface_hub, 'HfApi', RecordingApi), \
                patch.object(uploader, 'HfApi', RecordingApi):
            uploader.upload_batch('user/dataset', batch_dir, batch_dir.name)
        return sent


class TestRound1Pipeline(TemporaryTest):

    def test_C17_failed_circle_stage_exits_nonzero_and_ignores_stale_summary(self):
        """A sample whose circle stage failed in this run makes the run fail and is not formatted from an old summary."""
        pipeline = load_script('run_pipeline')
        out = self.root / 'out'
        metadata = self.root / 'samples.tsv'
        write_tsv(metadata, [{'srr_accession': 'S', 'organism': 'A'}])
        self.sniffles_table(out / 'sniffles_output', 'S', ['u'])
        stale = out / 'circle_output' / 'S'
        stale.mkdir(parents=True)
        (stale / 'S_circle_summary.tsv').write_text('is_uuid\tn_tail_head_reads\nold\t3\n')

        finder = MagicMock()
        finder.run_batch.return_value = {'S': None}  # every current circle attempt fails
        formatter = MagicMock()
        formatter.run_batch.return_value = {'S': None}
        with patch.object(pipeline, 'CircleFinder', return_value=finder), \
                patch.object(pipeline, 'ISFormatter', return_value=formatter):
            code = run_main(pipeline, ['--metadata', metadata, '--outdir', out, '--steps', 'circle', 'format'])

        self.assertEqual(finder.run_batch.call_count, 1)  # the circle stage was really attempted
        with self.subTest('exit status'):
            self.assertNotEqual(code, 0, 'run_pipeline exited 0 although circle detection failed for its only sample')
        with self.subTest('stale input'):
            handed = [call.kwargs.get('circle_results', {}).get('S') for call in formatter.run_batch.call_args_list]
            self.assertFalse(any(handed), f'formatter was given the old circle directory of a failed sample: {handed}')

    def test_C18_truncated_guide_json_does_not_mark_sample_complete(self):
        """A sample whose guide JSON cannot be parsed is still scheduled for recomputation."""
        write_tsv(self.root / 'batch_metadata' / 'batch_001.tsv', [{'srr_accession': 'sample_A', 'organism': 'A'}])
        guide = self.root / 'outputs' / 'batch_001' / 'is_formatter_output' / 'sample_A' / 'sample_A_is_records_guide.json'
        guide.parent.mkdir(parents=True)
        guide.write_text('{broken')
        proc, _ = self.submit('1')
        self.assertEqual(proc.returncode, 0, proc.stderr)
        # Assumes the fix keeps chunk_N_incomplete.tsv as the list of samples to recompute.
        chunks = sorted((self.root / 'outputs' / 'batch_001' / 'chunks').glob('chunk_*_incomplete.tsv'))
        incomplete = [sid for path in chunks for sid in pd.read_csv(path, sep='\t')['srr_accession']]
        self.assertIn('sample_A', incomplete, 'sample with an unparseable guide JSON was treated as complete')

    def test_C21_failed_prodigal_worker_fails_the_run(self):
        """When any Prodigal worker fails, is110_run_all.sh exits non-zero instead of finishing as Done."""
        chunk = self.root / 'outputs' / 'batch_006' / 'chunks' / 'chunk_0' / 'sniffles_output'
        self.sniffles_table(chunk, 'S', ['good1', 'failme1'])
        proc = self.run_is110_all()
        chunks = self.is110_out / 'prodigal_chunks'
        called = sorted(p.read_text() for p in chunks.glob('*.called'))
        self.assertEqual(called, ['failme1', 'good1'], proc.stderr)  # both workers ran; one of them exited 7
        self.assertNotEqual(proc.returncode, 0,
                            'is110_run_all.sh exited 0 although a Prodigal worker exited 7; stdout ends: '
                            + ' | '.join(proc.stdout.splitlines()[-4:]))

    def test_C22_stale_prodigal_chunk_is_not_merged(self):
        """Proteins from a chunk file left by an earlier run are not part of this run's protein set."""
        chunk = self.root / 'outputs' / 'batch_006' / 'chunks' / 'chunk_0' / 'sniffles_output'
        self.sniffles_table(chunk, 'S', ['good1'])
        old = self.is110_out / 'prodigal_chunks'
        old.mkdir(parents=True)
        (old / 'chunk_47.faa').write_text('>stale_1 # 1 # 300 # 1 # ID=1_1;partial=00\nMKV*\n')
        proc = self.run_is110_all()
        # Assumes the merged protein set stays at <out>/all_proteins.faa.
        proteins = self.is110_out / 'all_proteins.faa'
        self.assertTrue(proteins.is_file(), proc.stderr)
        names = [line[1:].split()[0] for line in proteins.read_text().splitlines() if line.startswith('>')]
        self.assertIn('good1_1', names, proc.stderr)  # the current input was processed
        self.assertNotIn('stale_1', names, 'protein from a previous run was merged into all_proteins.faa')

    def test_C23_existing_top_level_sniffles_dir_does_not_hide_chunk_results(self):
        """Chunk-local Sniffles calls are collected even when the batch also has a top-level sniffles_output."""
        batch = self.root / 'outputs' / 'batch_006'
        (batch / 'sniffles_output').mkdir(parents=True)  # legacy directory, empty
        self.sniffles_table(batch / 'chunks' / 'chunk_0' / 'sniffles_output', 'S', ['u'])
        with self.subTest('is110_build_json.py'):
            self.assertIn('u', self.build_is110_json(), 'chunk-local insertion u is missing from is110_all_records.json')
        with self.subTest('is110_run_all.sh step 1'):
            proc = self.run_is110_all()
            consensus = self.is110_out / 'all_consensus.fna'
            self.assertTrue(consensus.is_file(), proc.stderr)
            self.assertIn('>u\n', consensus.read_text(), 'chunk-local insertion u is missing from all_consensus.fna')

    def test_C24_top_level_circle_evidence_is_kept(self):
        """Tail-head reads recorded in a top-level circle_output summary are not reported as no circle evidence."""
        batch = self.root / 'outputs' / 'batch_006'
        self.sniffles_table(batch / 'sniffles_output', 'S', ['u'])
        write_tsv(batch / 'circle_output' / 'S' / 'S_circle_summary.tsv', [{
            'is_uuid': 'u', 'n_tail_head_reads': 3, 'n_genome_head_reads': 1,
            'n_tail_genome_reads': 1, 'n_total_mapped': 9}])
        records = self.build_is110_json()
        self.assertIn('u', records)
        self.assertNotEqual(records['u'].get('circle_type'), 'without_circle_evidence',
                            f"element with 3 tail-head reads exported as {records['u'].get('circle_evidence')}")

    def test_C25_symlinked_result_folders_are_uploaded(self):
        """Files in is_formatter_output, partial_circle_output and sniffles_output are part of the upload."""
        batch = self.root / 'outputs' / 'batch_099'
        files = {
            'is_formatter_output/S/S_is_records_guide.json': b'[{"is_id": "marker-guide"}]',
            'partial_circle_output/S/S_partial_circles.json': b'[{"is_id": "marker-partial"}]',
            'sniffles_output/S/S.table.txt': b'UUID\nmarker-sniffles\n',
            'chunks/chunk_0_all.tsv': b'srr_accession\nS\n',
        }
        for name, data in files.items():
            (batch / name).parent.mkdir(parents=True, exist_ok=True)
            (batch / name).write_bytes(data)
        sent = self.staged_upload(batch)
        self.assertIn(files['chunks/chunk_0_all.tsv'], sent.values())  # the replaced Hub call did list files
        # Compared by content, so the fix is free to choose the paths in the repo.
        missing = [name for name, data in files.items() if data not in sent.values()]
        self.assertEqual(missing, [], f'result files left out of the upload; sent only {sorted(sent)}')

    def test_C26_chunk_local_sniffles_calls_are_uploaded(self):
        """Sniffles tables and VCFs that exist only under chunks/chunk_N/sniffles_output are part of the upload."""
        batch = self.root / 'outputs' / 'batch_099'
        files = {
            'chunks/chunk_0/sniffles_output/S/S.table.txt': b'UUID\nmarker-chunk-table\n',
            'chunks/chunk_0/sniffles_output/S/S.vcf': b'##fileformat=VCFv4.2\nmarker-chunk-vcf\n',
            'chunks/chunk_0_all.tsv': b'srr_accession\nS\n',
        }
        for name, data in files.items():
            (batch / name).parent.mkdir(parents=True, exist_ok=True)
            (batch / name).write_bytes(data)
        sent = self.staged_upload(batch)
        self.assertIn(files['chunks/chunk_0_all.tsv'], sent.values())  # the replaced Hub call did list files
        # Assumes the calls are uploaded as they are or merged into one file, not re-encoded.
        blob = b'\n'.join(sent.values())
        missing = [m for m in (b'marker-chunk-table', b'marker-chunk-vcf') if m not in blob]
        self.assertEqual(missing, [], f'chunk-local Sniffles calls left out of the upload; sent only {sorted(sent)}')

    def test_p_a_pacbio_sample_is_not_aligned_with_nanopore_preset(self):
        """A PACBIO_SMRT row of the batch table is not aligned with minimap2 -x map-ont."""
        pipeline = load_script('run_pipeline')
        metadata = self.root / 'batch.tsv'
        write_tsv(metadata, [{'srr_accession': 'SRRPB1', 'organism': 'A', 'platform': 'PACBIO_SMRT',
                              'instrument': 'Sequel II', 'library_layout': 'SINGLE'}])
        fastq = self.root / 'fastq' / 'SRRPB1.fastq'
        fastq.parent.mkdir()
        fastq.write_text('@r\nACGT\n+\nIIII\n')
        fasta = self.root / 'ref' / 'GCF_1' / 'GCF_1_genomic.fna'
        fasta.parent.mkdir(parents=True)
        fasta.write_text('>chr\nACGT\n')
        write_tsv(self.root / 'ref' / 'resolve_status.tsv', [{'organism': 'A', 'accession': 'GCF_1', 'fasta': str(fasta)}])

        alignments = []
        def pipe(producer, consumer, **kwargs):
            alignments.append(producer)
            Path(consumer[consumer.index('-o') + 1]).write_bytes(b'bam')
        # Assumes read alignment stays in cycle.preprocess.aligner and uses run_checked_pipe.
        with patch('shutil.which', side_effect=fake_which), \
                patch('cycle.preprocess.aligner.run_checked_pipe', side_effect=pipe), \
                patch('cycle.preprocess.aligner.subprocess.run', side_effect=fake_run):
            run_main(pipeline, ['--metadata', metadata, '--fastq-dir', fastq.parent, '--ref-dir', self.root / 'ref',
                                '--align-dir', self.root / 'aln', '--steps', 'align'])
        self.assertEqual(len(alignments), 1, 'the sample was not aligned at all')
        self.assertEqual(alignments[0][0], 'minimap2')
        self.assertNotIn('map-ont', alignments[0], f"PACBIO_SMRT sample aligned with: {' '.join(alignments[0])}")

    def test_p_b_one_bad_sample_does_not_strand_the_batch(self):
        """A sample without a FASTQ neither stops its chunk's good sample before formatting nor blocks later phases."""
        write_tsv(self.root / 'batch_metadata' / 'batch_001.tsv', [
            {'srr_accession': 'sample_GOOD', 'organism': 'A'}, {'srr_accession': 'sample_BAD', 'organism': 'A'}])
        proc, jobs = self.submit('1')
        self.assertEqual(proc.returncode, 0, proc.stderr)
        phase1 = [job for job in jobs if '_p1_' in job['name']]
        self.assertEqual(len(phase1), 1, [job['name'] for job in jobs])

        # Replay the Phase 1 job as the shell would: run_pipeline.py calls joined by &&, stopping at
        # the first non-zero exit. Assumes the fix keeps that form of the job script.
        commands = []
        for part in phase1[0]['wrap'].replace('\\\n', '').split('&&'):
            tokens = shlex.split(part)
            if len(tokens) > 1 and tokens[0] == 'python' and tokens[1].endswith('run_pipeline.py'):
                commands.append(tokens[2:])
        self.assertGreaterEqual(len(commands), 2, phase1[0]['wrap'])

        pipeline = load_script('run_pipeline')
        def downloader(output_dir, **kwargs):
            def download_batch(frame, **kwargs):
                Path(output_dir).mkdir(parents=True, exist_ok=True)
                ok = frame['srr_accession'] != 'sample_BAD'  # the download of sample_BAD fails
                for sid in frame['srr_accession'][ok]:
                    (Path(output_dir) / f'{sid}.fastq').write_text('@r\nACGT\n+\nIIII\n')
                return frame.assign(fastq_path='', download_status=['ok' if v else 'failed' for v in ok])
            return MagicMock(download_batch=download_batch)
        def resolver(output_dir, **kwargs):
            def resolve_all(frame):
                fasta = Path(output_dir) / 'A.fna'
                fasta.parent.mkdir(parents=True, exist_ok=True)
                fasta.touch()
                return {org: {'accession': org, 'fasta': fasta} for org in frame['organism']}
            return MagicMock(resolve_all=resolve_all)
        aligner = MagicMock()
        aligner.align_batch.side_effect = lambda samples: [{**s, 'bam': self.root / 'ok.bam'} for s in samples]
        sniffles = MagicMock()
        sniffles.run_batch.side_effect = lambda metadata, ref_map, **kwargs: {sid: self.root / 'table' for sid in ref_map}
        finder = MagicMock()
        finder.run_batch.side_effect = lambda tldr_results, **kwargs: {
            sid: self.root / 'circle' for sid, table in tldr_results.items() if table}
        formatted = set()
        def format_batch(circle_results, **kwargs):
            formatted.update(sid for sid, found in circle_results.items() if found)
            return {sid: self.root / 'records' for sid in circle_results}
        formatter = MagicMock()
        formatter.run_batch.side_effect = format_batch

        job_failed = False
        with patch.object(pipeline, 'SRADownloader', side_effect=downloader), \
                patch.object(pipeline, 'ReferenceGenomeResolver', side_effect=resolver), \
                patch.object(pipeline, 'Aligner', side_effect=aligner_factory(aligner)), \
                patch.object(pipeline, 'SnifflesRunner', return_value=sniffles), \
                patch.object(pipeline, 'CircleFinder', return_value=finder), \
                patch.object(pipeline, 'ISFormatter', return_value=formatter):
            for argv in commands:
                if run_main(pipeline, argv) != 0:
                    job_failed = True
                    break

        with self.subTest('good sample in the same chunk'):
            self.assertIn('sample_GOOD', formatted, 'sample_GOOD never reached the format stage of its chunk')
        with self.subTest('later phases'):
            # A job is stranded when it waits with afterok on a job that failed or is itself stranded.
            dead = {phase1[0]['id']} if job_failed else set()
            stranded = []
            for job in jobs:
                for spec in filter(None, job['dependency'].split(',')):
                    kind, *ids = spec.split(':')
                    if kind == 'afterok' and dead & set(ids):
                        dead.add(job['id'])
                        stranded.append(job['name'])
                        break
            self.assertEqual(stranded, [], 'jobs that can never start after one sample lacked a FASTQ')

    def test_p_c_bam_from_another_reference_is_not_reused(self):
        """An existing BAM is not returned for a sample whose reference differs from the one it was aligned to."""
        fastq = self.root / 'reads.fastq'
        fastq.write_text('@r\nACGT\n+\nIIII\n')
        old_ref, new_ref = self.root / 'old_ref.fna', self.root / 'new_ref.fna'
        for ref in (old_ref, new_ref):
            ref.write_text('>chr\nACGT\n')
        alignments = []
        def pipe(producer, consumer, **kwargs):
            alignments.append(producer)
            Path(consumer[consumer.index('-o') + 1]).write_bytes(b'bam')
        with patch('shutil.which', side_effect=fake_which), \
                patch('cycle.preprocess.aligner.run_checked_pipe', side_effect=pipe), \
                patch('cycle.preprocess.aligner.subprocess.run', side_effect=fake_run):
            aligner = Aligner(output_dir=str(self.root / 'aln'))
            first = aligner.align(fastq, old_ref, 'sample')
            self.assertTrue(first and first.exists() and len(alignments) == 1)
            try:
                second = aligner.align(fastq, new_ref, 'sample')
            except Exception:
                return  # refusing the stale BAM with an error is also correct
        # Either the sample is realigned to the new reference (FASTA or its .mmi) or no BAM is returned.
        realigned = any(str(new_ref.with_suffix('')) in arg for cmd in alignments[1:] for arg in cmd)
        self.assertTrue(second is None or realigned,
                        f'{second} was aligned to old_ref.fna but returned for new_ref.fna without realignment')

    def test_p_d_leading_zero_batch_number_is_decimal(self):
        """submit_one_batch.sh 010 never submits batch 8."""
        write_tsv(self.root / 'batch_metadata' / 'batch_008.tsv', [{'srr_accession': 'sample_EIGHT', 'organism': 'A'}])
        write_tsv(self.root / 'batch_metadata' / 'batch_010.tsv', [{'srr_accession': 'sample_TEN', 'organism': 'A'}])
        proc, jobs = self.submit('010')
        # Running batch 10 and rejecting the argument are both correct; running batch 8 is not.
        wrong = [job['name'] for job in jobs if 'batch_008' in job['name'] or 'batch_008' in job['wrap']]
        self.assertEqual(wrong, [], 'argument 010 submitted jobs for batch_008: ' + proc.stdout[:120].replace('\n', ' '))
        self.assertFalse((self.root / 'outputs' / 'batch_008').exists(), 'argument 010 created batch_008 outputs')

    def test_p_e_truncated_or_paired_fastq_is_not_silently_used(self):
        """A truncated FASTQ is never handed to alignment, and a paired run is not reduced to its _1 file."""
        with self.subTest('truncated download'):
            pipeline = load_script('run_pipeline')
            metadata = self.root / 'batch.tsv'
            write_tsv(metadata, [{'srr_accession': 'SRRTR1', 'organism': 'A'}])
            fasta = self.root / 'ref' / 'GCF_1' / 'GCF_1_genomic.fna'
            fasta.parent.mkdir(parents=True)
            fasta.write_text('>chr\nACGT\n')
            write_tsv(self.root / 'ref' / 'resolve_status.tsv', [{'organism': 'A', 'accession': 'GCF_1', 'fasta': str(fasta)}])
            fastq_dir = self.root / 'fastq'
            truncated = fastq_dir / 'SRRTR1.fastq.gz'
            whole = gzip.compress(''.join(f'@r{i}\nACGTACGTAC\n+\nIIIIIIIIII\n' for i in range(5000)).encode())

            def run(cmd, **kwargs):
                if cmd[0] == 'kingfisher':  # reports success but leaves half a file
                    truncated.write_bytes(whole[:len(whole) // 2])
                return fake_run(cmd, **kwargs)
            aligned = []
            aligner = MagicMock()
            aligner.align_batch.side_effect = lambda samples: aligned.extend(samples) or [
                {**s, 'bam': self.root / 'ok.bam'} for s in samples]
            with patch('shutil.which', side_effect=fake_which), \
                    patch('cycle.download_manager.sra_download.subprocess.run', side_effect=run), \
                    patch.object(pipeline, 'Aligner', side_effect=aligner_factory(aligner)):
                run_main(pipeline, ['--metadata', metadata, '--fastq-dir', fastq_dir, '--ref-dir', self.root / 'ref',
                                    '--align-dir', self.root / 'aln', '--steps', 'download', 'align'])
            with self.assertRaises(EOFError):  # the file left by the download really is unreadable
                gzip.open(truncated).read()
            used = [str(s['fastq']) for s in aligned]
            self.assertNotIn(str(truncated), used, 'a gzip file cut in half was passed to alignment')

        with self.subTest('paired run'):
            paired = self.root / 'paired'
            paired.mkdir()
            for mate in ('1', '2'):
                (paired / f'SRRPE1_{mate}.fastq.gz').write_bytes(gzip.compress(b'@r\nACGT\n+\nIIII\n'))
            try:
                found = find_fastq(paired, 'SRRPE1')
            except ValueError:
                found = None  # refusing, as find_fastq already does for other multi-file runs, is correct
            # Returning both files in any form is also correct; returning only _1 is not.
            self.assertNotEqual(found, paired / 'SRRPE1_1.fastq.gz', 'SRRPE1_2.fastq.gz was dropped without notice')


if __name__ == '__main__':
    unittest.main()
