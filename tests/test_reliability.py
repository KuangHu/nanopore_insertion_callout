"""Regression tests for sample identity, pipeline exits and resumable jobs."""

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

from cycle.processes import run_checked_pipe
from cycle.utils import find_fastq
from cycle.partial_circle_state import input_signature, is_complete, mark_complete, output_paths
from cycle.preprocess.aligner import Aligner
from cycle.is_formatter.formatter import ISFormatter
from cycle.circle_detect.circle_finder import CircleFinder
from cycle.partial_circle.detector import PartialCircleDetector

ROOT = Path(__file__).resolve().parents[1]


def load_script(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / 'scripts' / f'{name}.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class TemporaryTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)


class TestFastq(TemporaryTest):
    def test_neighbor_accession_is_not_selected(self):
        (self.root / 'SRR1234.fastq.gz').touch()
        self.assertIsNone(find_fastq(self.root, 'SRR123'))

    def test_string_path_and_exact_file(self):
        expected = self.root / 'SRR123.fastq.gz'
        expected.touch()
        self.assertEqual(find_fastq(str(self.root), 'SRR123'), expected)

    def test_unfinished_download_and_directory_are_not_inputs(self):
        (self.root / 'SRR123.fastq.gz.part').touch()
        (self.root / 'SRR123.fastq.gz').mkdir()
        self.assertIsNone(find_fastq(self.root, 'SRR123'))

    def test_delimited_suffix_is_accepted(self):
        expected = self.root / 'SRR123_custom.fq.gz'
        expected.touch()
        self.assertEqual(find_fastq(self.root, 'SRR123'), expected)

    def test_multiple_fallback_shards_are_not_silently_truncated(self):
        for suffix in ('laneA', 'laneB'):
            (self.root / f'SRR123_{suffix}.fastq.gz').touch()
        with self.assertRaises(ValueError):
            find_fastq(self.root, 'SRR123')


class TestCheckedPipe(unittest.TestCase):
    def test_success_and_large_stderr_do_not_deadlock(self):
        run_checked_pipe(
            [sys.executable, '-c', "import sys; sys.stderr.write('x'*200000); print('ok')"],
            [sys.executable, '-c', "import sys; assert sys.stdin.read().strip() == 'ok'"],
            timeout=10,
        )

    def test_producer_failure_after_output_is_not_success(self):
        with self.assertRaises(subprocess.CalledProcessError) as raised:
            run_checked_pipe(
                [sys.executable, '-c', "import sys; print('partial'); sys.exit(9)"],
                [sys.executable, '-c', 'import sys; sys.stdin.read()'], timeout=10,
            )
        self.assertEqual(raised.exception.returncode, 9)

    def test_consumer_failure(self):
        with self.assertRaises(subprocess.CalledProcessError):
            run_checked_pipe(
                [sys.executable, '-c', "print('input')"],
                [sys.executable, '-c', 'import sys; sys.stdin.read(); sys.exit(3)'], timeout=10,
            )

    def test_timeout(self):
        with self.assertRaises(subprocess.TimeoutExpired):
            run_checked_pipe(
                [sys.executable, '-c', 'import time; time.sleep(10)'],
                [sys.executable, '-c', 'import sys; sys.stdin.read()'], timeout=0.1,
            )


class TestAlignmentPublication(TemporaryTest):
    def test_failed_producer_cannot_publish_or_index_any_mapping_stage(self):
        for cls, method, module in (
            (Aligner, 'align', 'cycle.preprocess.aligner'),
            (CircleFinder, '_map_reads', 'cycle.circle_detect.circle_finder'),
            (PartialCircleDetector, '_map_reads', 'cycle.partial_circle.detector'),
        ):
            with self.subTest(stage=module), patch('shutil.which', return_value='fake'):
                obj = cls(output_dir=str(self.root))
                output = self.root / f'{cls.__name__}.bam'
                def fail(producer, consumer, **kwargs):
                    Path(consumer[consumer.index('-o') + 1]).touch()
                    raise subprocess.CalledProcessError(9, producer)
                with patch(module + '.run_checked_pipe', side_effect=fail), patch(module + '.subprocess.run') as index:
                    if cls is Aligner:
                        result = obj.align(self.root / 'reads.fq', self.root / 'ref.fa', 'sample')
                    else:
                        result = getattr(obj, method)(self.root / 'reads.fq', self.root / 'ref.fa', output)
                self.assertIsNone(result)
                index.assert_not_called()
                self.assertFalse(list(self.root.glob('*.bam')))
                self.assertFalse(list(self.root.glob('.*.pending.bam')))

    def test_success_publishes_bam_and_index(self):
        with patch('shutil.which', return_value='fake'):
            obj = Aligner(output_dir=str(self.root))
        def align(producer, consumer, **kwargs):
            Path(consumer[consumer.index('-o') + 1]).write_bytes(b'complete')
        def index(cmd, **kwargs):
            Path(cmd[-1] + '.bai').write_bytes(b'index')
        with patch('cycle.preprocess.aligner.run_checked_pipe', side_effect=align), patch(
            'cycle.preprocess.aligner.subprocess.run', side_effect=index
        ):
            result = obj.align(self.root / 'reads.fq', self.root / 'ref.fa', 'sample')
        self.assertEqual(result.read_bytes(), b'complete')
        self.assertEqual(Path(str(result) + '.bai').read_bytes(), b'index')

    def test_index_failure_does_not_publish_bam(self):
        with patch('shutil.which', return_value='fake'):
            obj = Aligner(output_dir=str(self.root))
        def align(producer, consumer, **kwargs):
            Path(consumer[consumer.index('-o') + 1]).write_bytes(b'complete')
        with patch('cycle.preprocess.aligner.run_checked_pipe', side_effect=align), patch(
            'cycle.preprocess.aligner.subprocess.run', side_effect=subprocess.CalledProcessError(1, 'index')
        ):
            result = obj.align(self.root / 'reads.fq', self.root / 'ref.fa', 'sample')
        self.assertIsNone(result)
        self.assertFalse((self.root / 'sample.sorted.bam').exists())
        self.assertFalse(list(self.root.glob('.*.pending.bam')))


class TestCompleteBoundary(TemporaryTest):
    def locate(self, *alignments):
        obj = object.__new__(ISFormatter)
        obj.assembly_timeout = 1
        sam = '\n'.join('\t'.join([
            'query', str(flag), 'contig', '101', str(mapq), cigar, '*', '0', '0', '*', '*',
        ]) for flag, mapq, cigar in alignments)
        with patch('cycle.is_formatter.formatter.subprocess.run', return_value=MagicMock(
            returncode=0, stdout=sam.encode(),
        )):
            return obj._locate_is_in_contig('A' * 2000, 'A' * 1000, self.root)

    def test_clipped_alignment_cannot_define_full_element(self):
        for cigar in ('400S200M400S', '400H200M400H', '200M', '500M10N500M'):
            with self.subTest(cigar=cigar):
                self.assertIsNone(self.locate((0, 60, cigar)))

    def test_supplementary_does_not_override_complete_primary(self):
        self.assertEqual(self.locate((0, 20, '1000M'), (2048, 60, '400H200M400H')), (100, 1100, '+'))

    def test_full_reverse_and_internal_indel(self):
        self.assertEqual(self.locate((16, 60, '500M10I490M')), (100, 1090, '-'))

    def test_unknown_or_zero_mapping_confidence(self):
        for mapq in (0, 255):
            self.assertIsNone(self.locate((0, mapq, '1000M')))


class TestPartialCompletion(TemporaryTest):
    def setUp(self):
        super().setUp()
        self.sid = 'sample'
        self.records = self.root / 'records.json'
        self.records.write_text('[]')
        self.fastq = self.root / 'sample.fastq'
        self.fastq.touch()
        self.out = self.root / 'out' / self.sid
        self.out.mkdir(parents=True)
        self.signature = input_signature(self.records, self.fastq, {})

    def write_outputs(self):
        for path in output_paths(self.out, self.sid):
            path.write_text('[]' if path.suffix == '.json' else 'header\n')

    def test_directory_and_unmarked_outputs_do_not_skip(self):
        self.assertFalse(is_complete(self.out, self.sid, self.signature))
        self.write_outputs()
        self.assertFalse(is_complete(self.out, self.sid, self.signature))

    def test_success_empty_results_are_complete_but_changed_inputs_are_not(self):
        self.write_outputs()
        mark_complete(self.out, self.sid, self.signature)
        self.assertTrue(is_complete(self.out, self.sid, self.signature))
        self.records.write_text('[{}]')
        self.assertFalse(is_complete(self.out, self.sid, input_signature(self.records, self.fastq, {})))
        output_paths(self.out, self.sid)[0].unlink()
        self.assertFalse(is_complete(self.out, self.sid, self.signature))

    def test_invalid_summary_cannot_be_marked_complete(self):
        self.write_outputs()
        output_paths(self.out, self.sid)[-1].write_text('{broken')
        with self.assertRaises(ValueError):
            mark_complete(self.out, self.sid, self.signature)

    def test_failed_child_is_reported_and_retried_without_double_nesting(self):
        batch = load_script('run_partial_circle_batch')
        fmt = self.root / 'formatter' / self.sid
        fmt.mkdir(parents=True)
        (fmt / f'{self.sid}_is_records_guide.json').write_text('[]')
        calls = []
        def fail(cmd, **kwargs):
            calls.append(cmd)
            self.out.mkdir(exist_ok=True)
            raise subprocess.CalledProcessError(1, cmd)
        with patch.object(batch.subprocess, 'run', side_effect=fail):
            for _ in range(2):
                failures = batch.run_batch(fmt.parent, self.root, self.out.parent)
                self.assertEqual(failures, [self.sid])
        self.assertEqual(len(calls), 2)
        self.assertEqual(calls[0][calls[0].index('--output-dir') + 1], str(self.out.parent))

    def test_single_sample_rebuilds_incomplete_cache_then_skips_success(self):
        single = load_script('run_partial_circle')
        (self.out / f'{self.sid}_is_ref.fa').write_text('stale')
        (self.out / f'{self.sid}.partial.nsorted.bam').write_text('incomplete')
        detector = MagicMock()
        def run(**kwargs):
            self.assertFalse((self.out / f'{self.sid}_is_ref.fa').exists())
            self.assertFalse((self.out / f'{self.sid}.partial.nsorted.bam').exists())
            self.write_outputs()
            return self.out
        detector.run_sample.side_effect = run
        argv = ['run_partial_circle.py', '--is-records', str(self.records),
                '--fastq', str(self.fastq), '--sample-id', self.sid, '--output-dir', str(self.out.parent)]
        with patch.object(sys, 'argv', argv), patch.object(single, 'PartialCircleDetector', return_value=detector) as factory:
            single.main()
            single.main()
        self.assertEqual(factory.call_count, 1)


class TestChunkReferences(TemporaryTest):
    def test_actual_submission_commands_keep_interleaved_resolve_jobs_separate(self):
        bash = os.environ.get('BASH') or shutil.which('bash')
        if not bash:
            self.skipTest('bash is required to test generated SLURM commands')
        batch_dir = self.root / 'batch_metadata'
        batch_dir.mkdir()
        metadata = batch_dir / 'batch_001.tsv'
        pd.DataFrame([{'srr_accession':'sample_A', 'organism':'A'},
                      {'srr_accession':'sample_B', 'organism':'B'}]).to_csv(metadata, sep='\t', index=False)
        script = (ROOT / 'submit_one_batch.sh').read_text(encoding='utf-8')
        script = script.replace('CYCLE_DIR=/home/kuangh/tools/Cycle', 'CYCLE_DIR=' + shlex.quote(ROOT.as_posix()))
        script = script.replace('OUTROOT=/groups/rubin/projects/kuang/out/IS_cycle', 'OUTROOT=' + shlex.quote((self.root / 'outputs').as_posix()))
        script = script.replace('METADATA=$CYCLE_DIR/data/batches/${batch_name}.tsv', 'METADATA=' + shlex.quote(metadata.as_posix()))
        script = script.replace('LOGDIR=$HOME/logs', 'LOGDIR=' + shlex.quote((self.root / 'logs').as_posix()))
        script_path = self.root / 'submit.sh'
        script_path.write_text(script, encoding='utf-8', newline='\n')
        capture = self.root / 'args.bin'
        harness = ('python3() { ' + shlex.quote(Path(sys.executable).as_posix()) + ' "$@"; }; export -f python3\n'
                   'sbatch() { printf "%s\\0" "$@" >> ' + shlex.quote(capture.as_posix()) + '; echo 101; }; export -f sbatch\n'
                   + shlex.quote(Path(bash).as_posix()) + ' ' + shlex.quote(script_path.as_posix()) + ' 1 --nodes 2\n')
        result = subprocess.run([bash, '-c', harness], capture_output=True, text=True, encoding='utf-8')
        self.assertEqual(result.returncode, 0, result.stderr)
        args = capture.read_bytes().decode().split('\0')
        wraps = [a[len('--wrap='):] for a in args if a.startswith('--wrap=') and '--steps download resolve index' in a]
        pipeline = load_script('run_pipeline')
        calls = []
        aligner = MagicMock()
        aligner.align_batch.side_effect = lambda samples: calls.extend(samples) or [{**s, 'bam':self.root / 'ok.bam'} for s in samples]
        resolve_commands, align_commands = [], []
        for wrap in wraps:
            for part in wrap.replace('\\\n', '').split('&&'):
                tokens = shlex.split(part)
                if '--steps' not in tokens:
                    continue
                if 'resolve' in tokens:
                    begin = tokens.index('--steps')
                    end = tokens.index('--threads')
                    tokens[begin+1:end] = ['resolve']
                    resolve_commands.append(tokens[1:])
                elif 'align' in tokens:
                    align_commands.append(tokens[1:])
        self.assertEqual(len(resolve_commands), 2)
        def resolver_factory(output_dir, **kwargs):
            root = Path(output_dir)
            root.mkdir(parents=True, exist_ok=True)
            def resolve(frame):
                mapping = {}
                for _, row in frame.iterrows():
                    fasta = root / f'{row.organism}.fna'
                    fasta.touch()
                    mapping[row.organism] = {'accession':row.organism, 'fasta':fasta}
                return mapping
            obj = MagicMock()
            obj.resolve_all.side_effect = resolve
            return obj
        with patch.object(pipeline, 'ReferenceGenomeResolver', side_effect=resolver_factory), patch.object(pipeline, 'Aligner', return_value=aligner):
            # Both resolves complete before either align reads status from disk.
            for argv in resolve_commands:
                with patch.object(sys, 'argv', argv):
                    pipeline.main()
            for argv in align_commands:
                fastq = Path(argv[argv.index('--fastq-dir')+1])
                fastq.mkdir(parents=True, exist_ok=True)
                for sid in ('sample_A', 'sample_B'):
                    (fastq / f'{sid}.fastq').touch()
                Path(argv[argv.index('--align-dir')+1]).mkdir(parents=True, exist_ok=True)
                with patch.object(sys, 'argv', argv):
                    pipeline.main()
        self.assertEqual({s['sample_id'] for s in calls}, {'sample_A', 'sample_B'})

    def test_missing_alignment_inputs_exit_nonzero(self):
        pipeline = load_script('run_pipeline')
        metadata = self.root / 'samples.tsv'
        pd.DataFrame([{'srr_accession':'sample_A', 'organism':'A'}]).to_csv(metadata, sep='\t', index=False)
        argv = ['run_pipeline.py', '--metadata', str(metadata), '--outdir', str(self.root), '--steps', 'align']
        with patch.object(sys, 'argv', argv), patch.object(pipeline, 'Aligner'), self.assertRaises(SystemExit) as raised:
            pipeline.main()
        self.assertEqual(raised.exception.code, 1)


if __name__ == '__main__':
    unittest.main()
