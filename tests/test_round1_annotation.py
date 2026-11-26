"""Round-1 review tests for nesting, IS110 filter, system clustering and novelty annotation.

Each test fails on the code as reviewed (commit 93cda7b) because of one finding and should pass
once that finding is fixed. All inputs are synthetic; external tools are mocked.
"""

import json
import logging
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from cycle.is110_filter.filter import IS110Filter
from cycle.nesting_detector.detector import NestingDetector
from cycle.novelty_annotator import config as novelty_config
from cycle.novelty_annotator.annotator import NoveltyAnnotator
from cycle.system_clustering.clusterer import SystemClusterer


class TemporaryTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)


def write_records(formatter_dir, sample, kind, records):
    sample_dir = Path(formatter_dir) / sample
    sample_dir.mkdir(parents=True, exist_ok=True)
    (sample_dir / f'{sample}_is_records_{kind}.json').write_text(json.dumps(records))


class TestNestingDetector(TemporaryTest):
    def test_C12_minus_strand_cigar_blocks_pair_the_aligned_bases(self):
        """Sub-blocks of a minus-strand PAF record pair each target interval with the query bases aligned to it."""
        # PAF: query q (1500 bp) [100,1100) on the minus strand against target t (1300 bp) [0,1200).
        # The CIGAR runs along the target forward strand, so on a minus-strand record it walks the
        # query interval from its END: target [0,300) is aligned to query [800,1100) and, after the
        # 200 bp deletion, target [500,1200) is aligned to query [100,800).
        paf = self.root / 'minus.paf'
        paf.write_text('\t'.join([
            'q', '1500', '100', '1100', '-', 't', '1300', '0', '1200', '1000', '1200', '60',
            'cg:Z:300=200D700=',
        ]) + '\n')
        with patch('shutil.which', return_value='fake'):
            detector = NestingDetector(output_dir=str(self.root / 'out'))
        blocks = detector._parse_paf(paf)[('q', 't')]
        got = {(b.query_start, b.query_end, b.target_start, b.target_end) for b in blocks}

        # Design assumption: the fix has not chosen a coordinate frame yet. Accepted here are the
        # PAF frame (both forward, what the no-CIGAR branch already stores) and the two frames in
        # which one sequence is reverse-complemented (x -> length - x). The current output,
        # {(100, 400, 0, 300), (400, 1100, 500, 1200)}, is none of them.
        both_forward = {(800, 1100, 0, 300), (100, 800, 500, 1200)}
        query_flipped = {(400, 700, 0, 300), (700, 1400, 500, 1200)}
        target_flipped = {(800, 1100, 1000, 1300), (100, 800, 100, 800)}
        self.assertIn(
            got, (both_forward, query_flipped, target_flipped),
            'minus-strand sub-blocks pair target intervals with query bases that are not aligned to them',
        )


class TestIS110Filter(TemporaryTest):
    def test_C13_annotated_only_record_with_tail_head_reads_has_circle_evidence(self):
        """An IS110 record with tail-head reads counts as circle evidence even when it is only in the annotated JSON."""
        formatter_dir = self.root / 'is_formatter_output'
        write_records(formatter_dir, 'S1', 'annotated', [{
            'is_id': 'is-th', 'sample_id': 'S1',
            'circle_evidence': {'n_tail_head_reads': 3},
        }])
        with_circle, without_circle = IS110Filter().filter_records(str(formatter_dir), {'is-th'})
        self.assertEqual(
            ([r['is_id'] for r in with_circle], [r['is_id'] for r in without_circle]),
            (['is-th'], []),
            'record with 3 tail-head reads was not placed in the circle-evidence group',
        )

    def test_C14_orfs_separated_by_exactly_max_gap_bases_are_adjacent(self):
        """DEDD and Tnp20 ORFs with exactly max_orf_gap bases between them form an IS110 split pair."""
        # Prodigal coordinates are 1-based inclusive: [1,300] and [601,900] leave bases 301..600,
        # i.e. 300 bases, between the ORFs; [1,300] and [602,901] leave 301.
        hits = {
            'DEDD.hmm': ['at_limit__1_300_+', 'over_limit__1_300_+'],
            'Tnp20.hmm': ['at_limit__601_900_+', 'over_limit__602_901_+'],
        }

        def fake_hmmsearch(protein_fasta, output_tbl, hmm_path, cpus=8):
            Path(output_tbl).write_text(''.join(f'{p} - hmm - 1e-30 100.0 0.0\n' for p in hits[hmm_path]))
            return output_tbl

        filt = IS110Filter(dedd_hmm='DEDD.hmm', tnp20_hmm='Tnp20.hmm', max_orf_gap=300)
        with patch.object(filt, 'hmmsearch', side_effect=fake_hmmsearch):
            is110_ids, _ = filt.run_hmmsearch('unused.faa', str(self.root / 'hmm'))
        self.assertNotIn('over_limit', is110_ids, 'a 301-base gap must stay above max_orf_gap=300')
        self.assertIn('at_limit', is110_ids, 'a 300-base gap was rejected at max_orf_gap=300')

    def test_R2b_domain_hits_that_match_no_record_orf_are_not_dropped_silently(self):
        """Domain hits whose coordinates match no ORF of the record are still tagged, or the mismatch is reported."""
        # scripts/is110_run_all.sh calls the HMM proteins on the Sniffles consensus, while the ORFs
        # of a formatter record are called on is_element.sequence (a slice of the assembled contig
        # or of the longest read). Here the two differ by a 3 bp shift of the same ORF.
        protein = 'M' + 'A' * 332
        formatter_dir = self.root / 'is_formatter_output'
        write_records(formatter_dir, 'S1', 'guide', [{
            'is_id': 'is-a', 'sample_id': 'S1',
            'circle_evidence': {'n_tail_head_reads': 2},
            'orf_annotation': {'num_orfs': 1, 'orfs': [{
                'start': 104, 'end': 1105, 'strand': '+', 'length_nt': 1002, 'protein_sequence': protein,
            }]},
        }])
        consensus_orf = ('is-a__101_1102_+', 101, 1102, '+')
        trans_domains = {'is-a': {'DEDD': [consensus_orf], 'Tnp20': [consensus_orf]}}

        warnings = []
        handler = logging.Handler(level=logging.WARNING)
        handler.emit = warnings.append
        cycle_logger = logging.getLogger('cycle')
        cycle_logger.addHandler(handler)
        self.addCleanup(cycle_logger.removeHandler, handler)

        # Design assumption: the fix keeps the join inside IS110Filter and either matches the ORF
        # another way (overlap, protein sequence), raises, or logs a warning. A fix that only moves
        # the HMM search onto the record proteins, with no check here, would not satisfy this test.
        try:
            with_circle, _ = IS110Filter().filter_records(
                str(formatter_dir), {'is-a'}, trans_domains=trans_domains,
            )
        except Exception:
            return
        tagged = [o.get('domains') for o in with_circle[0]['orf_annotation']['orfs']]
        self.assertTrue(
            any(tagged) or warnings,
            'IS110 record came back with no domain-tagged ORF and no error or warning',
        )


class TestSystemClustering(TemporaryTest):
    def test_C19_forced_rerun_assigns_families_to_the_current_proteins(self):
        """A forced rerun clusters the proteins extracted now, not the ones in an older MMseqs2 result."""
        output_dir = self.root / 'system_clustering'
        stale = output_dir / 'mmseqs_clusters' / 'clust_cluster.tsv'
        stale.parent.mkdir(parents=True)
        stale.write_text('old__1_300_+\told__1_300_+\n')

        # Two elements carrying the same single protein: they share their only family.
        formatter_dir = self.root / 'is_formatter_output'
        write_records(formatter_dir, 'S1', 'guide', [
            {
                'is_id': is_id, 'sample_id': 'S1',
                'orf_annotation': {'orfs': [{
                    'start': 1, 'end': 300, 'strand': '+', 'length_nt': 300, 'protein_sequence': 'M' + 'K' * 98 + '*',
                }]},
            }
            for is_id in ('newA', 'newB')
        ])

        def fake_mmseqs(cmd, **kwargs):
            # Stand-in for `mmseqs easy-cluster <fasta> <prefix> <tmp>`: one family holding every
            # protein of the input FASTA. Assumes the fix still clusters with easy-cluster.
            i = cmd.index('easy-cluster')
            headers = [line[1:] for line in Path(cmd[i + 1]).read_text().splitlines() if line.startswith('>')]
            Path(f'{cmd[i + 2]}_cluster.tsv').write_text(''.join(f'{headers[0]}\t{h}\n' for h in headers))
            return subprocess.CompletedProcess(cmd, 0, '', '')

        with patch('shutil.which', return_value='fake'), patch('subprocess.run', side_effect=fake_mmseqs):
            results = SystemClusterer(output_dir).run([formatter_dir], skip_existing=False)
        clusters = sorted(sorted(c['members']) for c in json.loads(Path(results).read_text()))
        self.assertEqual(
            clusters, [['newA', 'newB']],
            'elements sharing their only protein family were not clustered together after a forced rerun',
        )


class TestNoveltyAnnotator(TemporaryTest):
    def annotator(self):
        with patch('shutil.which', return_value='fake'):
            return NoveltyAnnotator(output_dir=self.root / 'novelty')

    def test_C20_blast_database_is_rebuilt_for_a_different_isfinder_fasta(self):
        """A BLAST database built from one ISfinder FASTA is not reused for a FASTA with other content."""
        old_fasta = self.root / 'old.fna'
        new_fasta = self.root / 'new.fna'
        old_fasta.write_text('>ISold_IS3_IS3\n' + 'ACGT' * 50 + '\n')
        new_fasta.write_text('>ISnew_IS110_IS110\n' + 'TTGCA' * 50 + '\n')
        built_from = []

        def fake_makeblastdb(cmd, **kwargs):
            if cmd[0] == 'makeblastdb':
                built_from.append(cmd[cmd.index('-in') + 1])
                for suffix in ('.nhr', '.nin', '.nsq'):
                    Path(cmd[cmd.index('-out') + 1] + suffix).touch()
            return subprocess.CompletedProcess(cmd, 0, '', '')

        annotator = self.annotator()
        with patch('subprocess.run', side_effect=fake_makeblastdb):
            annotator._build_blastdb(old_fasta)
            annotator._build_blastdb(new_fasta)
        self.assertIn(
            str(new_fasta), built_from,
            'no BLAST database was built from the new ISfinder FASTA; the one built from old.fna was reused',
        )

    def test_R2a_blastn_task_can_find_hits_in_the_scored_identity_range(self):
        """The blastn command names a task that can find hits across the identity range the scoring uses."""
        # The divergence score is linear between PIDENT_CEILING and PIDENT_FLOOR (95 and 50 in
        # config.py; no 75 % threshold is written anywhere in the repo), so hits well below 95 %
        # identity have to be found. blastn's default task is megablast (word size 28).
        self.assertLess(novelty_config.PIDENT_FLOOR, 80)
        commands = []

        def fake_blastn(cmd, **kwargs):
            commands.append(cmd)
            return subprocess.CompletedProcess(cmd, 0, '', '')

        annotator = self.annotator()
        with patch('subprocess.run', side_effect=fake_blastn):
            annotator._run_blast(self.root / 'query.fna', self.root / 'isfinder_blastdb')
        cmd = [c for c in commands if c[0] == 'blastn'][0]
        task = cmd[cmd.index('-task') + 1] if '-task' in cmd else 'megablast (blastn default)'
        # Design assumption: sensitivity is set through -task. A fix that keeps megablast and only
        # lowers -word_size would need this assertion widened.
        self.assertIn(
            task, ('blastn', 'dc-megablast'),
            f'blastn task {task} cannot find hits near PIDENT_FLOOR={novelty_config.PIDENT_FLOOR}',
        )


if __name__ == '__main__':
    unittest.main()
