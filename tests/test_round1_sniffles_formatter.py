"""Round 1 review regressions: Sniffles insertion table, IS formatter, IS110 collectors.

Each test fails on the reviewed code because of one finding in
reviews/2026-11-26_codex_round1.md (C1, C2, C3, C4, C16) or one of two follow-up findings,
and should pass once that finding is fixed. Inputs are synthetic and built in a temp directory.
Sniffles and the assembler are replaced at the subprocess boundary; `minimap2 -a` runs for real.
"""

import csv
import importlib.util
import json
import os
from pathlib import Path
import random
import re
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import pysam

from cycle.is_formatter.formatter import ISFormatter
from cycle.preprocess.sniffles_runner import SnifflesRunner

ROOT = Path(__file__).resolve().parents[1]
REAL_RUN = subprocess.run
RNG = random.Random(20261126)


def rand_seq(n):
    return ''.join(RNG.choice('ACGT') for _ in range(n))


def revcomp(seq):
    return seq[::-1].translate(str.maketrans('ACGT', 'TGCA'))


def vcf_text(contig_len, records):
    header = ('##fileformat=VCFv4.2\n'
              f'##contig=<ID=chr1,length={contig_len}>\n'
              '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n')
    return header + ''.join(
        f'chr1\t{pos}\tsv{i}\t{ref}\t{alt}\t60\tPASS\t{info}\tGT\t0/1\n'
        for i, (pos, ref, alt, info) in enumerate(records))


def fake_tools(vcf=None, sniffles_rc=0, unitig=None):
    """subprocess.run stand-in: fake sniffles / overlap / miniasm, everything else runs for real."""
    def run(cmd, *args, **kwargs):
        tool = Path(str(cmd[0])).name
        if tool == 'sniffles':
            Path(cmd[cmd.index('--vcf') + 1]).write_text(vcf)
            return subprocess.CompletedProcess(cmd, sniffles_rc, '', 'killed' if sniffles_rc else '')
        if tool == 'minimap2' and 'ava-ont' in cmd:
            return subprocess.CompletedProcess(cmd, 0, b'a\t1\t0\t1\t+\tb\t1\t0\t1\t1\t1\t0\n', b'')
        if tool == 'miniasm':
            return subprocess.CompletedProcess(cmd, 0, f'S\tutg000001l\t{unitig}\n'.encode(), b'')
        if tool == 'minipolish':
            return subprocess.CompletedProcess(cmd, 1, b'', b'not available in this test')
        return REAL_RUN(cmd, *args, **kwargs)
    return run


def write_bam(path, references, reads):
    """reads: (name, sequence, reference index, start, mapq, cigartuples)."""
    header = {'HD': {'VN': '1.6', 'SO': 'coordinate'},
              'SQ': [{'SN': name, 'LN': length} for name, length in references]}
    with pysam.AlignmentFile(str(path), 'wb', header=header) as bam:
        for name, seq, ref_id, start, mapq, cigar in reads:
            read = pysam.AlignedSegment()
            read.query_name = name
            read.query_sequence = seq
            read.flag = 0
            read.reference_id = ref_id
            read.reference_start = start
            read.mapping_quality = mapq
            read.cigartuples = cigar
            read.query_qualities = pysam.qualitystring_to_array('I' * len(seq))
            bam.write(read)
    pysam.index(str(path))


def read_table(path):
    if path is None:
        return []
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


class TemporaryTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)


class TestSnifflesTable(TemporaryTest):
    def run_sniffles(self, bam, vcf, sniffles_rc=0, unitig=None, **kwargs):
        with patch('cycle.preprocess.sniffles_runner.shutil.which', return_value='fake'):
            runner = SnifflesRunner(output_dir=str(self.root / 'sniffles'), alignment_dir=str(self.root))
        ref = self.root / 'ref.fa'
        ref.touch()
        with patch('cycle.preprocess.sniffles_runner.subprocess.run',
                   side_effect=fake_tools(vcf, sniffles_rc, unitig)):
            return runner.run_sample(bam=bam, ref_fasta=ref, sample_id='S', **kwargs)

    def test_c1_symbolic_ins_consensus_excludes_flanks(self):
        """A symbolic <INS> must be stored as the inserted bases only, without reference flank."""
        genome = rand_seq(6000)
        insertion = rand_seq(500)
        # Three reads span the insertion point 3000 with 1500 bp of reference on each side.
        spanning = genome[1500:3000] + insertion + genome[3000:4500]
        bam = self.root / 'S.sorted.bam'
        write_bam(bam, [('chr1', 6000)],
                  [(name, spanning, 0, 1500, 60, [(0, 1500), (1, 500), (0, 1500)]) for name in 'abc'])
        vcf = vcf_text(6000, [(3000, 'N', '<INS>', 'PRECISE;SVTYPE=INS;SVLEN=500;SUPPORT=3;RNAMES=a,b,c')])
        # The assembler returns what miniasm returns for such reads: one unitig, flank + insertion + flank.
        table = self.run_sniffles(bam, vcf, unitig=spanning, min_size=100, max_size=1000)
        for row in read_table(table):
            consensus = row['Consensus'].upper()
            for flank in (genome[2900:3000], genome[3000:3100]):
                self.assertFalse(
                    flank in consensus or revcomp(flank) in consensus,
                    f'Consensus ({len(consensus)} bp) of a declared 500 bp <INS> contains 100 bp of '
                    'reference flank: the whole assembly unitig was stored as the insertion')
            self.assertLessEqual(len(consensus), 1000, 'Consensus is longer than max_size')

    def test_c2_pure_insertion_is_a_point_with_its_length_in_its_own_column(self):
        """A pure insertion is exported as Start = End = POS, with the inserted length in a Length column."""
        insertion = rand_seq(500)
        vcf = vcf_text(1200, [(1000, 'A', 'A' + insertion, 'PRECISE;SVTYPE=INS;SVLEN=500;SUPPORT=5')])
        rows = read_table(self.run_sniffles(self.root / 'S.sorted.bam', vcf))
        self.assertEqual(len(rows), 1, 'the resolved 500 bp insertion should give one table row')
        self.assertEqual(
            (int(rows[0]['Start']), int(rows[0]['End'])), (1000, 1000),
            'POS+SVLEN was exported as a reference coordinate')
        self.assertEqual(int(rows[0].get('Length', -1)), 500, 'inserted length is not in a Length column')
        # Start/End go unchanged into the formatter record, as CircleFinder reads them today.
        formatter = object.__new__(ISFormatter)
        record = formatter._format_record(
            {'is_uuid': rows[0]['UUID'], 'chrom': rows[0]['Chrom'],
             'start': int(rows[0]['Start']), 'end': int(rows[0]['End']),
             'n_tail_head_reads': 1, 'n_genome_head_reads': 0, 'n_tail_genome_reads': 0},
            insertion, '', '', 'S', 'longest_read')
        self.assertEqual((record['insertion_site']['pos_5p'], record['insertion_site']['pos_3p']), (1000, 1000))
        self.assertEqual((record['is_element']['start'], record['is_element']['end']), (1000, 1000))
        self.assertEqual(record['is_element']['length'], 500)

    def test_c16_failed_sniffles_vcf_not_reused_on_resume(self):
        """A VCF left by a failed Sniffles run must not become an insertion table on the next run."""
        vcf = vcf_text(6000, [(1000, 'A', 'A' + rand_seq(500), 'PRECISE;SVTYPE=INS;SVLEN=500;SUPPORT=5')])
        bam = self.root / 'S.sorted.bam'
        first = self.run_sniffles(bam, vcf, sniffles_rc=7)
        self.assertIsNone(first, 'a Sniffles run that exits 7 is a failure')
        # Sniffles fails again on the retry, so no successful run ever produced this VCF.
        second = self.run_sniffles(bam, vcf, sniffles_rc=7)
        self.assertIsNone(
            second,
            f'resume returned {second} ({len(read_table(second))} insertion rows) built from the '
            'VCF of a failed Sniffles run')

    def test_new1_bnd_breakend_alt_not_written_as_consensus(self):
        """A BND breakend ALT such as N[contig_2:1234[ is not a sequence and must not be a Consensus."""
        vcf = vcf_text(6000, [(1000, 'N', 'N[contig_2:1234[',
                               'PRECISE;SVTYPE=BND;SUPPORT=5;CHR2=contig_2;RNAMES=a,b,c')])
        for row in read_table(self.run_sniffles(self.root / 'S.sorted.bam', vcf)):
            self.assertRegex(
                row['Consensus'], r'^[ACGTNacgtn]+$',
                'breakend notation was written to the insertion table as a consensus sequence')


class TestFormatterRecords(TemporaryTest):
    def format_sample(self, entries, reads, unitig=None, tail_head=1):
        """entries: (uuid, start, consensus). reads: (name, sequence, bait uuid, mapq)."""
        for tool in ('minimap2', 'miniasm'):
            if not shutil.which(tool):
                self.skipTest(f'{tool} is required')
        circle = self.root / 'circle' / 'S'
        circle.mkdir(parents=True)
        columns = ['is_uuid', 'chrom', 'start', 'end', 'family', 'subfamily', 'consensus_length',
                   'consensus', 'n_tail_head_reads', 'n_genome_head_reads', 'n_tail_genome_reads',
                   'n_total_mapped', 'example_th_read']
        with open(circle / 'S_circle_summary.tsv', 'w', newline='') as fh:
            writer = csv.writer(fh, delimiter='\t')
            writer.writerow(columns)
            for uuid, start, consensus in entries:
                writer.writerow([uuid, 'chr1', start, start, 'Sniffles2_INS', 'S', len(consensus),
                                 consensus, tail_head, 0, 0, len(reads), ''])
        # Baits are the consensus twice, named {uuid}__th__j{length}, as CircleFinder builds them.
        baits = [(f'{uuid}__th__j{len(consensus)}', 2 * len(consensus)) for uuid, _, consensus in entries]
        uuids = [uuid for uuid, _, _ in entries]
        write_bam(circle / 'S.circle.sorted.bam', baits,
                  [(name, seq, uuids.index(bait), 0, mapq, [(4, 80), (0, len(seq) - 160), (4, 80)])
                   for name, seq, bait, mapq in reads])
        formatter = ISFormatter(output_dir=str(self.root / 'formatter'), threads=1)
        with patch('cycle.is_formatter.formatter.subprocess.run', side_effect=fake_tools(unitig=unitig)):
            out = formatter.run_sample(circle_dir=circle, sample_id='S')
        if out is None:
            return {}
        return {rec['is_id']: rec for rec in json.loads((out / 'S_is_records.json').read_text())}

    def test_c3_mapq0_read_from_other_locus_gives_flagged_record_without_flanks(self):
        """A locus whose only read may belong to another copy keeps its record, with empty flanks and a flag."""
        consensus = rand_seq(400)
        up_a, down_a = rand_seq(80), rand_seq(80)
        # The same element sits at chr1:1000 (locus A) and chr1:5000 (locus B). The read comes from A;
        # the two baits are identical, so the mapper put it on B's bait with MAPQ 0. Nothing in the
        # input ties the read to locus B.
        records = self.format_sample(
            [('locusA', 1000, consensus), ('locusB', 5000, consensus)],
            [('read_from_locusA', up_a + consensus + down_a, 'locusB', 0)])
        self.assertIn('locusB', records, 'the record of the ambiguous locus was dropped')
        record = records['locusB']
        self.assertEqual(
            (record['flanking_upstream']['sequence'], record['flanking_downstream']['sequence']), ('', ''),
            f"locusB record (assembly_method={record['assembly_method']}) carries flanks taken from "
            'the MAPQ 0 read read_from_locusA')
        # Field name chosen here; the fix must use the same one.
        self.assertEqual(record.get('flank_status'), 'ambiguous')

    def test_fallback_uses_a_shorter_read_when_the_longest_lacks_the_element(self):
        """With no assembly, the record comes from the longest read that holds the whole element."""
        consensus = rand_seq(400)
        up, down = rand_seq(80), rand_seq(80)
        records = self.format_sample(
            [('locusD', 1000, consensus)],
            [('long_read_without_element', rand_seq(3000), 'locusD', 60),
             ('spanning_read', up + consensus + down, 'locusD', 60)])
        record = records['locusD']
        self.assertEqual(record['assembly_method'], 'longest_read')
        self.assertEqual(record['is_element']['sequence'], consensus)
        self.assertEqual((record['flanking_upstream']['sequence'], record['flanking_downstream']['sequence']),
                         (up, down))

    def test_fallback_read_is_polished_with_the_other_reads(self):
        """The backbone read's own errors are removed by polishing with the other reads."""
        if not shutil.which('racon'):
            self.skipTest('racon is required')
        consensus = rand_seq(600)
        up, down = rand_seq(80), rand_seq(80)
        left, right = rand_seq(700), rand_seq(700)
        molecule = left + up + consensus + down + right
        noisy = list(molecule)
        for position in range(720, len(molecule) - 720, 53):   # substitutions inside flanks and element
            noisy[position] = 'A' if noisy[position] != 'A' else 'C'
        # mismatches are kept away from the element ends so that the consensus still aligns end to end
        reads = [('backbone', ''.join(noisy), 'locusE', 60)]
        reads += [(f'clean{i}', molecule[20 * i:len(molecule) - 20 * i], 'locusE', 60) for i in range(1, 7)]
        record = self.format_sample([('locusE', 1000, consensus)], reads)['locusE']
        self.assertEqual(record['assembly_method'], 'polished_read')
        self.assertEqual(record['is_element']['sequence'], consensus)
        self.assertEqual((record['flanking_upstream']['sequence'], record['flanking_downstream']['sequence']),
                         (up, down))

    def test_insertion_without_tail_head_reads_gets_a_record(self):
        """Every insertion gets a record; tail-head evidence is a field, not a filter."""
        consensus = rand_seq(400)
        up, down = rand_seq(80), rand_seq(80)
        records = self.format_sample(
            [('locusF', 1000, consensus)], [('spanning_read', up + consensus + down, 'locusF', 60)], tail_head=0)
        self.assertEqual(records['locusF']['evidence'], 'none')
        self.assertEqual(records['locusF']['is_element']['sequence'], consensus)

    def test_c4_reported_strand_matches_sequence_orientation(self):
        """A record must not say strand + when its sequence is the reverse complement of the consensus."""
        consensus = rand_seq(400)
        up, down = rand_seq(80), rand_seq(80)
        molecule = up + consensus + down
        # Three reads in bait orientation; the assembler returns the unitig on the other strand,
        # so the real `minimap2 -a` localisation of the consensus in it is on the minus strand.
        records = self.format_sample(
            [('locusC', 1000, consensus)],
            [(name, molecule, 'locusC', 60) for name in 'abc'],
            unitig=revcomp(molecule))
        record = records.get('locusC')
        if record is None:
            return
        element = record['is_element']
        # Weakest form: any fix may normalise the sequence, report "-", or leave the strand unset.
        self.assertFalse(
            element['strand'] == '+' and element['sequence'].upper() == revcomp(consensus),
            'is_element.sequence is the reverse complement of the consensus but strand is "+": '
            'the minus-strand localisation was discarded')


class TestIS110Collectors(TemporaryTest):
    def test_new2_collector_keeps_sniffles_sequence_when_formatter_sequence_empty(self):
        """An IS110 whose formatter record has no sequence still gets its Sniffles consensus."""
        consensus = rand_seq(600)
        for hmm in ('is110_all_006_026', 'is110_all_006_052'):
            (self.root / hmm).mkdir()
            (self.root / hmm / 'is110_ids.json').write_text('["u-1"]')
            for name in ('DEDD_hits.tbl', 'Tnp20_hits.tbl', 'all_proteins.faa'):
                (self.root / hmm / name).touch()
        # Formatter record as written when neither assembly nor the longest read located the element.
        fmt = self.root / 'batch_006' / 'is_formatter_output' / 'S'
        fmt.mkdir(parents=True)
        (fmt / 'S_is_records_guide.json').write_text(json.dumps([{
            'is_id': 'u-1', 'sample_id': 'S', 'organism': '',
            'is_element': {'sequence': '', 'length': 0, 'contig': 'chr1', 'start': 1000, 'end': 1600,
                           'strand': '+'},
            'flanking_upstream': {'sequence': '', 'length': 0},
            'flanking_downstream': {'sequence': '', 'length': 0},
            'insertion_site': {'contig': 'chr1', 'pos_5p': 1000, 'pos_3p': 1600},
            'circle_evidence': {'n_tail_head_reads': 3, 'n_genome_head_reads': 0,
                                'n_tail_genome_reads': 0, 'example_th_read': ''},
            'assembly_method': 'none', 'orf_annotation': None, 'guide_hits': [],
            'guide_summary': {'n_hits': 0, 'best_length': 0, 'has_revcomp_hit': False},
        }]))
        sniffles = self.root / 'batch_006' / 'sniffles_output' / 'S'
        sniffles.mkdir(parents=True)
        (sniffles / 'S.table.txt').write_text(
            'UUID\tChrom\tStart\tEnd\tFamily\tSubfamily\tConsensus\tSupport\n'
            f'u-1\tchr1\t1000\t1600\tSniffles2_INS\tS_chr1_1000\t{consensus}\t5\n')
        env = dict(os.environ, CYCLE_DATA_ROOT=str(self.root))
        for script, hmm in (('is110_build_json', 'is110_all_006_052'),
                            ('is110_collect_all', 'is110_all_006_026'),
                            ('is110_collect_simple', 'is110_all_006_026')):
            with self.subTest(script=script):
                if script == 'is110_collect_simple' and not importlib.util.find_spec('matplotlib'):
                    self.skipTest('is110_collect_simple imports the visualizer, which needs matplotlib')
                result = subprocess.run([sys.executable, str(ROOT / 'scripts' / f'{script}.py')],
                                        env=env, capture_output=True, text=True, timeout=300)
                self.assertEqual(result.returncode, 0, result.stderr[-2000:])
                records = [rec for path in sorted((self.root / hmm).rglob('is110*records.json'))
                           for rec in json.loads(path.read_text()) if rec['is_id'] == 'u-1']
                for path in (self.root / hmm).rglob('is110*records.json'):
                    path.unlink()
                self.assertTrue(records, 'u-1 is missing from the collector output')
                for rec in records:
                    self.assertTrue(
                        rec['is_element']['sequence'],
                        f'{script}: u-1 exported with an empty sequence (source '
                        f"{rec.get('_source', 'formatter')}) although its Sniffles table has a "
                        '600 bp consensus')


if __name__ == '__main__':
    unittest.main()
