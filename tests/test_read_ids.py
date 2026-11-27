"""Records name the reads they rest on, so the reads can be fetched again."""
import unittest

from cycle.circle_detect.circle_finder import CircleFinder
from cycle.is_formatter.formatter import ISFormatter


class ReadIds(unittest.TestCase):
    def test_query_position_of_the_junction(self):
        # 10 clipped bases, 50 aligned from reference 20: reference 40 is read base 30.
        self.assertEqual(CircleFinder._query_pos([(4, 10), (0, 50)], 20, 40), 30)
        # A deletion over the junction: the next read base after it.
        self.assertEqual(CircleFinder._query_pos([(0, 20), (2, 4), (0, 20)], 20, 41), 20)
        # An insertion before the junction shifts the read position.
        self.assertEqual(CircleFinder._query_pos([(0, 10), (1, 3), (0, 30)], 0, 20), 23)

    def test_record_carries_read_ids_and_junction(self):
        formatter = object.__new__(ISFormatter)
        record = formatter._format_record(
            {'is_uuid': 'u', 'chrom': 'c', 'start': 5, 'end': 5,
             'n_tail_head_reads': 2, 'n_genome_head_reads': 0, 'n_tail_genome_reads': 0,
             'example_th_read_id': 'SRR1.7', 'junction_seq': 'ACGT', 'target_duplication': '4',
             'tail_head_read_ids': ['SRR1.7', 'SRR1.9']},
            'AAAA', '', '', 'S', 'polished_read', 'SRR1.3', 12)
        self.assertEqual(record['sequence_read_id'], 'SRR1.3')
        self.assertEqual(record['n_element_reads'], 12)
        evidence = record['circle_evidence']
        self.assertEqual(evidence['tail_head_read_ids'], ['SRR1.7', 'SRR1.9'])
        self.assertEqual(evidence['example_th_read_id'], 'SRR1.7')
        self.assertEqual(evidence['junction_sequence'], 'ACGT')
        self.assertEqual(evidence['target_duplication'], '4')


if __name__ == '__main__':
    unittest.main()
