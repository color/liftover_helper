import collections
import os
import tempfile
import unittest

import vcf

from scripts import liftover


class RecordTest(unittest.TestCase):
    def test_convert_hg19_vcf_to_grch37_vcf(self):
        temp_dir = tempfile.gettempdir()
        output_file = os.path.join(temp_dir, 'output_grch37.vcf')
        liftover.convert_hg19_vcf_to_grch37_vcf('tests/hg19.vcf', output_file)

        with open(output_file) as fh:
            records = vcf.Reader(fh)
            expected_contigs = collections.OrderedDict([('1', liftover.contig_spec('1', 249250621)),
                                                        ('2', liftover.contig_spec('2', 243199373)),
                                                        ('3', liftover.contig_spec('3', 198022430))])
            self.assertEqual(records.contigs, expected_contigs)
            record = next(records)
            self.assertEqual(record.CHROM, '1')
            self.assertEqual(record.POS, 97915604)
            record = next(records)
            self.assertEqual(record.CHROM, '2')
            self.assertEqual(record.POS, 97915605)


    def test_record_overlaps_mismatch(self):
        records = list(vcf.Reader(filename='tests/grch38.vcf'))
        expected_result = 'chr2:21012603'
        for i in range(0, 5):
            self.assertEqual(liftover.find_overlapping_mismatch_site(records[i]), expected_result)
        for i in range(5, 7):
            self.assertFalse(liftover.find_overlapping_mismatch_site(records[i]))


    def test_update_record(self):
        records = list(vcf.Reader(filename='tests/grch38_2.vcf'))
        expected_results = [
            ['ATATG', 'ACATG,A', '1/2'],
            ['ATATG', 'A', '1/1'],
            ['AT', 'AC,A', '1/2'],
            ['AT', 'A', '1/1'],
            ['T', 'C,A', '1/2'],
            ['T', 'A', '1/1'],
            ['T', 'C', '0/1'],
            None,
            ['C', 'T', './.', True],
            ['C', 'T', '1/2', True],
            ['TATG', 'CATG,C', '1/2'],
            ['TATG', 'C', '1/1'],
            ['T', 'C,CAAT', '1/2'],
            ['T', 'CAAT', '0/1'],
            ['ATG', 'A', '0/1'],
            ['C', 'T', '0/1'],
        ]
        for i in range(len(records)):
            expected_record = expected_results[i]
            if expected_record and len(expected_record) == 4:
                self.assertRaises(ValueError, liftover.update_grch38_ref_to_grch37_for_record_if_needed, records[i])
                continue

            observed_record = liftover.update_grch38_ref_to_grch37_for_record_if_needed(records[i])
            if expected_record is None:
                self.assertIsNone(observed_record)
            else:
                self.assertEqual(observed_record.REF, expected_record[0])
                observed_alts = ','.join(map(str, observed_record.ALT))
                self.assertEqual(observed_alts, expected_record[1])
                self.assertEqual(observed_record.samples[0].data.GT, expected_record[2])

        # TODO: tests for is_anchor_base

if __name__ == '__main__':
    unittest.main()
