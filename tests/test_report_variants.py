import tempfile
import unittest
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from importlib import import_module
report = import_module('03_summarize_variants')

class ReportTests(unittest.TestCase):
    def test_classes(self):
        self.assertEqual(report.classify_variant('A','G'),'SNP')
        self.assertEqual(report.classify_variant('A','AT'),'insertion')
        self.assertEqual(report.classify_variant('AT','A'),'deletion')

    def test_parse(self):
        vcf = ('##fileformat=VCFv4.2\n'
               '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA12878\n'
               'chr19\t10\t.\tA\tG\t49\tPASS\tDP=21\tGT:DP\t0/1:21\n'
               'chr19\t20\t.\tA\tAT\t20\tLowQUAL\tDP=7\tGT:DP\t0/1:7\n')
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'test.vcf'
            path.write_text(vcf)
            rows = report.parse_vcf(path)
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0]['genotype'], '0/1')
        self.assertEqual(rows[1]['type'], 'insertion')
        self.assertEqual(report.summarize(rows)['pass_or_unfiltered'], 1)

if __name__ == '__main__':
    unittest.main()
