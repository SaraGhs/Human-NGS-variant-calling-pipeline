import sys
import unittest
from pathlib import Path
from importlib.util import spec_from_file_location, module_from_spec
SCRIPT=Path(__file__).resolve().parents[1]/'scripts/09_compare_fastqc.py'
spec=spec_from_file_location('fastqc_comp',SCRIPT)
mod=module_from_spec(spec);spec.loader.exec_module(mod)

class TestFastQCComparison(unittest.TestCase):
    def test_keys(self):
        stem='H06HDADXX130110.1.ATCACGAT.20k'
        self.assertEqual(mod.sample_key(stem+'_reads_1.fastq'),(stem,'1'))
        self.assertEqual(mod.sample_key(stem+'_R2_paired.fastq.gz',trimmed=True),(stem,'2'))
    def test_bad_name(self):
        with self.assertRaises(ValueError):mod.sample_key('my_sample.fastq')
    def test_counts(self):
        records=[{'raw_adapter_content':'warn','trimmed_adapter_content':'pass'}, {'raw_adapter_content':'pass','trimmed_adapter_content':'pass'}]
        self.assertEqual(mod.count_flags(records,'adapter_content','trimmed')['pass'],2)
    def test_mismatch(self):
        with self.assertRaises(ValueError):mod.merge({('a','1'):{}},{('b','1'):{}})

if __name__=='__main__':unittest.main()
