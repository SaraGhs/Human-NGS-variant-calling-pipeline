"""Small synthetic fixtures; not results of the user's actual FASTQ analysis."""
import importlib.util,tempfile,unittest
from pathlib import Path
s=Path(__file__).resolve().parents[1]/'scripts'/'05_validate_fastq_pairs.py'
spec=importlib.util.spec_from_file_location('fastq_pair_validator',s)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
class Tests(unittest.TestCase):
    def test_good_pair(self):
        with tempfile.TemporaryDirectory() as t:
            a,b=Path(t)/'x_reads_1.fastq',Path(t)/'x_reads_2.fastq'
            a.write_text('@r001/1\nACGT\n+\nIIII\n@r002/1\nNN\n+\n##\n')
            b.write_text('@r001/2\nTGCA\n+\nIIII\n@r002/2\nNN\n+\n##\n')
            self.assertEqual(m.examine_pair(a,b)[0],2)
    def test_mismatched_id(self):
        with tempfile.TemporaryDirectory() as t:
            a,b=Path(t)/'x_reads_1.fastq',Path(t)/'x_reads_2.fastq'
            a.write_text('@r001/1\nACGT\n+\nIIII\n')
            b.write_text('@r002/2\nACGT\n+\nIIII\n')
            with self.assertRaises(ValueError): m.examine_pair(a,b)
    def test_unequal_seq_qual(self):
        with tempfile.TemporaryDirectory() as t:
            a=Path(t)/'bad.fastq';a.write_text('@r001\nACGT\n+\nIII\n')
            with self.assertRaises(ValueError):list(m.record_iter(a))
if __name__=='__main__':unittest.main()
