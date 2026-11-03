"""Synthetic unit check of FASTQ pair counter."""
import importlib.util, tempfile, unittest, gzip
from pathlib import Path
s=Path(__file__).resolve().parents[1]/'scripts'/'08_summarize_trimming.py'
spec=importlib.util.spec_from_file_location('trim_summary',s)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
class Tests(unittest.TestCase):
    def test_gzip_read_count(self):
        with tempfile.TemporaryDirectory() as t:
            f=Path(t)/'test.fastq.gz'
            with gzip.open(f,'wt') as out:
                out.write('@r1\nAC\n+\nII\n@r2\nGG\n+\n##\n')
            self.assertEqual(m.reads(f),2)
if __name__=='__main__':unittest.main()
