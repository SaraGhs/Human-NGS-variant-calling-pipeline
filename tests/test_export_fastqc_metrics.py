"""Synthetic unit tests for FastQC ZIP parsing; not real-data results."""
import importlib.util
from pathlib import Path
import tempfile
import unittest
import zipfile

SCRIPT=Path(__file__).resolve().parents[1]/'scripts'/'04_export_fastqc_metrics.py'
spec=importlib.util.spec_from_file_location('qc_parser',SCRIPT)
module=importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)

FIXTURE="""##FastQC\t0.13.0
>>Basic Statistics\tpass
#Measure\tValue
Filename\ttest.fastq
Total Sequences\t10000
Sequence length\t250
%GC\t41
>>END_MODULE
>>Per base sequence quality\tfail
#Base\tMean
1\t30
>>END_MODULE
>>Per base N content\tfail
#Base\tN-Count
1\t0.0
249\t20.0
250\t99.0
>>END_MODULE
>>Sequence Duplication Levels\tpass
#Total Deduplicated Percentage\t99.5
>>END_MODULE
>>Adapter Content\twarn
>>END_MODULE
"""

class TestExporter(unittest.TestCase):
    def test_basic_fields_and_status(self):
        r=module.parse_fastqc_text(FIXTURE)
        self.assertEqual(r['total_sequences'],'10000')
        self.assertEqual(r['gc_percent'],'41')
        self.assertEqual(r['max_per_base_n_percent'],'99.0')
        self.assertEqual(r['max_n_position'],'250')
        self.assertEqual(r['adapter_content'],'warn')
        self.assertEqual(r['deduplicated_percent'],'99.5')
    def test_archive_read(self):
        with tempfile.TemporaryDirectory() as temp:
            p=Path(temp)/'test_fastqc.zip'
            with zipfile.ZipFile(p,'w') as z:
                z.writestr('test_fastqc/fastqc_data.txt',FIXTURE)
            r=module.extract_one(p)
            self.assertEqual(r['filename'],'test.fastq')
            self.assertEqual(r['fastqc_zip'],'test_fastqc.zip')

if __name__=='__main__':unittest.main()
