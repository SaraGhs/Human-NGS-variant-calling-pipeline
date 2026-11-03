import importlib.util,unittest
from pathlib import Path
s=Path(__file__).resolve().parents[1]/'scripts'/'06_build_qc_summary.py'
spec=importlib.util.spec_from_file_location('qc_summary',s)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
class Tests(unittest.TestCase):
    def test_does_not_invent(self):
        rows=[{'filename':f'file{i}.fastq','total_sequences':'123','sequence_length':'151','gc_percent':'41','per_base_quality':'fail','per_sequence_quality':'warn','adapter_content':'pass','per_base_n_content':'fail'} for i in range(6)]
        pairs=[{'group':f'group{i}','read_pairs':'123','min_length_bp':'151','max_length_bp':'151','mean_length_bp':'151'} for i in range(3)]
        t=m.build(rows,pairs)
        self.assertIn('file0.fastq',t);self.assertIn('group2',t)
        self.assertIn('requires human review',t)
if __name__=='__main__':unittest.main()
