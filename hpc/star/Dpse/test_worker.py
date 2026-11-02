import tempfile,gzip,hashlib,unittest,importlib.util
from pathlib import Path
s=importlib.util.spec_from_file_location('worker',Path(__file__).with_name('mapping_worker.py'));w=importlib.util.module_from_spec(s);s.loader.exec_module(w)
RAW=b'@READ1\nACGT\n+\nIIII\n@READ2\nTT\n+\nII\n'
class WorkerTests(unittest.TestCase):
 def setUp(self):self.t=tempfile.TemporaryDirectory();self.r=Path(self.t.name)
 def tearDown(self):self.t.cleanup()
 def make(self,name,depth,raw=RAW):
  data=raw
  for _ in range(depth):data=gzip.compress(data,mtime=0)
  p=self.r/name;p.write_bytes(data);return p
 def test_nested_wrapper_removal_preserves_reads(self):
  a=self.make('a.gz',2);before=a.read_bytes();paths,e=w.normalize_fastq_compression([a],[2],self.r)
  self.assertEqual(gzip.decompress(Path(paths[0]).read_bytes()),RAW);self.assertEqual(a.read_bytes(),before)
  self.assertEqual(e[0]['derived_gzip_sha256'],hashlib.sha256(gzip.decompress(before)).hexdigest())
  self.assertEqual(w.check_fastq(paths,'SE')['checked_records_per_file'],2)
 def test_one_layer_unchanged(self):
  a=self.make('a.gz',1);paths,e=w.normalize_fastq_compression([a],[1],self.r);self.assertEqual(paths,[a]);self.assertEqual(e[0]['transformation'],'none')
 def test_mixed_compression_pair(self):
  a=self.make('a.gz',2);b=self.make('b.gz',1);paths,_=w.normalize_fastq_compression([a,b],[2,1],self.r);self.assertEqual(w.check_fastq(paths,'PE')['checked_records_per_file'],2)
 def test_unexpected_depth(self):
  a=self.make('a.gz',2)
  with self.assertRaises(ValueError):w.normalize_fastq_compression([a],[1],self.r)
 def test_triple_gzip_rejected(self):
  a=self.make('a.gz',3)
  with self.assertRaises(ValueError):w.normalize_fastq_compression([a],[2],self.r)
 def test_truncated_outer_rejected(self):
  a=self.make('a.gz',2);a.write_bytes(a.read_bytes()[:-5])
  with self.assertRaises((EOFError,OSError)):w.normalize_fastq_compression([a],[2],self.r)
 def test_mate_mismatch_rejected(self):
  a=self.make('a.gz',1);b=self.make('b.gz',1,RAW.replace(b'READ2',b'READ3'))
  with self.assertRaises(ValueError):w.check_fastq([a,b],'PE')
 def test_gtf_no_exclusions(self):
  p=self.r/'original.gtf';p.write_text('chr1\tNCBI\texon\t1\t10\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n')
  result=w.prepare_rsem_gtf(p,self.r/'genes.gtf',[]);self.assertEqual(result['retained_transcripts'],1);self.assertEqual(result['excluded_transcripts'],[])
 def test_gtf_mixed_strand_requires_review(self):
  p=self.r/'original.gtf';p.write_text('chr1\tNCBI\texon\t1\t10\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\nchr1\tNCBI\texon\t20\t30\t.\t-\t.\tgene_id "g1"; transcript_id "t1";\n')
  with self.assertRaises(ValueError):w.prepare_rsem_gtf(p,self.r/'genes.gtf',[])
if __name__=='__main__':unittest.main(verbosity=2)
