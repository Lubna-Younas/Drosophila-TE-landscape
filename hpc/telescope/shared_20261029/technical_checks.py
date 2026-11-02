"""Small synthetic checks for shared-library grouping and genomic coverage."""
from pathlib import Path
import json,tempfile,unittest
import annotation
from landscape_annotation import make_models,revcomp,interval_union,landscape_tables

class AnnotationChecks(unittest.TestCase):
 def test_exact_sequence_aliases_keep_class_conflicts_separate(self):
  seq='ACGTTTAAAC'
  df=[('DF000000001.1#DNA/hAT','ACGTTTAAAC'),('DF000000002.1#rRNA','ACGT'),('IS1.None#Unknown','ACGT')]
  hb=[('same f hAT DNA',revcomp(seq)),('candidate f ? ?',seq)]
  models,excluded,n=make_models(df,hb,{'DF000000001.1':'named_model'})
  self.assertEqual(n,3);self.assertEqual(len(models),2);self.assertEqual(len(excluded),2)
  self.assertEqual(len(models['DF000000001.1']['aliases']),2)
  self.assertEqual(models['HB2018_candidate']['class'],'Unknown')
 def test_invalid_consensus_rejected(self):
  with self.assertRaises(ValueError):make_models([('DF000000001.1#DNA/hAT','ACGTX')],[],{})
 def test_interval_union_preserves_gaps(self):
  self.assertEqual(interval_union([(10,15),(1,3),(3,5),(7,8)]),[(1,5),(7,8),(10,15)])
 def test_nested_repeat_blocks_and_star_policy(self):
  text='100 0.0 0.0 0.0 chr 1 3 (7) + DF000000001.1 DNA/hAT 1 3 (0) 1\n'
  text+='100 0.0 0.0 0.0 chr 7 9 (1) + DF000000001.1 DNA/hAT 4 6 (0) 1\n'
  text+='80 0.0 0.0 0.0 chr 2 4 (6) + DF000000001.1 DNA/hAT 1 3 (0) 2 *\n'
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/'test.out';p.write_text(text)
   models,stats=annotation.repeatmasker(p,{'chr':10},'TEST')
  self.assertEqual(len(models),1);self.assertEqual(next(iter(models.values()))['blocks'],[(1,3),(7,9)])
  self.assertEqual(stats['source_rows'],3);self.assertEqual(stats['excluded_non_TE_rows']['LOWER_PRIORITY_OVERLAPPING_MATCH'],1)
 def test_non_N_coverage_uses_global_union(self):
  genome={'chr':'ACGTNNACGT'}
  models={'one':{'family':'a','class':'DNA','chrom':'chr','blocks':[(1,8)]},'two':{'family':'b','class':'LTR','chrom':'chr','blocks':[(7,10)]}}
  with tempfile.TemporaryDirectory() as d:
   landscape_tables(models,genome,Path(d));data=json.loads((Path(d)/'genomic_coverage_summary.json').read_text())
  self.assertEqual(data['assembly_non_N_bases'],8);self.assertEqual(data['covered_non_N_bases_union'],8);self.assertEqual(data['percent_assembly_non_N_bases'],100)

if __name__=='__main__':unittest.main()
