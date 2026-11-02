"""Audited two-line SE strand correction in the isolated Telescope 1.0.4.1 runtime."""
import os
from pathlib import Path
from common import ROOT,PLAN,sha,write_json
ORIGINAL='5bf65c71c4e84d8e95371819aa869f0bfa1b9e78ff82e08125a35bf0d594847e'
PATCHED='efb4eb672cce427ce9699f005f2b28478eb22dacc09b396fbea45ce369741b8f'

def apply_patch(env):
 p=Path(env)/'lib/python3.6/site-packages/telescope/utils/model.py'
 folder=ROOT/'runtime/patches';folder.mkdir(exist_ok=True)
 original=folder/'model.py.original'
 if sha(p)==ORIGINAL:
  source=p.read_text()
  if not original.exists():
   with original.open('x') as f:f.write(source)
  assert sha(original)==ORIGINAL
  a="frag_strand = '+' if self.opts.stranded_mode[0] == 'F' else '-'"
  b="frag_strand = '-' if self.opts.stranded_mode[0] == 'F' else '+'"
  assert source.count(a)==1 and source.count(b)==1
  revised=source.replace(a,'SE_REVERSE_FIX').replace(b,a).replace('SE_REVERSE_FIX',b)
  tmp=p.with_suffix('.py.repair');tmp.write_text(revised)
  assert sha(tmp)==PATCHED
  # Replace the inode; never modify a Conda package-cache hard link in place.
  os.replace(tmp,p)
 assert sha(p)==PATCHED and sha(original)==ORIGINAL,'Unexpected Telescope source'
 receipt={'patch':'SE_STRAND_FIX_V1','original_sha256':ORIGINAL,'patched_sha256':PATCHED,
          'path':str(p),'scope':'Correct SE F/R orientation only; paired-end and unstranded logic unchanged'}
 write_json(folder/'receipt.json',receipt)
 return receipt

def verify_patch():
 p=Path(PLAN['telescope_env'])/'lib/python3.6/site-packages/telescope/utils/model.py'
 if sha(p)!=PATCHED:raise ValueError('Telescope SE strand patch missing or changed')
