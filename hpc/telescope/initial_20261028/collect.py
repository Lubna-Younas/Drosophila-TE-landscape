"""Compile one species on a compute node after every expected sample is complete."""
import csv,fcntl,json,os,shutil,sys,tempfile
from pathlib import Path
from common import ROOT,PLAN,canonical,sha,verify_file,write_json,verify_package,compute_only
def collect(species):
 compute_only();rows=json.loads((ROOT/'manifests'/f'{species}.json').read_text())
 base=ROOT/'results'/species
 if not all((base/r['sample_id']/'receipt.json').exists() for r in rows):return False
 dest=ROOT/'matrices'/species;dest.parent.mkdir(exist_ok=True)
 with (dest.parent/(species+'.lock')).open('a') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX)
  if dest.exists():return True
  counts={};cpms={};keys=None;receipts={}
  for r in rows:
   d=base/r['sample_id'];rp=d/'receipt.json';receipt=json.loads(rp.read_text())
   if receipt['state']!='TE_COUNTS_COMPLETE' or receipt['package_sha256']!=verify_package():raise ValueError('Incomplete or mismatched sample')
   for x in receipt['files']:verify_file(dict(x,path=str(d/x['name'])))
   receipts[r['sample_id']]=sha(rp)
   with (d/'family_counts.tsv').open() as f:values=list(csv.DictReader(f,delimiter='\t'))
   data={(v['family'],v['class']):v for v in values}
   if keys is None:keys=sorted(data)
   if set(keys)!=set(data):raise ValueError('Mismatched species annotation families')
   counts[r['sample_id']]=[data[k]['estimated_fragments'] for k in keys]
   cpms[r['sample_id']]=[data[k]['CPM_total_input_fragments'] for k in keys]
  stage=Path(tempfile.mkdtemp(prefix='.'+species+'_',dir=dest.parent))
  for name,data in [('family_estimated_counts.tsv',counts),('family_CPM_total_input.tsv',cpms)]:
   with (stage/name).open('w') as f:
    w=csv.writer(f,delimiter='\t');w.writerow(['family','class']+[r['sample_id'] for r in rows])
    for i,key in enumerate(keys):w.writerow(list(key)+[data[r['sample_id']][i] for r in rows])
  shutil.copy2(ROOT/'manifests'/f'{species}.tsv',stage/'samples.tsv')
  write_json(stage/'receipt.json',{'state':'ALL_EXPECTED_TE_COUNTS_COLLECTED','species':species,'samples':len(rows),
      'scientific_QC':'NOT_YET_REVIEWED','package_sha256':verify_package(),'sample_receipt_hashes':receipts,
      'files':{p.name:sha(p) for p in stage.iterdir()}})
  stage.rename(dest);return True
if __name__=='__main__':
 if not collect(sys.argv[1]):raise SystemExit('Some expected samples are unfinished; no partial matrix published.')
