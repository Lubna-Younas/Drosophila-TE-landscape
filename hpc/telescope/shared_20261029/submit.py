#!/usr/bin/env python3
"""Metadata-only launcher. No arrays, dependencies, throttling, or fixed partitions."""
import argparse,fcntl,json,os,shlex,subprocess
from datetime import datetime,timezone
from pathlib import Path
from common import ROOT,PLAN,canonical,metadata,verify_package,write_json

def scripts(species):
 rows=json.loads((ROOT/'manifests'/f'{species}.json').read_text())
 commands=[('prepare',32,['prepare',species])]+[(r['sample_id'],r['memory_gb'],['sample',species,r['sample_id']]) for r in rows]
 result={}
 for key,mem,args in commands:
  modules=['SAMtools/1.21-GCC-13.3.0']
  if key=='prepare':modules+=['RepeatMasker/4.2.6-foss-2024a']
  text='\n'.join(['#!/usr/bin/env bash',
   '#SBATCH --job-name=TEH_'+species+'_'+key[:60], '#SBATCH --ntasks=1', '#SBATCH --cpus-per-task=4',
   '#SBATCH --mem='+str(mem)+'G','#SBATCH --time='+PLAN['walltime'],
   '#SBATCH --output='+str(ROOT/'logs'/species/(key+'_%j.out')),
   '#SBATCH --error='+str(ROOT/'logs'/species/(key+'_%j.err')),
   'set -euo pipefail','module purge','module load '+' '.join(modules),
   ('unset PYTHONHOME' if key=='prepare' else 'unset PYTHONPATH PYTHONHOME'),'export PYTHONNOUSERSITE=1',
   'export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1',
   'cd '+shlex.quote(str(ROOT)),
   'exec python3 '+shlex.quote(str(ROOT/'pipeline.py'))+' '+shlex.join(args),''])
  result[key]=text
 return rows,result

def eligible_samples(species,rows):
 # New annotations must never reuse counts produced against previous GTFs.
 eligible=[]
 for row in rows:
  final=ROOT/'results'/species/row['sample_id']/'receipt.json'
  if final.exists():
   d=json.loads(final.read_text())
   if d['state']!='TE_COUNTS_COMPLETE' or d['package_sha256']!=verify_package():raise ValueError('Existing result identity mismatch')
   continue
  eligible.append(row['sample_id'])
 return eligible

def main():
 ap=argparse.ArgumentParser();ap.add_argument('phase',choices=['prepare','samples']);ap.add_argument('species',choices=['Dmel','Dmir','Dpse']);ap.add_argument('--submit',action='store_true');a=ap.parse_args()
 digest=verify_package()
 if a.species not in PLAN['execution_species']:raise ValueError('Species outside approved controller scope')
 rows,generated=scripts(a.species)
 if len(rows)!=PLAN['species'][a.species]['sample_count']:raise ValueError('Cohort mismatch')
 for row in rows:
  if row['run_accession'] in PLAN['excluded_runs']:raise ValueError('Excluded sample in manifest')
  for kind in ('bam','bai'):
   current=metadata(row[kind]['path'])
   if any(current[k]!=row[kind][k] for k in ('bytes','mtime_ns')):
    raise ValueError('Input metadata changed; compute-node checksum revalidation is required: '+row[kind]['path'])
  receipt=json.loads(Path(row['mapping_receipt']).read_text())
  if canonical(receipt)!=row['mapping_receipt_sha256']:raise ValueError('Input receipt changed')
 if a.phase=='samples':
  from pipeline import validate_annotation
  validate_annotation(a.species)
 names=['prepare'] if a.phase=='prepare' else eligible_samples(a.species,rows)
 if not a.submit:
  print(json.dumps({'species':a.species,'samples':len(rows),'phase':a.phase,'planned_jobs':len(names),
   'mode':'DRY_RUN_NO_SUBMISSION','sample_jobs':'independent; submitted only after annotation and runtime checks',
   'walltime':PLAN['walltime'],'package_sha256':digest},indent=2));return
 for directory in ('logs','jobs','submissions'):(ROOT/directory/a.species).mkdir(parents=True,exist_ok=True)
 with (ROOT/'submissions'/a.species/'launcher.lock').open('a') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX)
  for name in names:
   record=ROOT/'submissions'/a.species/(name+'.json')
   if record.exists():
    old=json.loads(record.read_text())
    if old['state']=='SUBMITTED':print('Already submitted',name,old['job_id']);continue
    raise RuntimeError('Previous submission is uncertain/failed. Inspect '+str(record)+' before retrying.')
   script=ROOT/'jobs'/a.species/(name+'.sbatch')
   with script.open('x') as f:f.write(generated[name])
   value={'state':'SUBMITTING','species':a.species,'sample':name,'script':str(script),
          'package_sha256':digest,'created_utc':datetime.now(timezone.utc).isoformat()}
   write_json(record,value)
   try:r=subprocess.run(['sbatch','--parsable',str(script)],text=True,capture_output=True,timeout=30)
   except Exception as exc:
    write_json(record,dict(value,state='SUBMISSION_UNCERTAIN',error=str(exc)));raise
   job=r.stdout.strip().split(';')[0]
   if r.returncode or not job.isdigit():
    write_json(record,dict(value,state='SUBMISSION_UNCERTAIN',stdout=r.stdout,stderr=r.stderr));raise RuntimeError('Check '+str(record))
   write_json(record,dict(value,state='SUBMITTED',job_id=job,stdout=r.stdout,stderr=r.stderr));print(name,job)
if __name__=='__main__':main()
