#!/usr/bin/env python3
"""Submit one Dpse prerequisite and 53 independent afterok jobs, with receipts."""
import argparse,fcntl,json,os,shlex,subprocess
from datetime import datetime,timezone
from pathlib import Path
from common import ROOT,PLAN,canonical,metadata,verify_package,write_json,sha

def scripts():
 rows=json.loads((ROOT/'manifests/Dpse.json').read_text())
 target=Path(PLAN['root']);result={}
 tasks=[('prepare',32,['prepare','Dpse'])]+[(r['sample_id'],r['memory_gb'],['sample','Dpse',r['sample_id']]) for r in rows]
 for key,mem,args in tasks:
  modules=['SAMtools/1.21-GCC-13.3.0']
  if key=='prepare':modules.append(PLAN['repeatmasker_module'])
  result[key]='\n'.join(['#!/usr/bin/env bash',
   '#SBATCH --job-name=TE3_Dpse_'+key[:60],'#SBATCH --ntasks=1','#SBATCH --cpus-per-task=4',
   '#SBATCH --mem='+str(mem)+'G','#SBATCH --time='+PLAN['walltime'],
   '#SBATCH --output='+str(target/'logs/Dpse'/(key+'_%j.out')),
   '#SBATCH --error='+str(target/'logs/Dpse'/(key+'_%j.err')),
   'set -euo pipefail','module purge','module load '+' '.join(modules),
   'unset PYTHONHOME' if key=='prepare' else 'unset PYTHONPATH PYTHONHOME',
   'export PYTHONNOUSERSITE=1',
   'export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1',
   'cd '+shlex.quote(str(target)),
   'exec python3 '+shlex.quote(str(target/'pipeline.py'))+' '+shlex.join(args),''])
 return rows,result

def validate_metadata():
 digest=verify_package();rows,expected=scripts()
 if len(rows)!=53 or PLAN['species']['Dpse']['sample_count']!=53:raise ValueError('Cohort mismatch')
 if len(set(r['sample_id'] for r in rows))!=53:raise ValueError('Duplicate sample IDs')
 for key,text in expected.items():
  if (ROOT/'jobs/Dpse'/(key+'.sbatch')).read_text()!=text:raise ValueError('Script mismatch: '+key)
 for row in rows:
  if row['species']!='Dpse' or row['run_accession'] in PLAN['excluded_runs']:raise ValueError('Unexpected sample')
  for kind in ('bam','bai'):
   current=metadata(row[kind]['path'])
   if any(current[k]!=row[kind][k] for k in ('bytes','mtime_ns')):raise ValueError('Input metadata changed: '+row[kind]['path'])
  receipt=json.loads(Path(row['mapping_receipt']).read_text())
  if canonical(receipt)!=row['mapping_receipt_sha256']:raise ValueError('Mapping receipt changed')
 source=PLAN['species']['Dpse']['annotation_source']
 if sha(source['path'])!=source['sha256']:raise ValueError('Source library changed')
 return rows,digest

def submit_one(key,digest,dependency=None):
 folder=ROOT/'submissions/Dpse';record=folder/(key+'.json');script=ROOT/'jobs/Dpse'/(key+'.sbatch')
 if record.exists():
  prior=json.loads(record.read_text())
  if prior['state']=='SUBMITTED' and prior['package_sha256']==digest and prior['dependency_job_id']==dependency:
   return prior
  raise RuntimeError('Prior submission requires inspection before retry: '+str(record))
 command=['sbatch','--parsable','--comment=Dpse_unknown_labels_v1']
 if dependency:command+=['--dependency=afterok:'+dependency,'--kill-on-invalid-dep=yes']
 command.append(str(script))
 value={'state':'SUBMITTING','species':'Dpse','sample':key,'script':str(script),
        'package_sha256':digest,'dependency_job_id':dependency,'command':command,
        'created_utc':datetime.now(timezone.utc).isoformat()}
 write_json(record,value)
 try:r=subprocess.run(command,text=True,capture_output=True,timeout=30)
 except Exception as exc:
  write_json(record,dict(value,state='SUBMISSION_UNCERTAIN',error=str(exc)));raise
 job=r.stdout.strip().split(';')[0]
 if r.returncode or not job.isdigit():
  write_json(record,dict(value,state='SUBMISSION_UNCERTAIN',stdout=r.stdout,stderr=r.stderr))
  raise RuntimeError('Submission not confirmed: '+str(record))
 value=dict(value,state='SUBMITTED',job_id=job,stdout=r.stdout,stderr=r.stderr)
 write_json(record,value);print(key,job,flush=True)
 return value

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--generate',action='store_true');ap.add_argument('--submit',action='store_true');a=ap.parse_args()
 if a.generate:
  if a.submit:raise ValueError('Cannot combine generate and submit')
  rows,generated=scripts();folder=ROOT/'jobs/Dpse';folder.mkdir(parents=True,exist_ok=True)
  for key,value in generated.items():
   with (folder/(key+'.sbatch')).open('x') as f:f.write(value)
  print('Generated',len(generated),'scripts; no jobs submitted');return
 if ROOT!=Path(PLAN['root']):raise ValueError('Run metadata checks and submission only in deployed target directory')
 rows,digest=validate_metadata()
 if not a.submit:
  print(json.dumps({'state':'METADATA_VALIDATED_NO_SUBMISSION','sample_count':len(rows),
   'jobs':54,'package_sha256':digest,'bam_bytes':sum(r['bam']['bytes'] for r in rows)},indent=2));return
 for d in ('logs','submissions'):(ROOT/d/'Dpse').mkdir(parents=True,exist_ok=True)
 with (ROOT/'submissions/Dpse/launcher.lock').open('a') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX)
  queue=subprocess.check_output(['squeue','-u',os.environ['USER'],'-h','-o','%i|%j|%T'],text=True)
  if any('|TE2_Dpse_' in x or '|TE_Dpse_' in x for x in queue.splitlines()):
   raise RuntimeError('An older Dpse job is active; inspect before creating duplicates')
  prep=submit_one('prepare',digest)
  children=[submit_one(row['sample_id'],digest,prep['job_id']) for row in rows]
  receipt={'state':'ALL_54_JOBS_SUBMITTED','created_utc':datetime.now(timezone.utc).isoformat(),
    'controller':str(ROOT),'package_sha256':digest,'preparation_job_id':prep['job_id'],
    'sample_jobs':[{'sample_id':r['sample'],'job_id':r['job_id'],'dependency_job_id':r['dependency_job_id']} for r in children],
    'submission_semantics':'53 independent jobs; afterok on preparation; invalid dependency cancels rather than runs'}
  write_json(ROOT/'submission_receipt.json',receipt)
  print(json.dumps({'state':receipt['state'],'preparation':prep['job_id'],'sample_count':len(children)},indent=2))
if __name__=='__main__':main()
