#!/usr/bin/env python3
"""Verify selected assembly inputs and duplicate aliases on a compute node."""
import hashlib,json,os,subprocess
from pathlib import Path
from datetime import datetime,timezone
root=Path(__file__).resolve().parent
if not os.environ.get('SLURM_JOB_ID'):raise SystemExit('Use a compute job')
c=json.loads((root/'config.template.json').read_text())
if root!=Path(c['run_dir']):raise SystemExit('Unexpected staged directory')
if (root/'config.json').exists():raise SystemExit('A resolved configuration already exists; inspect previous attempt')
(root/'reference_preparation.lock').mkdir()
def digests(path):
 md5,sha=hashlib.md5(),hashlib.sha256()
 with Path(path).open('rb') as f:
  for block in iter(lambda:f.read(8*1024**2),b''):md5.update(block);sha.update(block)
 return {'md5':md5.hexdigest(),'sha256':sha.hexdigest()}
checks=[]
for key,source in c['reference_sources'].items():
 path=Path(source['path']);path.parent.mkdir(parents=True,exist_ok=True)
 if source['download']:
  partial=path.with_suffix(path.suffix+'.partial')
  if path.exists() or partial.exists():raise SystemExit('Prior reference source exists; review before reuse')
  with (root/'logs'/('download_'+key+'.log')).open('wb') as log:
   subprocess.run(['curl','--fail','--location','--retry','3','--max-time','900','--output',str(partial),source['url']],stdout=log,stderr=log,check=True,timeout=1200)
  hashes=digests(partial)
  if hashes['md5']!=source['md5']:raise ValueError('NCBI source checksum mismatch: '+key)
  partial.rename(path)
 else:
  hashes=digests(path)
  if hashes['md5']!=source['md5'] or hashes['sha256']!=source['sha256']:raise ValueError('Existing reference checksum mismatch: '+key)
 c[key+'_sha256']=hashes['sha256'];checks.append({'key':key,'path':str(path),'bytes':path.stat().st_size,**hashes})
alias_checks=[]
for item in c['duplicate_aliases']:
 one,two=Path(item['canonical']),Path(item['alias'])
 if one.stat().st_size!=two.stat().st_size:raise ValueError('Alias size mismatch; do not deduplicate')
 a,b=digests(one),digests(two)
 if a!=b:raise ValueError('Alias contents differ; sample grouping needs review')
 alias_checks.append({**item,**a,'bytes':one.stat().st_size,'state':'BYTE_IDENTICAL_ALIAS'})
(root/'reference_source_verification.json').write_text(json.dumps({'time_utc':datetime.now(timezone.utc).isoformat(),'job_id':os.environ['SLURM_JOB_ID'],'sources':checks,'aliases':alias_checks},indent=2)+'\n')
with (root/'config.json').open('x') as f:json.dump(c,f,indent=2);f.write('\n')
print(json.dumps({'reference_sources_verified':c['assembly'],'duplicate_alias_files_verified':len(alias_checks)}),flush=True)
