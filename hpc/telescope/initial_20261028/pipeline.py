#!/usr/bin/env python3
"""Compute-node annotation preparation and independent sample TE quantification."""
import argparse,collections,csv,fcntl,json,math,os,re,shutil,socket,subprocess,tempfile,traceback
from datetime import datetime,timezone
from pathlib import Path
from common import ROOT,PLAN,sha,canonical,metadata,write_json,verify_file,verify_package,compute_only,run,dictionary_from_bam,isolated_env
import annotation
from runtime_patch import apply_patch as patch_runtime,verify_patch

def utc():return datetime.now(timezone.utc).isoformat()
def runtime(scratch):
 env=ROOT/'runtime/telescope-1.0.4.1'
 (ROOT/'runtime').mkdir(exist_ok=True)
 for key in ('XDG_CACHE_HOME','CONDA_PKGS_DIRS','CONDA_ENVS_PATH','PIP_CACHE_DIR'):
  Path(isolated_env()[key]).mkdir(parents=True,exist_ok=True)
 with (ROOT/'runtime/setup.lock').open('a') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX)
  ready=ROOT/'runtime/READY.json'
  if ready.exists():
   data=json.loads(ready.read_text())
   if data['package_sha256']!=verify_package():raise ValueError('Runtime belongs to different package')
   verify_patch()
   return env
  log=ROOT/'runtime/setup.log'
  if not env.exists():
   run([PLAN['conda'],'create','--yes','--solver=libmamba','--override-channels','--strict-channel-priority',
        '-c','conda-forge','-c','bioconda','--prefix',env,'python=3.6',
        'telescope='+PLAN['telescope_version']+'='+PLAN['telescope_build'],'pandas=1.1.3'],log)
  patch_receipt=patch_runtime(env)
  telescope=env/'bin/telescope'
  version=subprocess.check_output([str(telescope),'--version'],text=True,env=isolated_env()).strip()
  if '1.0.4.1' not in version:raise ValueError('Unexpected Telescope version: '+version)
  run([env/'bin/python',ROOT/'runtime_smoke.py',scratch/'runtime_smoke'],log)
  shutil.copy2(scratch/'runtime_smoke/checks.json',ROOT/'runtime/synthetic_validation.json')
  spec=subprocess.check_output([PLAN['conda'],'list','--explicit','--prefix',str(env)],text=True,env=isolated_env())
  (ROOT/'runtime/explicit_environment.txt').write_text(spec)
  write_json(ready,{'state':'INSTALLED_AND_SYNTHETICALLY_VALIDATED','version':version,
      'package_sha256':verify_package(),'environment_spec_sha256':sha(ROOT/'runtime/explicit_environment.txt'),
      'created_utc':utc(),'job_id':os.environ['SLURM_JOB_ID'],'runtime_patch':patch_receipt})
 return env

def validate_annotation(species):
 folder=ROOT/'annotations'/species;r=json.loads((folder/'receipt.json').read_text())
 if r['state']!='ANNOTATION_VALIDATED' or r['assembly']!=PLAN['species'][species]['assembly']:
  raise ValueError('Annotation is not validated for this assembly')
 if r['package_sha256']!=verify_package():raise ValueError('Annotation package changed')
 for x in r['files']:verify_file(dict(x,path=str(folder/x['name'])))
 return folder,r

def prepare(species,scratch):
 cfg=PLAN['species'][species];dest=ROOT/'annotations'/species
 runtime(scratch)
 if dest.exists():validate_annotation(species);return
 log=scratch/'prepare.log';work=scratch/'annotation';work.mkdir()
 rr=json.loads(Path(cfg['reference_receipt']).read_text())
 if canonical(rr)!=cfg['reference_receipt_sha256']:raise ValueError('Mapping reference receipt changed')
 genome_meta=verify_file(cfg['genome'],always_hash=True)
 genome=annotation.read_genome(cfg['genome']['path']);lengths={k:len(v) for k,v in genome.items()}
 samples=json.loads((ROOT/'manifests'/f'{species}.json').read_text())
 for row in samples:
  verify_file(row['bam']);verify_file(row['bai'])
  if dictionary_from_bam(row['bam']['path'])!=lengths:raise ValueError('BAM/genome dictionary mismatch: '+row['sample_id'])
 src=cfg['annotation_source'];source=Path(src['path'])
 if sha(source)!=src['sha256']:raise ValueError('Annotation source checksum changed')
 if species=='Dmel':models,stats=annotation.flybase(source,genome)
 elif species=='Dmir':models,stats=annotation.repeatmasker(source,lengths,species)
 else:
  library=work/'library.fa';catalog=annotation.normalize_library(source,library)
  query=work/'genome.fa';shutil.copyfile(cfg['genome']['path'],query)
  # RepeatMasker -pa 1 can use four RMBlast threads, fitting the four-CPU allocation.
  run(['RepeatMasker','-engine','rmblast','-pa','1','-nolow','-no_is','-lib',library,'-dir',work,query],log)
  models,stats=annotation.repeatmasker(work/'genome.fa.out',lengths,species,catalog)
  shutil.copyfile(work/'genome.fa.out',work/'Dpse_exact_genome_repeatmasker.out')
 annotation.export(models,work)
 if metadata(cfg['genome']['path'])!={k:genome_meta[k] for k in ('path','bytes','mtime_ns')}:
  raise ValueError('Genome changed during annotation preparation')
 write_json(work/'genome_dictionary.json',lengths)
 write_json(work/'annotation_validation.json',stats)
 keep=['TE.gtf','locus_catalog.tsv','genome_dictionary.json','annotation_validation.json']
 if species=='Dpse':keep+=['Dpse_exact_genome_repeatmasker.out','library.fa']
 parent=dest.parent;parent.mkdir(exist_ok=True)
 staging=Path(tempfile.mkdtemp(prefix='.'+species+'_',dir=parent))
 for name in keep:shutil.copy2(work/name,staging/name)
 if log.exists():shutil.copy2(log,staging/'prepare.log');keep+=['prepare.log']
 receipt={'state':'ANNOTATION_VALIDATED','scientific_QC':'ANNOTATION_REVIEW_PENDING',
    'species':species,'assembly':cfg['assembly'],'created_utc':utc(),'genome':genome_meta,
    'annotation_source':src,'validation':stats,'package_sha256':verify_package(),
    'job_id':os.environ['SLURM_JOB_ID'],'files':[dict(name=n,**{k:v for k,v in metadata(staging/n).items() if k!='path'},sha256=sha(staging/n)) for n in keep]}
 write_json(staging/'receipt.json',receipt)
 if dest.exists():raise FileExistsError(str(dest))
 staging.rename(dest)

def summarize(report,catalog,output,row):
 with open(catalog) as f:loci={r['locus']:r for r in csv.DictReader(f,delimiter='\t')}
 counts={};totals=collections.defaultdict(lambda:[0.0,0.0])
 if report.exists():
  with report.open() as f:
   for r in csv.DictReader((x for x in f if not x.startswith('#')),delimiter='\t'):
    key=r['transcript']
    if key=='__no_feature':continue
    if key not in loci or key in counts:raise ValueError('Unexpected or duplicate Telescope locus')
    vals=[float(r['final_count']),float(r['unique_count'])]
    if any(not math.isfinite(v) or v<0 for v in vals):raise ValueError('Invalid TE count')
    counts[key]=vals
 for key,m in loci.items():
  vals=counts.get(key,[0.,0.]);totals[(m['family'],m['class'])][0]+=vals[0];totals[(m['family'],m['class'])][1]+=vals[1]
 with (output/'locus_counts.tsv').open('w') as f:
  w=csv.writer(f,delimiter='\t');w.writerow(['locus','family','class','estimated_fragments','model_unique_fragments'])
  for key,m in sorted(loci.items()):w.writerow([key,m['family'],m['class']]+counts.get(key,[0.,0.]))
 with (output/'family_counts.tsv').open('w') as f:
  w=csv.writer(f,delimiter='\t');w.writerow(['family','class','estimated_fragments','model_unique_fragments','CPM_total_input_fragments'])
  for (family,cl),vals in sorted(totals.items()):w.writerow([family,cl]+vals+[vals[0]*1e6/row['total_input_fragments']])
 return {'loci':len(loci),'family_class_groups':len(totals),'estimated_TE_fragments':sum(v[0] for v in totals.values()),
         'model_unique_TE_fragments':sum(v[1] for v in totals.values())}

def quantify(species,sid,scratch):
 records=json.loads((ROOT/'manifests'/f'{species}.json').read_text())
 row=next(r for r in records if r['sample_id']==sid)
 if row['run_accession'] in PLAN['excluded_runs']:raise ValueError('Excluded sample')
 dest=ROOT/'results'/species/sid
 if dest.exists():raise FileExistsError('Existing results preserved: '+str(dest))
 cfg=PLAN['species'][species];ann,annreceipt=validate_annotation(species)
 receipt=json.loads(Path(row['mapping_receipt']).read_text())
 if canonical(receipt)!=row['mapping_receipt_sha256']:raise ValueError('Mapping receipt changed')
 bam_before=verify_file(row['bam']);verify_file(row['bai'])
 if dictionary_from_bam(row['bam']['path'])!=json.loads((ann/'genome_dictionary.json').read_text()):raise ValueError('BAM annotation mismatch')
 log=scratch/'quantification.log';run(['samtools','quickcheck','-v',row['bam']['path']],log)
 # Collation preserves secondary alignments. Do not use samtools collate -f.
 required=3*row['bam']['bytes']+10*1024**3
 if shutil.disk_usage(scratch).free<required:raise RuntimeError('Insufficient temporary space: need '+str(required)+' bytes')
 collated=scratch/'all_alignments.collated.bam'
 run(['samtools','collate','-@','3','-o',collated,row['bam']['path'],scratch/'collate_tmp'],log)
 run(['samtools','quickcheck','-v',collated],log)
 env=ROOT/'runtime/telescope-1.0.4.1'
 ready=json.loads((ROOT/'runtime/READY.json').read_text())
 if ready['package_sha256']!=verify_package():raise ValueError('Runtime validation is missing or stale')
 verify_patch()
 output=scratch/'counts';output.mkdir()
 command=[env/'bin/telescope','assign',collated,ann/'TE.gtf','--attribute','locus','--ncpu','1',
          '--stranded_mode',row['telescope_stranded_mode'],'--reassign_mode','average',
          '--overlap_mode','threshold','--overlap_threshold','0.2','--max_iter','1000',
          '--outdir',output,'--exp_tag',sid,'--tempdir',scratch]
 run(command,log)
 report=output/(sid+'-telescope_report.tsv');text=log.read_text()
 no_overlap='No alignments overlapping annotation' in text
 if not no_overlap and (not report.exists() or 'EM converged after ' not in text):
  raise RuntimeError('Telescope did not produce a converged result; inspect the preserved log.')
 if 'telescope assign complete' not in text:raise RuntimeError('Missing completion message')
 summary=summarize(report,ann/'locus_catalog.tsv',output,row)
 if metadata(row['bam']['path'])!={k:bam_before[k] for k in ('path','bytes','mtime_ns')}:raise ValueError('Input BAM changed during counting')
 parent=dest.parent;parent.mkdir(parents=True,exist_ok=True)
 stage=Path(tempfile.mkdtemp(prefix='.'+sid+'_',dir=parent))
 # Only final counts, model report, logs and provenance leave temporary space.
 for p in output.glob('*.tsv'):shutil.copy2(p,stage/p.name)
 shutil.copy2(log,stage/'quantification.log')
 result={'state':'TE_COUNTS_COMPLETE','scientific_QC':'NOT_YET_REVIEWED','sample':row,
         'annotation_receipt_sha256':canonical(annreceipt),'bam_verified':bam_before,'summary':summary,
         'zero_overlap':no_overlap,'created_utc':utc(),'host':socket.gethostname(),'job_id':os.environ['SLURM_JOB_ID'],
         'telescope_version':ready['version'],'runtime_patch':ready['runtime_patch'],'package_sha256':verify_package(),'command':list(map(str,command)),
         'count_rounding':'Telescope average-mode report stores two decimal places; family sums use those values.',
         'expression_scope':'TE-overlapping RNA signal; autonomous TE transcription is not established.'}
 result['files']=[dict(name=p.name,**{k:v for k,v in metadata(p).items() if k!='path'},sha256=sha(p)) for p in stage.iterdir()]
 write_json(stage/'receipt.json',result)
 if dest.exists():raise FileExistsError(str(dest))
 stage.rename(dest)

def main():
 ap=argparse.ArgumentParser();ap.add_argument('action',choices=['runtime','prepare','sample']);ap.add_argument('species',choices=['Dmel','Dmir','Dpse']);ap.add_argument('sample',nargs='?');a=ap.parse_args()
 compute_only();verify_package()
 for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
 key=a.species+'_'+(a.sample or a.action);(ROOT/'locks').mkdir(exist_ok=True)
 with (ROOT/'locks'/(key+'.lock')).open('a') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
  scratch=Path(tempfile.mkdtemp(prefix='TE_'+key+'_',dir=os.environ.get('SLURM_TMPDIR') or '/tmp'))
  try:
   if a.action=='runtime':
    runtime(scratch)
   elif a.action=='prepare':
    prepare(a.species,scratch)
    # Sample jobs enter Slurm only after validated prerequisites exist.
    subprocess.run(['python3',str(ROOT/'submit.py'),'samples',a.species,'--submit'],check=True)
   else:
    quantify(a.species,a.sample,scratch)
    from collect import collect
    collect(a.species)
  except BaseException as exc:
   failure=ROOT/'failures'/(key+'_'+os.environ['SLURM_JOB_ID']);failure.mkdir(parents=True,exist_ok=True)
   for p in scratch.rglob('*.log'):shutil.copy2(p,failure/(p.parent.name+'_'+p.name))
   for p in scratch.rglob('*telescope_report.tsv'):shutil.copy2(p,failure/p.name)
   write_json(failure/'failure.json',{'state':'FAILED','error':str(exc),'traceback':traceback.format_exc(),'job_id':os.environ['SLURM_JOB_ID'],'time':utc()})
   raise
  finally:shutil.rmtree(scratch)
if __name__=='__main__':main()
