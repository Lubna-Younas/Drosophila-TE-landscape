"""Frozen shared reference-model library and assembly-specific TE coordinates."""
from pathlib import Path
import bisect,collections,csv,datetime,fcntl,gzip,json,os,re,shutil,subprocess,tempfile
import annotation
from common import ROOT,PLAN,sha,canonical,metadata,write_json,verify_package,verify_file,compute_only,run,dictionary_from_bam

NON_TE=set(PLAN['non_TE_classes'])
def utc():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def revcomp(seq):return seq.translate(str.maketrans('ACGTRYMKSWBDHVN','TGCAYRKMSWVHDBN'))[::-1]
def interval_union(intervals):
 result=[]
 for a,b in sorted(intervals):
  if result and a<=result[-1][1]+1:result[-1]=(result[-1][0],max(result[-1][1],b))
  else:result.append((a,b))
 return result
def normalize_class(value):
 return value.replace('?/?','Unknown').replace('Unknown/Unknown','Unknown')
def make_models(dfam_records,obscura_records,names):
 records=[];excluded=[]
 for header,seq in dfam_records:
  token=header.split()[0];identifier,sep,cl=token.partition('#')
  if not re.fullmatch(r'DF\d+\.\d+',identifier):
   excluded.append({'source':'Dfam','header':header,'reason':'NOT_CURATED_DF_ACCESSION'});continue
  if not sep:raise ValueError('Dfam class label missing: '+header)
  cl=normalize_class(cl)
  if cl.split('/')[0] in NON_TE:
   excluded.append({'source':'Dfam','header':header,'reason':'NON_TE_CLASS','class':cl});continue
  records.append({'id':identifier,'source':'Dfam4.0','source_id':identifier,'display_name':names.get(identifier,identifier),'class':cl,'sequence':seq})
 for header,seq in obscura_records:
  fields=header.split()
  if len(fields)<4:raise ValueError('Unexpected obscura header: '+header)
  name,family,subclass,cl=fields[:4];cl=normalize_class(cl+'/'+subclass)
  if cl.split('/')[0] in NON_TE:
   excluded.append({'source':'Hill_Betancourt_2018','header':header,'reason':'NON_TE_CLASS','class':cl});continue
  identifier='HB2018_'+re.sub(r'[^A-Za-z0-9_.-]','_',name)
  records.append({'id':identifier,'source':'Hill_Betancourt_2018','source_id':name,'display_name':family,'class':cl,'sequence':seq})
 if len({r['id'] for r in records})!=len(records):raise ValueError('Reference-model ID collision')
 groups=collections.defaultdict(list)
 for row in records:
  seq=row['sequence']
  if not seq or re.fullmatch(r'[ACGTRYMKSWBDHVN]+',seq) is None:raise ValueError('Invalid consensus sequence')
  groups[(min(seq,revcomp(seq)),row['class'])].append(row)
 result={}
 for (_,cl),members in groups.items():
  members.sort(key=lambda r:(r['source']!='Dfam4.0',r['id']))
  representative=members[0];key=representative['id']
  result[key]={'family':key,'class':cl,'sequence':representative['sequence'],
               'aliases':[{k:v for k,v in m.items() if k!='sequence'} for m in members],
               'interpretation':'CANDIDATE_UNKNOWN' if cl.split('/')[0] in ('Unknown','Novel') else 'CLASSIFIED_REFERENCE_MODEL'}
 return result,excluded,len(records)

def library(scratch):
 folder=ROOT/'shared_library';(ROOT/'locks').mkdir(exist_ok=True)
 with (ROOT/'locks/shared_library.lock').open('a') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX)
  if (folder/'receipt.json').exists():
   receipt=json.loads((folder/'receipt.json').read_text())
   if receipt['state']!='SHARED_LIBRARY_FROZEN' or receipt['package_sha256']!=verify_package():raise ValueError('Shared library identity mismatch')
   for f in receipt['files']:verify_file(dict(f,path=str(folder/f['name'])))
   return folder,json.loads((folder/'model_catalog.json').read_text()),receipt
  cfg=PLAN['shared_library'];work=scratch/'shared_library';work.mkdir()
  before=[]
  for spec in cfg['database_files']:
   m=metadata(spec['path'])
   if any(m[k]!=spec[k] for k in ['bytes','mtime_ns']):raise ValueError('Installed Dfam partition changed')
   before.append(dict(m,sha256=sha(spec['path'])))
  source=cfg['obscura_source']
  if sha(source['path'])!=source['sha256']:raise ValueError('Published obscura library changed')
  prefix=[cfg['famdb_executable'],'-i',cfg['database_directory'],'-e','ch,uc,uh']
  command=prefix+['families','-ad','-c','-f','fasta_acc','--include-class-in-name',str(cfg['taxon_id'])]
  with (work/'Dfam4_Drosophila_export.fa').open('w') as f:
   p=subprocess.run(command,stdout=f,stderr=subprocess.PIPE,text=True)
  (work/'famdb_export.stderr').write_text(p.stderr)
  if p.returncode:raise RuntimeError('FamDB export failed: '+p.stderr)
  summary=subprocess.check_output(prefix+['families','-ad','-c','-f','summary',str(cfg['taxon_id'])],text=True)
  (work/'Dfam4_Drosophila_model_summary.txt').write_text(summary)
  names={m.group(1):m.group(2) for m in re.finditer(r"^(\S+) '([^']+)':",summary,re.M)}
  dfam=list(annotation.fasta(work/'Dfam4_Drosophila_export.fa'))
  dfcount=sum(bool(re.fullmatch(r'DF\d+\.\d+',h.split()[0].split('#')[0])) for h,s in dfam)
  if dfcount!=cfg['expected_dfam_curated_DF_records']:raise ValueError('Unexpected Dfam curated export count: '+str(dfcount))
  obscura=list(annotation.fasta(source['path']))
  if len(obscura)!=303:raise ValueError('Unexpected published obscura record count')
  models,excluded,input_count=make_models(dfam,obscura,names)
  if not models:raise ValueError('Shared library empty')
  with (work/'library.fa').open('w') as f:
   for key,value in sorted(models.items()):f.write('>'+key+'#'+value['class']+'\n'+value['sequence']+'\n')
  smoke=work/'smoke';smoke.mkdir()
  exemplar=next(v['sequence'] for v in models.values() if len(v['sequence'])>=1500 and v['class'].split('/')[0]!='Unknown')
  (smoke/'query.fa').write_text('>known_consensus\n'+exemplar+'\n')
  run(['python3','-c','import h5py; print(h5py.__version__)'],smoke/'smoke.log')
  run(['RepeatMasker','-v'],smoke/'smoke.log')
  run(['rmblastn','-version'],smoke/'smoke.log')
  run(['RepeatMasker','-engine','rmblast','-pa','1','-no_is','-lib',work/'library.fa','-dir',smoke,smoke/'query.fa'],smoke/'smoke.log')
  smoke_models,_=annotation.repeatmasker(smoke/'query.fa.out',{'known_consensus':len(exemplar)},'SMOKE',{k:{'family':v['family'],'class':v['class']} for k,v in models.items()})
  if not smoke_models:raise ValueError('Shared-library RepeatMasker self-match failed')
  shutil.copy2(smoke/'smoke.log',work/'shared_library_smoke.log')
  shutil.copy2(smoke/'query.fa.out',work/'shared_library_smoke.out')
  write_json(work/'model_catalog.json',models)
  with (work/'model_catalog.tsv').open('w') as f:
   w=csv.writer(f,delimiter='\t');w.writerow(['model_group','class','interpretation','source_ids','source_names','consensus_bases','family_equivalence_status'])
   for key,v in sorted(models.items()):w.writerow([key,v['class'],v['interpretation'],';'.join(a['source_id'] for a in v['aliases']),';'.join(a['display_name'] for a in v['aliases']),len(v['sequence']),'REFERENCE_MODEL_GROUP_NOT_AUTOMATIC_BIOLOGICAL_FAMILY'])
  shutil.copy2(source['path'],work/'Hill_Betancourt_2018_original.fa')
  write_json(work/'excluded_library_records.json',excluded)
  for spec,record in zip(cfg['database_files'],before):
   if metadata(spec['path'])!={k:record[k] for k in ['path','bytes','mtime_ns']}:raise ValueError('Dfam source changed during export')
  files=[dict(name=p.name,bytes=p.stat().st_size,mtime_ns=p.stat().st_mtime_ns,sha256=sha(p)) for p in work.iterdir() if p.is_file()]
  receipt={'state':'SHARED_LIBRARY_FROZEN','created_utc':utc(),'package_sha256':verify_package(),'job_id':os.environ['SLURM_JOB_ID'],
   'dfam_version':cfg['dfam_version'],'dfam_source_partitions':before,'export_command':command,'dfam_export_records':len(dfam),'dfam_curated_DF_records':dfcount,
   'obscura_source':source,'obscura_records':len(obscura),'retained_models_before_exact_deduplication':input_count,
   'shared_model_groups':len(models),'exact_sequence_duplicate_models_collapsed':input_count-len(models),
   'repeatmasker_library_smoke':'PASSED_FULL_CONSENSUS_QUERY',
   'unknown_candidate_model_groups':sum(v['interpretation']=='CANDIDATE_UNKNOWN' for v in models.values()),
   'interpretation':cfg['biological_family_resolution'],'files':files}
  write_json(work/'receipt.json',receipt)
  if folder.exists():raise FileExistsError('Unfinished shared library requires inspection')
  shutil.copytree(work,folder)
  return folder,models,receipt

def landscape_tables(models,genome,work):
 total=sum(len(s) for s in genome.values());non_n=sum(len(s)-s.count('N') for s in genome.values())
 groups=collections.defaultdict(lambda:collections.defaultdict(list));classes=collections.defaultdict(lambda:collections.defaultdict(list));all_blocks=collections.defaultdict(list)
 loci=collections.Counter();class_loci=collections.Counter()
 for value in models.values():
  key=(value['family'],value['class']);loci[key]+=1;class_loci[value['class']]+=1
  for a,b in value['blocks']:
   groups[key][value['chrom']].append((a,b));classes[value['class']][value['chrom']].append((a,b));all_blocks[value['chrom']].append((a,b))
 def bases(bychrom):
  intervals={c:interval_union(v) for c,v in bychrom.items()}
  return sum(sum((b-a+1)-genome[c][a-1:b].count('N') for a,b in v) for c,v in intervals.items())
 def table(path,kind,values,counts):
  with path.open('w') as f:
   w=csv.writer(f,delimiter='\t');w.writerow([kind,'class','locus_models','covered_non_N_bases','percent_assembly_non_N_bases'])
   for key,blocks in sorted(values.items()):
    model,cl=key if isinstance(key,tuple) else ('ALL_'+key,key)
    bp=bases(blocks);w.writerow([model,cl,counts[key],bp,100*bp/non_n])
 table(work/'genomic_model_coverage.tsv','model_group',groups,loci)
 table(work/'genomic_class_coverage.tsv','class_group',classes,class_loci)
 overall=bases(all_blocks)
 write_json(work/'genomic_coverage_summary.json',{'assembly_bases':total,'assembly_non_N_bases':non_n,'covered_non_N_bases_union':overall,'percent_assembly_non_N_bases':100*overall/non_n,
  'annotation_locus_models':len(models),'unit':'known-library matching locus models; not validated intact insertion counts',
  'scope':'reference assembly; not the TE insertion repertoire of every RNA sample strain','class_model_coverage_note':'Union within each group; sums across groups can overlap and are not a global coverage denominator.'})

def flybase_comparison(models,work):
 prior=Path(PLAN['base_controller'])/'annotations/Dmel/locus_catalog.tsv'
 with prior.open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
 bychrom=collections.defaultdict(list)
 for m in models.values():
  for a,b in m['blocks']:bychrom[m['chrom']].append((a,b,m['locus'],m['family']))
 index={}
 for chrom,values in bychrom.items():
  values.sort();starts=[];maximum=[];end=0
  for a,b,_,_ in values:starts.append(a);end=max(end,b);maximum.append(end)
  index[chrom]=(values,starts,maximum)
 # Coordinates establish overlap only; cross-source family homology is not inferred.
 with (work/'FlyBase_overlap.tsv').open('w') as f:
  w=csv.writer(f,delimiter='\t');w.writerow(['FlyBase_copy','FlyBase_family','chrom','covered_bases','annotated_bases','overlap_fraction','overlapping_locus_models','reference_model_groups','family_equivalence_status'])
  for row in rows:
   hits=[];ident=set();families=set()
   for block in row['blocks'].split(','):
    a,b=map(int,block.split('-'))
    values,starts,maximum=index.get(row['chrom'],([],[],[]))
    first=bisect.bisect_left(maximum,a);last=bisect.bisect_right(starts,b)
    for x,y,locus,model in values[first:last]:
     if x<=b and y>=a:hits.append((max(a,x),min(b,y)));ident.add(locus);families.add(model)
   covered=sum(y-x+1 for x,y in interval_union(hits));size=int(row['annotated_bases'])
   w.writerow([row['locus'],row['family'],row['chrom'],covered,size,covered/size,len(ident),';'.join(sorted(families)),'OVERLAP_ONLY_NOT_PROOF_OF_FAMILY_EQUIVALENCE'])

def prepare_annotation(species,scratch):
 compute_only();cfg=PLAN['species'][species];folder,library_models,libreceipt=library(scratch)
 work=scratch/'annotation';work.mkdir();log=work/'prepare.log'
 genome_meta=verify_file(cfg['genome'],always_hash=True)
 rr=json.loads(Path(cfg['reference_receipt']).read_text())
 if canonical(rr)!=cfg['reference_receipt_sha256']:raise ValueError('Mapping reference receipt changed')
 genome=annotation.read_genome(cfg['genome']['path']);lengths={k:len(v) for k,v in genome.items()}
 rows=json.loads((ROOT/'manifests'/f'{species}.json').read_text())
 for row in rows:
  for key in ['bam','bai']:verify_file(row[key])
  if dictionary_from_bam(row['bam']['path'])!=lengths:raise ValueError('BAM dictionary mismatch: '+row['sample_id'])
  if canonical(json.loads(Path(row['mapping_receipt']).read_text()))!=row['mapping_receipt_sha256']:raise ValueError('Mapping receipt changed')
 query=work/'genome.fa';shutil.copyfile(cfg['genome']['path'],query)
 command=['RepeatMasker']+PLAN['repeatmasker_options']+['-lib',folder/'library.fa','-dir',work,query]
 run(command,log)
 native=work/'genome.fa.out'
 if not native.exists() or not (work/'genome.fa.tbl').exists() or not (work/'genome.fa.align').exists():raise ValueError('Required native RepeatMasker output missing')
 catalogue={k:{'family':v['family'],'class':v['class']} for k,v in library_models.items()}
 models,stats=annotation.repeatmasker(native,lengths,species,catalogue)
 annotation.export(models,work);landscape_tables(models,genome,work)
 if species=='Dmel':flybase_comparison(models,work)
 write_json(work/'genome_dictionary.json',lengths);write_json(work/'annotation_validation.json',stats)
 shutil.copy2(native,work/(species+'_repeatmasker.out'))
 shutil.copy2(work/'genome.fa.tbl',work/(species+'_repeatmasker.tbl'))
 with (work/'genome.fa.align').open('rb') as src,gzip.open(work/(species+'_repeatmasker.align.gz'),'wb') as dest:shutil.copyfileobj(src,dest)
 gff=work/'genome.fa.out.gff'
 if gff.exists():shutil.copy2(gff,work/(species+'_repeatmasker.gff'))
 keep=['TE.gtf','locus_catalog.tsv','genome_dictionary.json','annotation_validation.json','prepare.log','genomic_model_coverage.tsv','genomic_class_coverage.tsv','genomic_coverage_summary.json',species+'_repeatmasker.out',species+'_repeatmasker.tbl',species+'_repeatmasker.align.gz']
 if gff.exists():keep.append(species+'_repeatmasker.gff')
 if species=='Dmel':keep.append('FlyBase_overlap.tsv')
 if metadata(cfg['genome']['path'])!={k:genome_meta[k] for k in ['path','bytes','mtime_ns']}:raise ValueError('Genome changed during annotation')
 parent=ROOT/'annotations';parent.mkdir(exist_ok=True);stage=Path(tempfile.mkdtemp(prefix='.'+species+'_',dir=parent))
 for name in keep:shutil.copy2(work/name,stage/name)
 receipt={'state':'ANNOTATION_VALIDATED','scientific_QC':'KNOWN_LIBRARY_ANNOTATION_AND_FAMILY_REVIEW_PENDING','created_utc':utc(),'species':species,'assembly':cfg['assembly'],
  'package_sha256':verify_package(),'job_id':os.environ['SLURM_JOB_ID'],'genome':genome_meta,'annotation_source':{'kind':'SHARED_FROZEN_PAN_DROSOPHILA_LIBRARY','path':str(folder/'library.fa'),'sha256':sha(folder/'library.fa')},
  'shared_library_receipt_sha256':canonical(libreceipt),'command':list(map(str,command)),'validation':stats,'annotation_policy':PLAN['annotation_policy'],
  'files':[dict(name=n,bytes=(stage/n).stat().st_size,mtime_ns=(stage/n).stat().st_mtime_ns,sha256=sha(stage/n)) for n in keep]}
 write_json(stage/'receipt.json',receipt)
 dest=parent/species
 if dest.exists():raise FileExistsError('Existing annotation preserved')
 stage.rename(dest)
