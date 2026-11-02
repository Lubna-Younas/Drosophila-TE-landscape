"""Exact-reference TE models, without conversion of the original inputs."""
import collections,csv,gzip,hashlib,json,re
from pathlib import Path

def fasta(path):
 opener=gzip.open if str(path).endswith('.gz') else open
 header=None;seq=[]
 with opener(path,'rt') as f:
  for line in f:
   if line.startswith('>'):
    if header is not None:yield header,''.join(seq).upper()
    header=line[1:].strip();seq=[]
   elif line.strip():seq.append(line.strip())
  if header is not None:yield header,''.join(seq).upper()

def read_genome(path):
 result={}
 for header,seq in fasta(path):
  name=header.split()[0]
  if name in result or not seq:raise ValueError('Empty or duplicate genome record: '+name)
  result[name]=seq
 return result

def flybase(path,genome):
 models={};trans=str.maketrans('ACGTRYMKBDHVN','TGCAYRKMVHDBN')
 for header,seq in fasta(path):
  copy=header.split()[0]; attrs=dict(re.findall(r'(\w+)=([^;]+);',header))
  if copy in models or not copy.startswith('FBti'):raise ValueError('Invalid FlyBase copy identifier')
  if attrs.get('release')!='r6.69':raise ValueError('Wrong FlyBase release')
  m=re.fullmatch(r'([^:]+):(complement\()?([0-9]+)\.\.([0-9]+)\)?',attrs['loc'])
  if not m:raise ValueError('Unsupported coordinate syntax: '+attrs['loc'])
  chrom,rev,a,b=m.groups();a,b=int(a),int(b);strand='-' if rev else '+'
  if chrom not in genome or not 1<=a<=b<=len(genome[chrom]):raise ValueError('TE outside genome')
  expected=genome[chrom][a-1:b]
  if rev:expected=expected.translate(trans)[::-1]
  if expected!=seq or len(seq)!=int(attrs['length']) or hashlib.md5(seq.encode()).hexdigest()!=attrs['MD5']:
   raise ValueError('TE sequence or header checksum mismatch: '+copy)
  models[copy]={'locus':copy,'chrom':chrom,'strand':strand,'blocks':[(a,b)],
                'family':attrs['name'].split('{}')[0],'class':'UNCLASSIFIED','source_id':copy}
 return models,{'sequence_verified_copies':len(models)}

def normalize_library(path,out):
 catalog={};provenance=[];source_records=list(fasta(path))
 for header,seq in source_records:
  fields=header.split();name=fields[0]
  if len(fields)<4 or name in catalog or not re.fullmatch('[ACGTRYMKSWBDHVN]+',seq):
   raise ValueError('Invalid consensus record: '+name)
  family,subclass,cl=fields[1:4]
  original_class=cl+'/'+subclass
  if cl=='?' and subclass=='?':normalized='Unknown'
  else:normalized=original_class
  if not re.fullmatch(r'[A-Za-z0-9_.-]+(?:/[A-Za-z0-9_.-]+)?',normalized):
   raise ValueError('Unsafe or unsupported classification: '+original_class)
  catalog[name]={'family':family,'class':normalized,'source_class':original_class}
  provenance.append({'sequence_id':name,'family':family,'original_class':cl,
    'original_subclass':subclass,'original_combined_class':original_class,
    'repeatmasker_class':normalized,'sequence_length':len(seq),
    'sequence_sha256':hashlib.sha256(seq.encode()).hexdigest(),'original_header':header})
 if len(catalog)!=303:raise ValueError('Unexpected consensus library record count')
 if sum(x['original_combined_class']=='?/?' for x in provenance)!=74:
  raise ValueError('Unexpected count of unknown source classifications')
 out=Path(out)
 with out.open('w') as f:
  for (header,seq),row in zip(source_records,provenance):
   f.write('>'+row['sequence_id']+'#'+row['repeatmasker_class']+'\n'+seq+'\n')
 normalized_records=list(fasta(out))
 if len(normalized_records)!=len(source_records):raise ValueError('Sequence record loss')
 for (old_header,old_seq),(new_header,new_seq) in zip(source_records,normalized_records):
  if old_header.split()[0]!=new_header.split('#')[0] or old_seq!=new_seq:
   raise ValueError('Consensus sequence/identifier changed')
 with (out.parent/'library_classification_provenance.tsv').open('w') as f:
  w=csv.DictWriter(f,list(provenance[0]),delimiter='\t');w.writeheader();w.writerows(provenance)
 validation={'state':'ALL_SEQUENCES_PRESERVED','records':len(catalog),'unknown_labels_corrected':74,
   'known_labels_unchanged':229,'sequence_bases':sum(len(s) for _,s in source_records),
   'source_file_sha256':hashlib.sha256(Path(path).read_bytes()).hexdigest(),
   'normalized_file_sha256':hashlib.sha256(out.read_bytes()).hexdigest(),
   'sequence_identity_check':'Ordered IDs and full nucleotide strings identical after FASTA parsing'}
 (out.parent/'library_validation.json').write_text(json.dumps(validation,indent=2)+'\n')
 return catalog

def repeatmasker(path,lengths,species,library_catalog=None):
 models={};excluded=collections.Counter();rows=0
 opener=gzip.open if str(path).endswith('.gz') else open
 non_te={'Simple_repeat','Low_complexity','Satellite','RNA','rRNA','scRNA','snRNA','srpRNA','tRNA'}
 with opener(path,'rt') as f:
  for line in f:
   v=line.split()
   if not v or not v[0].isdigit():continue
   if len(v) not in (15,16) or (len(v)==16 and v[-1]!='*'):raise ValueError('Unexpected RepeatMasker row')
   chrom,a,b,strand,name,cl,ident=v[4],int(v[5]),int(v[6]),v[8],v[9],v[10],v[14]
   rows+=1
   if chrom not in lengths or not 1<=a<=b<=lengths[chrom]:raise ValueError('RepeatMasker/reference mismatch')
   if cl.split('/')[0] in non_te:excluded[cl]+=1;continue
   if strand not in ('C','+'):raise ValueError('Unknown RepeatMasker strand')
   strand='-' if strand=='C' else '+'
   family=name
   if library_catalog is not None:
    if name not in library_catalog:raise ValueError('Unknown custom consensus name: '+name)
    family=library_catalog[name]['family'];cl=library_catalog[name]['class']
   # One RepeatMasker ID can have disjoint blocks; nested gaps are not filled.
   key=(chrom,ident,name,strand)
   locus=species+'_RM_'+hashlib.sha256('|'.join(key).encode()).hexdigest()[:20]
   model=models.setdefault(locus,{'locus':locus,'chrom':chrom,'strand':strand,'blocks':[],
       'family':family,'class':cl,'source_id':ident})
   model['blocks'].append((a,b))
 if not rows or not models:raise ValueError('No usable TE annotations')
 for m in models.values():
  merged=[]
  for a,b in sorted(m['blocks']):
   if merged and a<=merged[-1][1]+1:merged[-1]=(merged[-1][0],max(merged[-1][1],b))
   else:merged.append((a,b))
  m['blocks']=merged
 return models,{'source_rows':rows,'excluded_non_TE_rows':dict(excluded),'retained_TE_loci':len(models)}

def export(models,out):
 if not models:raise ValueError('Empty annotation')
 out=Path(out)
 with (out/'TE.gtf').open('w') as g,(out/'locus_catalog.tsv').open('w') as c:
  fields=['locus','family','class','chrom','strand','source_id','annotated_bases','blocks']
  w=csv.DictWriter(c,fields,delimiter='\t');w.writeheader()
  for m in sorted(models.values(),key=lambda x:(x['chrom'],x['blocks'][0][0],x['locus'])):
   if any('"' in str(m[k]) or '\t' in str(m[k]) for k in ('locus','family','class')):raise ValueError('Unsafe GTF attribute')
   for a,b in m['blocks']:
    attr='gene_id "{0}"; transcript_id "{0}"; locus "{0}"; family "{1}"; class "{2}";'.format(m['locus'],m['family'],m['class'])
    g.write('\t'.join(map(str,[m['chrom'],'TE_catalog','exon',a,b,'.',m['strand'],'.',attr]))+'\n')
   w.writerow({**{k:m[k] for k in fields if k in m},'annotated_bases':sum(b-a+1 for a,b in m['blocks']),
               'blocks':','.join(str(a)+'-'+str(b) for a,b in m['blocks'])})
