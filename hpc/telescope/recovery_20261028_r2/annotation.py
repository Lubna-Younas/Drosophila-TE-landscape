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
 catalog={}
 with open(out,'w') as f:
  for header,seq in fasta(path):
   fields=header.split();name=fields[0]
   if len(fields)<4 or name in catalog or not re.fullmatch('[ACGTRYMKSWBDHVN]+',seq):raise ValueError('Invalid consensus record')
   family,subclass,cl=fields[1:4]
   catalog[name]={'family':family,'class':cl+'/'+subclass}
   f.write('>'+name+'#'+cl+'/'+subclass+'\n'+seq+'\n')
 if len(catalog)!=303:raise ValueError('Unexpected consensus library record count')
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
