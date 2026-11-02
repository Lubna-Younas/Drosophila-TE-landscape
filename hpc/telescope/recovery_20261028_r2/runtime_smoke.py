"""Tiny counting regressions for both strands, layouts, and sense/antisense controls."""
import csv,itertools,json,subprocess,sys
from pathlib import Path
import pysam
out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=True)
checks=[]
for layout,mode in [('SE','F'),('SE','R'),('SE','None'),('PE','FR'),('PE','RF'),('PE','None')]:
 for strand,antisense in itertools.product(['+','-'],[False,True]):
  name=layout+'_'+mode+'_'+('plus' if strand=='+' else 'minus')+'_'+('antisense' if antisense else 'sense')
  folder=out/name;folder.mkdir(exist_ok=True);bam=folder/'toy.bam';gtf=folder/'toy.gtf'
  gtf.write_text('toy\ttest\texon\t101\t200\t.\t'+strand+'\t.\tlocus "A";\n'+'toy\ttest\texon\t301\t400\t.\t'+strand+'\t.\tlocus "B";\n')
  reverse=(mode in ('R','RF')) ^ (strand=='-') ^ antisense
  with pysam.AlignmentFile(str(bam),'wb',header={'HD':{'VN':'1.6','SO':'queryname'},'SQ':[{'SN':'toy','LN':1000}]}) as f:
   for qname,positions in [('one',[120]),('three',[120,320]),('two',[320])]:
    for i,pos in enumerate(positions):
     for mate in ((1,2) if layout=='PE' else (1,)):
      r=pysam.AlignedSegment();r.query_name=qname;r.query_sequence='A'*20
      if layout=='PE':
       r.flag=((83 if mate==1 else 163) if reverse else (99 if mate==1 else 147))+(256 if i else 0)
       r.reference_id=0;r.reference_start=pos+(40 if ((mate==1)==reverse) else 0)
       r.next_reference_id=0;r.next_reference_start=pos+(0 if ((mate==1)==reverse) else 40)
       r.template_length=60 if r.reference_start<r.next_reference_start else -60
      else:
       r.flag=(16 if reverse else 0)+(256 if i else 0);r.reference_id=0;r.reference_start=pos
      r.mapping_quality=255;r.cigar=[(0,20)];r.set_tag('AS',20);r.set_tag('NH',len(positions));r.set_tag('HI',i+1);f.write(r)
  with (folder/'smoke.log').open('w') as log:
   subprocess.run([str(Path(sys.executable).parent/'telescope'),'assign',str(bam),str(gtf),
                   '--outdir',str(folder),'--exp_tag','toy','--stranded_mode',mode,
                   '--reassign_mode','average','--max_iter','1000'],check=True,stdout=log,stderr=subprocess.STDOUT)
  expected=0 if mode!='None' and antisense else 3
  report=folder/'toy-telescope_report.tsv';logtext=(folder/'smoke.log').read_text();counts={}
  assert 'telescope assign complete' in logtext
  if expected:
   with report.open() as f:rows=list(csv.DictReader((x for x in f if not x.startswith('#')),delimiter='\t'))
   counts={r['transcript']:float(r['final_count']) for r in rows if r['transcript']!='__no_feature'}
   assert set(counts)=={'A','B'} and all(abs(v-1.5)<0.01 for v in counts.values()),(name,counts)
   assert 'EM converged' in logtext
  else:
   assert not report.exists() and 'No alignments overlapping annotation' in logtext,name
  checks.append({'layout':layout,'mode':mode,'annotation_strand':strand,'antisense_control':antisense,
                 'expected_total':expected,'counts':counts,'passed':True})
(out/'checks.json').write_text(json.dumps(checks,indent=2))
print(json.dumps({'passed_cases':len(checks),'details':checks}))
