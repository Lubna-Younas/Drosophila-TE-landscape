"""Numerical regression: unchanged epsilon/priors, more iterations, valid checkpoint resume."""
import csv,json,re,subprocess,sys
from pathlib import Path
import pysam
out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=True)
gtf=out/'slow.gtf';gtf.write_text('toy\ttest\texon\t101\t200\t.\t+\t.\tlocus "A";\n'+'toy\ttest\texon\t301\t400\t.\t+\t.\tlocus "B";\n')
bam=out/'slow.bam'
with pysam.AlignmentFile(str(bam),'wb',header={'HD':{'VN':'1.6','SO':'queryname'},'SQ':[{'SN':'toy','LN':1000}]}) as f:
 for i in range(1003):
  positions=[120,320] if i<1000 else ([120] if i==1000 else [320])
  for j,pos in enumerate(positions):
   r=pysam.AlignedSegment();r.query_name='r'+str(i);r.query_sequence='A'*20;r.flag=256 if j else 0
   r.reference_id=0;r.reference_start=pos;r.mapping_quality=255;r.cigar=[(0,20)]
   r.set_tag('AS',20);r.set_tag('NH',len(positions));r.set_tag('HI',j+1);f.write(r)
exe=str(Path(sys.executable).parent/'telescope')
common=['--outdir',str(out),'--exp_tag','slow','--reassign_mode','average','--em_epsilon','1e-7','--pi_prior','0','--theta_prior','200000']
with (out/'baseline.log').open('w') as log:subprocess.run([exe,'assign',str(bam),str(gtf),'--stranded_mode','None','--max_iter','1000']+common,check=True,stdout=log,stderr=subprocess.STDOUT)
baseline=(out/'baseline.log').read_text();assert 'EM terminated after 1000 iterations' in baseline,baseline[-1800:]
checkpoints=list(out.glob('*checkpoint.npz'));assert len(checkpoints)==1
with (out/'extended.log').open('w') as log:subprocess.run([exe,'resume',str(checkpoints[0]),'--max_iter','100000']+common,check=True,stdout=log,stderr=subprocess.STDOUT)
text=(out/'extended.log').read_text();match=re.search(r'EM converged after (\d+) iterations',text);assert match,text[-2000:]
assert int(match[1])>1000
with (out/'slow-telescope_report.tsv').open() as f:
 rows=list(csv.DictReader((x for x in f if not x.startswith('#')),delimiter='\t'))
total=sum(float(x['final_count']) for x in rows if x['transcript']!='__no_feature')
assert abs(total-1003)<0.02,total
result={'state':'NUMERICAL_RECOVERY_VALIDATED','baseline_iterations':1000,'baseline_converged':False,'extended_converged_iterations':int(match[1]),'em_epsilon':1e-7,'pi_prior':0,'theta_prior':200000,'total_fragments':total,'checkpoint_resume_passed':True}
(out/'checks.json').write_text(json.dumps(result,indent=2));print(json.dumps(result))
