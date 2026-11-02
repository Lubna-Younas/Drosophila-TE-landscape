import hashlib,json,os,socket,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parent
PLAN=json.loads((ROOT/'plan.json').read_text())
def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def canonical(obj):return hashlib.sha256(json.dumps(obj,sort_keys=True).encode()).hexdigest()
def metadata(path):
 p=Path(path);s=p.stat();return {'path':str(p),'bytes':s.st_size,'mtime_ns':s.st_mtime_ns}
def write_json(path,obj):
 path=Path(path);tmp=path.with_suffix(path.suffix+'.tmp');tmp.write_text(json.dumps(obj,indent=2)+'\n');os.replace(tmp,path)
def verify_file(spec,always_hash=False):
 current=metadata(spec['path'])
 if current['bytes']!=spec['bytes']:raise ValueError('Size changed: '+spec['path'])
 if always_hash or current['mtime_ns']!=spec['mtime_ns']:
  if sha(spec['path'])!=spec['sha256']:raise ValueError('Checksum changed: '+spec['path'])
 return dict(current,sha256=spec['sha256'])
def compute_only():
 if not os.environ.get('SLURM_JOB_ID') or socket.gethostname().startswith('login'):
  raise SystemExit('Run this worker through sbatch on a compute node.')
def isolated_env():
 env=dict(os.environ)
 for key in ('PYTHONPATH','PYTHONHOME','LD_LIBRARY_PATH'):env.pop(key,None)
 env['PYTHONNOUSERSITE']='1'
 # Keep workflow caches and Conda bookkeeping off the quota-limited home.
 # HOME and the user's shared Conda configuration are deliberately unchanged.
 cache=ROOT/'runtime/cache'
 env.update(XDG_CACHE_HOME=str(cache/'xdg'),
            CONDA_PKGS_DIRS=str(cache/'conda-pkgs'),
            CONDA_ENVS_PATH=str(ROOT/'runtime/envs'),
            PIP_CACHE_DIR=str(cache/'pip'),
            CONDA_REGISTER_ENVS='false',CONDA_NUMBER_CHANNEL_NOTICES='0',
            CONDA_NOTIFY_OUTDATED_CONDA='false')
 return env
def verify_package():
 data=json.loads((ROOT/'package_checksums.json').read_text())
 for name,digest in data.items():
  if sha(ROOT/name)!=digest:raise ValueError('Package changed: '+name)
 return canonical(data)
def run(args,log,cwd=None):
 with Path(log).open('a') as f:
  f.write('\nCOMMAND '+json.dumps(list(map(str,args)))+'\n');f.flush()
  executable=str(args[0])
  env=isolated_env() if '/runtime/' in executable or executable==PLAN['conda'] else None
  subprocess.run(list(map(str,args)),check=True,stdout=f,stderr=subprocess.STDOUT,cwd=cwd,env=env)
def dictionary_from_bam(path):
 output=subprocess.check_output(['samtools','view','-H',str(path)],text=True)
 result={}
 for line in output.splitlines():
  if line.startswith('@SQ\t'):
   row=dict(c.split(':',1) for c in line.split('\t')[1:]);result[row['SN']]=int(row['LN'])
 if not result:raise ValueError('Empty BAM reference dictionary')
 return result
