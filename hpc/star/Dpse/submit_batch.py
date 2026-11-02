#!/usr/bin/env python3
"""Submit one reference and independent sample jobs, with durable receipts."""
import hashlib, json, os, re, subprocess
from datetime import datetime, timezone
from pathlib import Path

root = Path(__file__).resolve().parent
c = json.loads((root / 'config.template.json').read_text())
samples = json.loads((root / 'sample_manifest.json').read_text())
assert root == Path(c['run_dir'])
assert len({s['sample_id'] for s in samples}) == len(samples)
receipts = root / 'submissions'
receipts.mkdir(exist_ok=True)

def save_new(path, obj):
    with path.open('x') as f:
        json.dump(obj, f, indent=2)
        f.write('\n'); f.flush(); os.fsync(f.fileno())

def stamp():
    return datetime.now(timezone.utc).isoformat()

for s in samples:
    for item in s['input_metadata']:
        st = Path(item['path']).stat()
        assert (st.st_size, st.st_mtime_ns) == (item['bytes'], item['mtime_ns']), item['path']
    assert not (Path(c['results_dir']) / s['sample_id']).exists()
    assert not (root / 'worker_locks' / (s['sample_id'] + '.lock')).exists()
assert not Path(c['reference_dir']).exists()
save_new(receipts / 'batch_intent.json', {
    'time_utc': stamp(), 'species': c['species'], 'assembly': c['assembly'],
    'samples': len(samples), 'authorization': 'User explicitly requested individual Dpse and Dmir mapping jobs',
    'worker_sha256': hashlib.sha256((root/'mapping_worker.py').read_bytes()).hexdigest(),
    'config_template_sha256': hashlib.sha256((root/'config.template.json').read_bytes()).hexdigest()})

def submit(label, dependency=None):
    script = root / 'jobs' / (label + '.sbatch')
    argv = ['sbatch', '--parsable']
    if dependency:
        argv += ['--dependency=afterok:' + dependency, '--kill-on-invalid-dep=yes']
    argv.append(str(script))
    intent = {'time_utc': stamp(), 'label': label, 'command': argv,
              'script_sha256': hashlib.sha256(script.read_bytes()).hexdigest()}
    save_new(receipts / (label + '.intent.json'), intent)
    p = subprocess.run(argv, text=True, capture_output=True, timeout=60)
    receipt = dict(intent, returncode=p.returncode, stdout=p.stdout, stderr=p.stderr,
                   finished_utc=stamp())
    if p.returncode == 0 and re.fullmatch(r'\d+(;[^\s]+)?', p.stdout.strip()):
        receipt['job_id'] = p.stdout.strip().split(';')[0]
    save_new(receipts / (label + '.receipt.json'), receipt)
    if 'job_id' not in receipt:
        raise RuntimeError('Submission failed or outcome uncertain; inspect saved receipt before retry: '+label)
    print(json.dumps({'label': label, 'job_id': receipt['job_id'], 'dependency': dependency}), flush=True)
    return receipt

all_receipts = [submit('reference')]
for s in samples:
    all_receipts.append(submit(s['sample_id'], all_receipts[0]['job_id']))
save_new(receipts / 'batch_receipt.json', {'finished_utc':stamp(), 'species':c['species'],
    'assembly':c['assembly'], 'jobs':all_receipts})
