#!/usr/bin/env python3
"""Run existing small tests without private inputs, installed NGS tools, or jobs."""
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
DMEL_CASES = [
    'test_annotation_compatibility_preserves_rows',
    'test_strand_classification_and_ambiguity',
    'test_spliced_alignment_blocks',
    'test_read_pair_validation',
    'test_no_login_node_execution',
    'test_publication_checksums_and_no_overwrite',
    'test_submission_dry_run_and_duplicate_prevention',
]


def run(directory, arguments):
    env = dict(os.environ)
    # The no-login-node test must see an ordinary local environment.
    env.pop('SLURM_JOB_ID', None)
    env['PYTHONDONTWRITEBYTECODE'] = '1'
    subprocess.run([sys.executable] + arguments, cwd=directory, env=env, check=True)


def main():
    run(ROOT / 'hpc/star/Dmel', ['-m', 'unittest', '-v'] + [
        'test_mapping_code.MappingTests.' + name for name in DMEL_CASES])
    run(ROOT / 'hpc/star/Dmir', ['test_worker.py'])
    run(ROOT / 'hpc/star/Dpse', ['-m', 'unittest', '-v', 'test_worker', 'test_genome_only'])
    # Imported annotation modules expect plan.json. Supply the example only in
    # temporary storage; these tests call pure synthetic helpers, not workers.
    with tempfile.TemporaryDirectory(prefix='te_source_checks_') as directory:
        tmp = Path(directory)
        source = ROOT / 'hpc/telescope/shared_20261029'
        for name in ('common.py', 'annotation.py', 'landscape_annotation.py', 'technical_checks.py'):
            shutil.copy2(source / name, tmp / name)
        shutil.copy2(source / 'plan.example.json', tmp / 'plan.json')
        run(tmp, ['technical_checks.py'])
    print('PASS: 32 existing synthetic tests; no scientific executables or HPC jobs run.')


if __name__ == '__main__':
    main()
