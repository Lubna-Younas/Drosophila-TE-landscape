#!/usr/bin/env python3
"""Static/publication checks only: never import HPC workers or submit jobs."""
import ast
import collections
import csv
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[1]
SKIP_DIRS = {'.git', '__pycache__', '.venv'}
FIGURE_PACKAGE = 'analysis/gene_TE_density_circos_20261103'


def repository_files():
    return sorted(
        p for p in ROOT.rglob('*')
        if p.is_file() and not any(part in SKIP_DIRS for part in p.relative_to(ROOT).parts)
        and not p.name.startswith('._') and p.name != '.DS_Store'
    )


def read_tsv(relative):
    with (ROOT / relative).open(encoding='utf-8', newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def check_sources():
    rows = read_tsv('metadata/source_manifest.tsv')
    paths = [r['public_path'] for r in rows]
    assert len(paths) == len(set(paths)), 'Duplicate source manifest paths'
    for row in rows:
        path = ROOT / row['public_path']
        assert path.is_file(), 'Missing source: ' + row['public_path']
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        assert digest == row['published_sha256'], 'Published hash changed: ' + row['public_path']
        assert re.fullmatch('[0-9a-f]{64}', row['source_sha256']), 'Invalid original source hash'
        if row['transformation'] == 'unchanged':
            assert row['source_sha256'] == digest, 'Unrecorded source transformation'
    verification = json.loads((ROOT / 'metadata/source_verification.json').read_text())
    expected = {'Dmel': 49, 'Dmir': 22, 'Dpse': 53}
    assert {r['species'] for r in verification['species']} == set(expected)
    for item in verification['species']:
        sp = item['species']
        assert item['matching_selected_sample_receipts'] == expected[sp]
        assert item['all_selected_receipts_match_worker'] is True
        worker = next(r for r in rows if r['public_path'] == 'hpc/star/' + sp + '/mapping_worker.py')
        assert worker['source_sha256'] == item['mapping_worker_sha256']
    return len(rows)


def check_metadata():
    rows = read_tsv('metadata/samples.tsv')
    keys = [(r['species'], r['sample_id']) for r in rows]
    assert len(keys) == len(set(keys)), 'Duplicate sample identities'
    assert collections.Counter(r['species'] for r in rows) == {'Dmel': 50, 'Dmir': 22, 'Dpse': 53}
    selected = [r for r in rows if r['included_in_TE_cohort'] == 'True']
    assert collections.Counter(r['species'] for r in selected) == {'Dmel': 49, 'Dmir': 22, 'Dpse': 53}
    assert len(selected) == 124
    excluded = read_tsv('metadata/excluded_samples.tsv')
    assert {r['run_accession'] for r in excluded} == {'SRR6056856'}
    assert 'SRR6056856' not in {r['run_accession'] for r in selected}
    for row in rows:
        assert row['layout'] in ('PE', 'SE'), 'Missing verified layout'
        assert bool(row['read1_basename'])
        assert bool(row['read2_basename']) == (row['layout'] == 'PE')
        assert '/' not in row['read1_basename'] and '/' not in row['read2_basename']
    for row in selected:
        assert row['telescope_stranded_mode'] in ('None', 'RF', 'FR', 'F', 'R')
    refs = read_tsv('metadata/references.tsv')
    assert len(refs) == 3
    dpse = next(r for r in refs if r['species'] == 'Dpse')
    assert dpse['assembly'] == 'GCA_004329205.1_UCBerk_Dpse_1.0'
    assert not dpse['gene_gtf_basename'] and not dpse['gene_gtf_sha256']
    return len(rows), len(selected)


def check_figure_package(files):
    """Allow only explicitly inventoried, checksum-verified figure artifacts."""
    manifest_path=ROOT/FIGURE_PACKAGE/'provenance/publication_manifest.json'
    manifest=json.loads(manifest_path.read_text())
    assert manifest['state']=='VALIDATED_PUBLICATION'
    records={r['file']:r for r in manifest['files']}
    assert len(records)==len(manifest['files'])
    observed={str(p.relative_to(ROOT/FIGURE_PACKAGE)) for p in files
              if (ROOT/FIGURE_PACKAGE) in p.parents and p!=manifest_path}
    assert set(records)==observed,'Figure package inventory mismatch'
    allowed={}
    for rel,rec in records.items():
        path=ROOT/FIGURE_PACKAGE/rel
        assert (ROOT/FIGURE_PACKAGE).resolve() in path.resolve().parents
        assert path.stat().st_size==rec['bytes']<20_000_000
        assert hashlib.sha256(path.read_bytes()).hexdigest()==rec['published_sha256'],rel
        allowed[str(path.relative_to(ROOT))]=rec
    return allowed


def check_files(files, figure_artifacts):
    counts = collections.Counter()
    forbidden_suffixes = ('.fastq', '.fastq.gz', '.fq', '.fq.gz', '.bam', '.bai',
                          '.cram', '.crai', '.sra', '.h5', '.hdf5', '.npz',
                          '.fasta', '.fasta.gz', '.fna', '.fna.gz', '.fa', '.fa.gz',
                          '.gtf', '.gtf.gz', '.gff', '.gff.gz', '.align', '.align.gz')
    # Generic patterns avoid embedding private site details in this checker.
    sensitive = [
        re.compile(r'/(?:Users|Volumes)/[^\s/]+/'),
        re.compile(r'/lisc/data/'),
        re.compile(r'/gpfs/data/'),
        re.compile(r'-----BEGIN (?:RSA |OPENSSH |EC )?PRIVATE KEY-----'),
        re.compile(r'\bgh[pousr]_[A-Za-z0-9]{20,}\b'),
        re.compile(r'\bgithub_pat_[A-Za-z0-9_]{20,}\b'),
    ]
    for path in files:
        rel = str(path.relative_to(ROOT))
        assert path.stat().st_size < 2_000_000 or rel in figure_artifacts, 'Unexpected large publication file: ' + rel
        assert not path.name.endswith(forbidden_suffixes), 'Data file in publication: ' + rel
        if path.suffix in ('.pdf','.png'):
            assert rel in figure_artifacts,'Uninventoried binary artifact: '+rel
            raw=path.read_bytes()
            if path.suffix=='.pdf':
                assert raw.startswith(b'%PDF-') and b'%%EOF' in raw[-1024:]
            else:
                assert raw.startswith(b'\x89PNG\r\n\x1a\n')
            counts['verified_figure_binary']+=1
            continue
        text = path.read_text(encoding='utf-8')
        # The patterns themselves are not private values.
        if path != Path(__file__).resolve():
            for pattern in sensitive:
                assert not pattern.search(text), 'Possible private path/credential: ' + rel
        if path.suffix == '.py':
            ast.parse(text, filename=rel)
            counts['python'] += 1
        elif path.suffix in ('.sh', '.sbatch'):
            subprocess.run(['bash', '-n', str(path)], check=True, capture_output=True)
            counts['shell'] += 1
        elif path.suffix == '.json':
            json.loads(text)
            counts['json'] += 1
        elif path.suffix == '.md':
            for target in re.findall(r'\[[^\]]*\]\(([^\s)]+)\)', text):
                if urlsplit(target).scheme or target.startswith('#'):
                    continue
                link = unquote(target.split('#', 1)[0]).strip('<>')
                assert (path.parent / link).exists(), 'Broken link in ' + rel + ': ' + target
            counts['markdown'] += 1
    return counts


def main():
    sources = check_sources()
    inventory, selected = check_metadata()
    files = repository_files()
    figure_artifacts=check_figure_package(files)
    counts = check_files(files, figure_artifacts)
    print(json.dumps({
        'state': 'PASS_STATIC_AND_METADATA_CHECKS',
        'files': len(files), 'source_hashes_verified': sources,
        'figure_artifact_hashes_verified':len(figure_artifacts),
        'inventory_libraries': inventory, 'selected_TE_libraries': selected,
        'syntax_and_documentation': dict(counts),
        'scientific_execution': 'NOT_RUN; no HPC submission or biological QC claim',
    }, indent=2))


if __name__ == '__main__':
    try:
        main()
    except (AssertionError, ValueError, OSError, subprocess.CalledProcessError) as exc:
        print('VALIDATION FAILED: ' + str(exc), file=sys.stderr)
        raise SystemExit(1)
