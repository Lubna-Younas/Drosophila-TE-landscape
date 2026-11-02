"""Compute-node validation for unknown repeat labels and RepeatMasker completion."""
import json, os, re, shutil
from pathlib import Path
import annotation
from common import run, sha, write_json, verify_package

def validate_repeatmasker_output(output, text):
    fatal = re.search(r'Quantifier follows nothing|Traceback \(most recent call last\)|'
                      r'Unmatched .* in regex|Execution of .* aborted', text)
    if fatal:
        raise RuntimeError('RepeatMasker internal error: ' + fatal.group(0))
    if not output.is_file() or output.stat().st_size == 0:
        raise RuntimeError('RepeatMasker did not produce a nonempty output: ' + str(output))

def checked_repeatmasker(library, query, work, log):
    start = log.stat().st_size if log.exists() else 0
    run(['RepeatMasker','-engine','rmblast','-pa','1','-no_is',
         '-lib',library,'-dir',work,query], log)
    with log.open('rb') as f:
        f.seek(start)
        text = f.read().decode(errors='replace')
    output = work/(query.name+'.out')
    validate_repeatmasker_output(output, text)
    return output

def unknown_label_smoke(source, work, destination):
    work.mkdir()
    catalog = annotation.normalize_library(source, work/'library.fa')
    records = list(annotation.fasta(work/'library.fa'))
    unknown = [(h,s) for h,s in records if h.endswith('#Unknown')]
    known = [(h,s) for h,s in records if not h.endswith('#Unknown')][:1]
    assert len(unknown)==74 and len(known)==1
    query = work/'unknown_queries.fa'
    lengths = {}
    with query.open('w') as f:
        for header,seq in unknown+known:
            name='query_'+header.split('#')[0]
            lengths[name]=len(seq)
            f.write('>'+name+'\n'+seq+'\n')
    output = checked_repeatmasker(work/'library.fa',query,work,work/'smoke.log')
    models,stats=annotation.repeatmasker(output,lengths,'Dpse_smoke',catalog)
    unknown_models=[v for v in models.values() if v['class']=='Unknown']
    if not unknown_models:
        raise RuntimeError('Unknown-label smoke produced no retained Unknown annotation')
    annotation.export(models,work)
    receipt={'state':'UNKNOWN_LABEL_SMOKE_PASSED','unknown_query_records':74,
             'known_query_controls':1,'query_bases':sum(lengths.values()),
             'retained_unknown_loci':len(unknown_models),'annotation_statistics':stats,
             'package_sha256':verify_package(),'job_id':os.environ['SLURM_JOB_ID'],
             'library_sha256':sha(work/'library.fa'),'query_sha256':sha(query),
             'output_sha256':sha(output)}
    destination.mkdir(parents=True,exist_ok=False)
    for name in ('smoke.log','unknown_queries.fa.out','TE.gtf','locus_catalog.tsv',
                 'library_classification_provenance.tsv','library_validation.json'):
        shutil.copy2(work/name,destination/name)
    write_json(destination/'receipt.json',receipt)
    return receipt
