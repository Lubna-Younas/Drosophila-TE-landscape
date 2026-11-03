#!/usr/bin/env python3
"""Calculate exclusive TE-class coverage from the frozen TE locus blocks.

Small local annotation-table computation; no genome/read scanning, reannotation,
or read mapping. Cross-class overlap is retained as Mixed rather than double counted.
"""
import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
CLASSES=['LTR','LINE','DNA','RC','Unknown','Mixed']
WINDOW=100000

def read(path):
    with path.open(newline='') as f:return list(csv.DictReader(f,delimiter='\t'))

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--species',nargs='+',choices=['Dmel','Dpse','Dmir'],default=['Dmel','Dpse','Dmir'])
    args=p.parse_args()
    receipt=json.loads((ROOT/'source/export_receipt.json').read_text())
    cross={r['original_library_class']:r['TE_class'] for r in read(ROOT/'source/classification_crosswalk.tsv')}
    records=[]
    for sp in args.species:
        cat=ROOT/'source'/f'{sp}_locus_catalog.tsv'
        expected=next(r for r in receipt['species'] if r['species']==sp)['source_checks']['locus_catalog.tsv']['sha256']
        assert cat.exists(),f'Matching {sp} catalogue required: {cat}'
        assert sha(cat)==expected, f'{sp} catalogue does not match current genome/TE annotation receipt'
        windows=read(ROOT/'source'/f'{sp}_windows_100kb.tsv')
        seqs={r['accession']:r for r in read(ROOT/'source'/f'{sp}_sequences.tsv')}
        counts={(r['accession'],int(r['start0'])):Counter() for r in windows}
        bychrom=defaultdict(list)
        loci=read(cat)
        assert len({r['locus'] for r in loci})==len(loci)
        for r in loci:bychrom[r['chrom']].append(r)
        assert set(bychrom)<=set(seqs)
        totals=defaultdict(Counter)
        for chrom,rows in bychrom.items():
            events=defaultdict(Counter)
            for r in rows:
                cl=cross[r['class']]
                blocks=[tuple(map(int,b.split('-'))) for b in r['blocks'].split(',')]
                assert sum(b-a+1 for a,b in blocks)==int(r['annotated_bases'])
                assert all(1<=a<=b<=int(seqs[chrom]['span_bp']) for a,b in blocks)
                for a,b in blocks:
                    events[a-1][cl]+=1
                    events[b][cl]-=1
            active=Counter(); prev=0
            for pos,delta in sorted(events.items()):
                positive=[c for c,n in active.items() if n>0]
                if positive and pos>prev:
                    cl=positive[0] if len(positive)==1 else 'Mixed'
                    totals[chrom][cl]+=pos-prev
                    for i in range(prev//WINDOW,(pos-1)//WINDOW+1):
                        amount=min(pos,(i+1)*WINDOW)-max(prev,i*WINDOW)
                        counts[chrom,i*WINDOW][cl]+=amount
                active.update(delta)
                assert all(v>=0 for v in active.values())
                prev=pos
            assert not any(active.values())
        original_class_totals=defaultdict(Counter)
        for r in read(ROOT/'source'/f'{sp}_TE_superfamilies.tsv'):
            original_class_totals[r['accession']][r['TE_class']]+=int(r['exclusive_bp'])
        for chrom in seqs:
            assert totals[chrom]==original_class_totals[chrom],(sp,chrom,'class totals mismatch')
        output=[]
        for r in windows:
            cv=counts[r['accession'],int(r['start0'])]
            assert sum(cv.values())==int(r['TE_bp']), (sp,r['accession'],r['start0'])
            row={k:r[k] for k in ['species','accession','start0','end0','span_bp','TE_bp']}
            row.update({f'{c}_bp':cv[c] for c in CLASSES})
            output.append(row)
        dest=ROOT/'tables'/f'{sp}_TE_classes_100kb.tsv'
        with dest.open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(output[0]),delimiter='\t');w.writeheader();w.writerows(output)
        info={'species':sp,'state':'VALIDATED','catalogue_sha256':expected,
              'classification_sha256':sha(ROOT/'source/classification_crosswalk.tsv'),
              'windows':len(windows),'TE_loci':len(loci),'output_sha256':sha(dest),
              'every_window_matches_existing_TE_union':True,
              'every_chromosome_class_total_matches_existing_export':True,
              'classes':CLASSES,'cross_class_overlap':'Mixed; counted once',
              'ring_denominator':'TE-covered bases in each window; white if no TE coverage'}
        (ROOT/'provenance'/f'{sp}_class_window_receipt.json').write_text(json.dumps(info,indent=2)+'\n')
        records.append(info)
    print(json.dumps(records,indent=2))

if __name__=='__main__':main()
