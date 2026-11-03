#!/usr/bin/env python3
"""Assemble both already-rendered scientific figures into comparison artifacts."""
import argparse
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--images',action='store_true')
    p.add_argument('--pdf',action='store_true')
    args=p.parse_args()
    both=not (args.images or args.pdf)
    if args.images or both:
        from PIL import Image
        panels=[]
        for name in ['A_composition_preview.png','B_dominant_preview.png']:
            with Image.open(ROOT/'figures'/name) as img:
                img=img.convert('RGB')
                panels.append(img.resize((1280,round(img.height*1280/img.width)),Image.Resampling.LANCZOS))
        canvas=Image.new('RGB',(sum(x.width for x in panels),max(x.height for x in panels)),'white')
        x=0
        for panel in panels:
            canvas.paste(panel,(x,0));x+=panel.width
        canvas.save(ROOT/'figures/compare_both_designs.png')
    if args.pdf or both:
        import pypdf
        import pypdf.filters
        # Bounded decompression for these trusted, locally generated vector plots.
        pypdf.filters.ZLIB_MAX_OUTPUT_LENGTH=150_000_000
        writer=pypdf.PdfWriter()
        validation=[]
        for name in ['A_composition_circos.pdf','B_dominant_circos.pdf']:
            path=ROOT/'figures'/name
            reader=pypdf.PdfReader(path)
            assert len(reader.pages)==1
            page=reader.pages[0]
            stream=page['/Contents'].get_object().get_data()
            assert 0<len(stream)<150_000_000
            assert list(page.mediabox)==[0,0,1296,1404]
            validation.append({'file':name,'pages':1,'decoded_content_bytes':len(stream),
                               'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
            writer.add_page(page)
        writer.add_metadata({'/Title':'Gene and TE density: two solid-class ring designs',
                             '/Subject':'A: chromosome composition. B: dominant class by position. Shared 0-1 density scales.'})
        out=ROOT/'figures/both_circos_designs.pdf'
        with out.open('wb') as f:writer.write(f)
        assert len(pypdf.PdfReader(out).pages)==2
        (ROOT/'provenance/pdf_validation.json').write_text(json.dumps(
            {'state':'PASS','individual_figures':validation,'comparison_pages':2,
             'comparison_sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'pypdf':pypdf.__version__},indent=2)+'\n')
    print('Comparison artifact assembly complete')

if __name__=='__main__':main()
