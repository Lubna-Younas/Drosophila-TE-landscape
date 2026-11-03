#!/usr/bin/env python3
"""Render two twelve-track Circos variants from audited 100 kb density tables.

No remapping, interval recomputation, smoothing or scaffold-order inference.
Raw values stay unchanged; densities are explicitly rescaled to common 0-1 axes.
"""
import argparse
import csv
import hashlib
import json
import math
import shutil
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SPECIES = ['Dmir', 'Dpse', 'Dmel']
NAMES = {'Dmel': 'D. melanogaster', 'Dpse': 'D. pseudoobscura', 'Dmir': 'D. miranda'}
COLORS = {'chromosome': '#B8CCD5', 'gene': '#D8B365', 'TE': '#5AA89A',
          'ink': '#283238', 'grid': '#DFE5E8', 'boundary': '#7C858B'}
ORDER = {'Dmel': ['X', '3L', '2L', '2R', '3R', '4', 'Y'],
         'Dpse': ['X (XL + XR)', '4', '3', '2', 'Dot', 'Y'],
         'Dmir': ['XL', 'XR', '4', 'neo-X (3)', '2', '5', 'Y scaffolds']}
LABELS = {'Dmel': ['X · A', '3L · D', '2L · B', '2R · C', '3R · E', '4 · F', 'Y'],
          'Dpse': ['X · A+D', '4 · B', '3 · C', '2 · E', 'Dot · F', 'Y'],
          'Dmir': ['XL · A', 'XR · D*', '4 · B', 'neo-X · C', '2 · E', '5 · F', 'Y scaffolds†']}
# Four tracks per species; shared density scales and shared class palette.
RADII = {'Dmir': {'chr': (1.82,1.845),'gene':(1.695,1.775),'TE':(1.575,1.655),'class':(1.47,1.53)},
         'Dpse': {'chr': (1.34,1.365),'gene':(1.215,1.295),'TE':(1.095,1.175),'class':(.99,1.05)},
         'Dmel': {'chr': (.86,.885),'gene':(.735,.815),'TE':(.615,.695),'class':(.51,.57)}}
CLASS_COLORS={'LTR':'#CC6677','LINE':'#88CCEE','DNA':'#44AA99','RC':'#DDCC77','Unknown':'#BBBBBB','Mixed':'#555555'}
GAP_DEG = 42.0
SECTOR_GAP_DEG = 1.6
GENE_MAX = 140.0
TE_MAX = 100.0
VARIANTS = {
    'composition': ('A', 'TE-class composition per chromosome',
                    'Solid segments show each class\'s share of TE-covered bases; segment positions are not genomic coordinates.'),
    'dominant': ('B', 'Dominant TE class along each chromosome',
                 'Solid colours show the largest class by covered bases in each 100 kb window; minority classes are not displayed.')}
EMPTY_COLOR = '#FFFFFF'
TIE_COLOR = '#F2E8DE'
PHYLOGENY_URL = 'https://link.springer.com/article/10.1186/1471-2148-7-87'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_tsv(path):
    with path.open(newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def write_tsv(path, rows):
    assert rows
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)


def ingest():
    """Verify the self-contained frozen inputs shipped with this public package."""
    manifest = json.loads((ROOT/'provenance/input_manifest.json').read_text())
    for record in manifest:
        path = ROOT/record['frozen_file']
        assert path.is_file(), f"Missing frozen input: {record['frozen_file']}"
        assert sha(path) == record['sha256'], f"Input checksum mismatch: {record['frozen_file']}"
        assert path.stat().st_size == record['bytes']
    receipt = json.loads((ROOT/'source/export_receipt.json').read_text())
    assert receipt['state'] == 'DENSITIES_VALIDATED'
    return receipt


def assemble(allow_missing_classes=False):
    labels = read_tsv(ROOT/'source/chromosome_labels_and_Muller_elements.tsv')
    groups = defaultdict(list)
    for r in labels:
        groups[r['species'], r['chromosome_or_group']].append(r['accession'])
    all_sectors, all_windows, scaffold_rows, excluded_rows, summaries = {}, [], [], [], []
    for sp in SPECIES:
        seqs = {r['accession']: r for r in read_tsv(ROOT/f'source/{sp}_sequences.tsv')}
        class_path=ROOT/f'tables/{sp}_TE_classes_100kb.tsv'
        class_rows=None
        if class_path.exists():
            checked=json.loads((ROOT/f'provenance/{sp}_class_window_receipt.json').read_text())
            assert checked['state']=='VALIDATED' and sha(class_path)==checked['output_sha256']
            class_rows={(r['accession'],int(r['start0'])):r for r in read_tsv(class_path)}
        else:
            assert allow_missing_classes,f'Validated class windows are required for {sp}'
        windows = defaultdict(list)
        for r in read_tsv(ROOT/f'source/{sp}_windows_100kb.tsv'):
            windows[r['accession']].append(r)
        total_bp = sum(int(seqs[a]['span_bp']) for group in ORDER[sp] for a in groups[sp, group])
        available = 360-GAP_DEG-SECTOR_GAP_DEG*(len(ORDER[sp])-1)
        angle = GAP_DEG/2
        sectors = []
        selected = set()
        for group, label in zip(ORDER[sp], LABELS[sp]):
            accessions = sorted(groups[sp, group], key=lambda a: (-int(seqs[a]['span_bp']), a))
            selected.update(accessions)
            span = sum(int(seqs[a]['span_bp']) for a in accessions)
            theta0, theta1 = np.deg2rad([angle, angle+available*span/total_bp])
            sector = {'species': sp, 'group': group, 'label': label, 'span_bp': span,
                      'theta0': float(theta0), 'theta1': float(theta1), 'windows': [],
                      'boundaries': [], 'accessions': accessions}
            offset = 0
            gene_sum = te_sum = 0
            for rank, acc in enumerate(accessions, 1):
                length = int(seqs[acc]['span_bp'])
                wrows = sorted(windows[acc], key=lambda r: int(r['start0']))
                cursor = 0
                for w in wrows:
                    start, end = int(w['start0']), int(w['end0'])
                    assert start == cursor and end > start
                    assert end-start == int(w['span_bp']) and end <= length
                    assert 0 <= int(w['TE_bp']) <= end-start
                    assert math.isclose(float(w['genes_per_100kb_span']), int(w['gene_count'])*1e5/(end-start))
                    assert math.isclose(float(w['TE_pct_span']), int(w['TE_bp'])*100/(end-start), abs_tol=1e-10)
                    assert float(w['genes_per_100kb_span']) <= GENE_MAX
                    assert 0 <= float(w['TE_pct_span']) <= TE_MAX
                    out = dict(w)
                    out.update(display_sector=group, display_start0=offset+start, display_end0=offset+end,
                               scaffold_display_rank=rank)
                    out['gene_density_0_1'] = float(w['genes_per_100kb_span']) / GENE_MAX
                    out['TE_density_0_1'] = float(w['TE_pct_span']) / TE_MAX
                    for cl in CLASS_COLORS:
                        out[cl+'_bp']=int(class_rows[acc,start][cl+'_bp']) if class_rows is not None else ''
                    out['class_data_available']=int(class_rows is not None)
                    if class_rows is not None:
                        assert sum(out[cl+'_bp'] for cl in CLASS_COLORS)==int(w['TE_bp'])
                        top_bp=max(out[cl+'_bp'] for cl in CLASS_COLORS)
                        leaders=[cl for cl in CLASS_COLORS if out[cl+'_bp']==top_bp]
                        out['dominant_class'] = ('No TE' if top_bp==0 else
                                                 leaders[0] if len(leaders)==1 else 'Tie')
                        out['dominant_class_bp'] = top_bp
                        out['dominant_fraction_of_TE'] = top_bp/int(w['TE_bp']) if int(w['TE_bp']) else ''
                        out['tied_classes'] = ';'.join(leaders) if top_bp>0 and len(leaders)>1 else ''
                    all_windows.append(out)
                    sector['windows'].append(out)
                    cursor = end
                assert cursor == length
                assert sum(int(w['gene_count']) for w in wrows) == int(seqs[acc]['genes'])
                assert sum(int(w['TE_bp']) for w in wrows) == int(seqs[acc]['TE_bp'])
                gene_sum += int(seqs[acc]['genes'])
                te_sum += int(seqs[acc]['TE_bp'])
                if group == 'Y scaffolds':
                    scaffold_rows.append({'species': sp, 'display_sector': group, 'display_rank': rank,
                                          'accession': acc, 'span_bp': length, 'display_start0': offset,
                                          'display_end0': offset+length,
                                          'order_basis': 'descending span, then accession; display only'})
                if rank > 1:
                    sector['boundaries'].append(offset)
                offset += length
            assert offset == span
            summaries.append({'species': sp, 'display_sector': group, 'figure_label': label,
                              'sequence_count': len(accessions), 'span_bp': span,
                              'gene_count': gene_sum, 'TE_union_bp': te_sum,
                              'TE_percent_span': te_sum/span*100,
                              'angle_start_deg': float(np.rad2deg(theta0)),
                              'angle_end_deg': float(np.rad2deg(theta1))})
            angle += available*span/total_bp+SECTOR_GAP_DEG
            sectors.append(sector)
        assert math.isclose(angle-SECTOR_GAP_DEG, 360-GAP_DEG/2, abs_tol=1e-8)
        for acc, s in seqs.items():
            if acc not in selected:
                excluded_rows.append({'species': sp, 'accession': acc, 'span_bp': s['span_bp'],
                                      'reason': 'mitochondrion or unplaced non-Y scaffold; retained in frozen source tables'})
        all_sectors[sp] = sectors
    assert len(scaffold_rows) == 62 and sum(int(r['span_bp']) for r in scaffold_rows) == 110504817
    write_tsv(ROOT/'tables/plotted_windows_100kb.tsv', all_windows)
    write_tsv(ROOT/'tables/Dmir_Y_scaffold_display_offsets.tsv', scaffold_rows)
    write_tsv(ROOT/'tables/plotted_chromosomes.tsv', summaries)
    write_tsv(ROOT/'tables/sequences_outside_circular_display.tsv', excluded_rows)
    composition=[]
    for sp in SPECIES:
        for s in all_sectors[sp]:
            total=sum(int(w['TE_bp']) for w in s['windows'])
            assert total>0, 'Zero-TE chromosome requires an explicit empty composition sector'
            cursor=s['theta0']
            for cl in CLASS_COLORS:
                bp=sum(int(w[cl+'_bp']) for w in s['windows'])
                end=cursor+(s['theta1']-s['theta0'])*bp/total
                composition.append(dict(species=sp,display_sector=s['group'],TE_class=cl,
                    class_bp=bp,total_TE_bp=total,chromosome_span_bp=s['span_bp'],
                    fraction_of_TE=bp/total,fraction_of_chromosome=bp/s['span_bp'],
                    composition_start_deg=float(np.rad2deg(cursor)),composition_end_deg=float(np.rad2deg(end)),
                    position_meaning='composition only; not genomic coordinates'))
                cursor=end
            assert math.isclose(cursor,s['theta1'],abs_tol=1e-10)
    write_tsv(ROOT/'tables/chromosome_TE_class_composition.tsv',composition)
    return all_sectors, summaries, all_windows


def arc(ax, theta0, theta1, radius, **kwargs):
    t = np.linspace(theta0, theta1, max(4, int((theta1-theta0)*150)))
    ax.plot(t, np.full_like(t, radius), **kwargs)


def tangent_rotation(theta):
    deg = -np.rad2deg(theta)
    while deg < -90:
        deg += 180
    while deg > 90:
        deg -= 180
    return deg


def draw_tracks(ax, sp, sectors, variant):
    radii = RADII[sp]
    rchr0, rchr1 = radii['chr']
    for s in sectors:
        t0, t1, span = s['theta0'], s['theta1'], s['span_bp']
        width = t1-t0
        ax.bar((t0+t1)/2, rchr1-rchr0, width=width, bottom=rchr0,
               color=COLORS['chromosome'], linewidth=0, zorder=3)
        for kind, value_key in [('gene', 'gene_density_0_1'), ('TE', 'TE_density_0_1')]:
            lo, hi = radii[kind]
            for level in [0, .5, 1]:
                arc(ax, t0, t1, lo+(hi-lo)*level, color=COLORS['grid'], lw=.45, zorder=0)
            w = s['windows']
            starts = np.array([int(x['display_start0']) for x in w])
            ends = np.array([int(x['display_end0']) for x in w])
            heights = np.array([float(x[value_key]) for x in w]) * (hi-lo)
            ax.bar(t0+(starts+ends)/2/span*width, heights, width=(ends-starts)/span*width,
                   bottom=lo, align='center', color=COLORS[kind], edgecolor='none',
                   linewidth=0, antialiased=False, zorder=2)
        # A separate, constant-thickness ring; both variants use the same data.
        lo,hi=radii['class']
        if not s['windows'][0]['class_data_available']:
            ax.bar((t0+t1)/2,hi-lo,width=width,bottom=lo,facecolor='#FAFAFA',
                   edgecolor='#B0B0B0',hatch='////',linewidth=.35,zorder=2)
        else:
            w=s['windows']
            starts=np.array([int(x['display_start0']) for x in w])
            ends=np.array([int(x['display_end0']) for x in w])
            if variant=='composition':
                total=sum(int(x['TE_bp']) for x in w)
                cursor=t0
                for cl,color in CLASS_COLORS.items():
                    class_width=width*sum(int(x[cl+'_bp']) for x in w)/total
                    if class_width>0:
                        ax.bar(cursor+class_width/2,hi-lo,width=class_width,bottom=lo,
                               color=color,edgecolor='none',linewidth=0,antialiased=False,zorder=2)
                    cursor+=class_width
                assert math.isclose(cursor,t1,abs_tol=1e-10)
            else:
                assert variant=='dominant'
                palette={**CLASS_COLORS,'No TE':EMPTY_COLOR,'Tie':TIE_COLOR}
                ax.bar(t0+(starts+ends)/2/span*width,np.full(len(w),hi-lo),
                       width=(ends-starts)/span*width,bottom=lo,
                       color=[palette[x['dominant_class']] for x in w],
                       edgecolor='none',linewidth=0,antialiased=False,zorder=2)
        arc(ax,t0,t1,lo,color=COLORS['grid'],lw=.35,zorder=3)
        arc(ax,t0,t1,hi,color=COLORS['grid'],lw=.35,zorder=3)
        # Y boundaries are annotations; no scaffold windows are merged or smoothed.
        for boundary in s['boundaries']:
            angle = t0+boundary/span*width
            ax.plot([angle,angle], [rchr0-.006,rchr1+.006], color=COLORS['boundary'],
                    lw=.45, alpha=.8, zorder=4)
        # Coordinates are in Mb, except Y-scaffold sector is cumulative displayed span.
        tick_mb = 20 if sp == 'Dmir' else (20 if sp == 'Dpse' else 10)
        tick_positions = np.arange(0, span, tick_mb*1_000_000)
        if span < 5_000_000:  # tiny sectors: name only, exact length in audit table
            tick_positions = []
        for bp in tick_positions:
            angle = t0+bp/span*width
            ax.plot([angle,angle], [rchr1+.004,rchr1+.014], color=COLORS['boundary'], lw=.5)
            ax.text(angle, rchr1+.031, str(int(bp/1e6)), ha='center', va='center',
                    rotation=tangent_rotation(angle), rotation_mode='anchor', fontsize=7.8,
                    color=COLORS['ink'])
        mid = (t0+t1)/2
        # Fine leaders separate adjacent tiny Dmel 4 and Y sectors.
        label_theta = mid
        label_r = rchr1+.085 if sp == 'Dmir' else rchr1+.060
        if sp == 'Dmel' and s['group'] == '4':
            label_theta -= np.deg2rad(2.8)
            label_r += .016
        if sp == 'Dmel' and s['group'] == 'Y':
            label_theta += np.deg2rad(2)
        if s['group'] in ['4', 'Y'] and sp == 'Dmel':
            ax.plot([mid, label_theta], [rchr1+.014, label_r-.018], color=COLORS['boundary'], lw=.5)
        ax.text(label_theta, label_r, s['label'], ha='center', va='center',
                rotation=tangent_rotation(label_theta), rotation_mode='anchor',
                fontsize=10.2 if sp == 'Dmir' else (9.2 if sp == 'Dpse' else 8.8),
                fontweight='semibold', color=COLORS['ink'])
    # Species names occupy the centre of the shared top gap at each group.
    name_r = (radii['class'][0]+rchr1)/2
    ax.text(0, name_r+.015, NAMES[sp], ha='center', va='center', fontsize=11.5 if sp=='Dmel' else 12.7,
            fontweight='bold', fontstyle='italic', color=COLORS['ink'])
    ax.text(0, name_r-.038, sp, ha='center', va='center', fontsize=9.2, color='#626C72')
    # Small colour keys at the gap end identify track types, not extra data.
    for kind, key in [('chr', 'chromosome'), ('gene','gene'), ('TE','TE')]:
        lo,hi = radii[kind]
        ax.plot([np.deg2rad(19),np.deg2rad(19)], [lo,hi], color=COLORS[key], lw=2)
        if kind in ['gene','TE']:
            for level in [0,.5,1]:
                r=lo+(hi-lo)*level
                ax.plot(np.deg2rad([-19,-18.5]),[r,r],color=COLORS['boundary'],lw=.5)
                ax.text(np.deg2rad(-18.8),r,f'{level:g}',ha='right',va='center',
                        fontsize=6.2,color='#58636A')
    lo,hi=radii['class']
    for i,(cl,color) in enumerate(CLASS_COLORS.items()):
        ax.plot([np.deg2rad(19),np.deg2rad(19)],
                [lo+(hi-lo)*i/6,lo+(hi-lo)*(i+1)/6],color=color,lw=2)


def draw_tree(fig, bounds):
    """Ultrametric chronogram scaled to explicit published-model display ages."""
    t = fig.add_axes(bounds)
    ages=json.loads((ROOT/'source/divergence_ages.json').read_text())['nodes']
    root=ages['Dmel_Dpse_Dmir']['plotted_age_Ma']; recent=ages['Dpse_Dmir']['plotted_age_Ma']
    assert root>recent>0
    t.set(xlim=(38,-12),ylim=(-.65,3.4)); t.axis('off')
    t.text(13,3.2,'DIVERGENCE TIME',ha='center',va='center',fontsize=10,fontweight='bold',color=COLORS['ink'])
    for xs,ys in [([root,root],[.65,1.925]),([root,0],[.65,.65]),
                  ([root,recent],[1.925,1.925]),([recent,recent],[1.55,2.3]),
                  ([recent,0],[2.3,2.3]),([recent,0],[1.55,1.55])]:
        t.plot(xs,ys,color=COLORS['ink'],lw=1.45,solid_capstyle='round')
    for label,y in [('Dmir',2.3),('Dpse',1.55),('Dmel',.65)]:
        t.text(-1.4,y,label,ha='left',va='center',fontsize=10.2,fontstyle='italic',color=COLORS['ink'])
    for name,y in [('Dmel_Dpse_Dmir',1.2875),('Dpse_Dmir',1.925)]:
        lo,hi=ages[name]['reported_range_Ma']
        t.plot([hi,lo],[y,y],color='#999999',lw=5,alpha=.5,solid_capstyle='butt')
        t.scatter([ages[name]['plotted_age_Ma']],[y],s=11,color=COLORS['ink'],zorder=5)
    t.text(root,2.2,'32.5 Ma',ha='center',va='bottom',fontsize=8.4)
    t.annotate('4.35 Ma',xy=(recent,1.925),xytext=(14,2.7),ha='center',fontsize=8.4,
               arrowprops={'arrowstyle':'-','color':'#999999','lw':.6})
    t.plot([35,0],[.1,.1],color='#777777',lw=.65)
    for age in [30,20,10,0]:
        t.plot([age,age],[.1,.02],color='#777777',lw=.65)
        t.text(age,-.08,str(age),ha='center',va='top',fontsize=7.6)
    t.text(17.5,-.45,'Million years ago',ha='center',va='center',fontsize=8.1,color='#626C72')


def render(sectors, variant):
    matplotlib.rcParams.update({'font.family':'DejaVu Sans','pdf.fonttype':42,
                                'ps.fonttype':42,'svg.fonttype':'none'})
    fig = plt.figure(figsize=(18,19.5), facecolor='white')
    fig.text(.5,.975,'Gene and TE density across three Drosophila genomes',
             ha='center',va='top',fontsize=22,fontweight='bold',color=COLORS['ink'])
    fig.text(.5,.946,'Dmir → Dpse → Dmel, outside to inside · four tracks per species',
             ha='center',va='top',fontsize=12.5,color='#626C72')
    letter,title,explanation=VARIANTS[variant]
    fig.text(.5,.925,f'{letter}  |  {title}',ha='center',va='top',fontsize=13,
             fontweight='semibold',color=COLORS['ink'])
    ax = fig.add_axes([.035,.14,.93,.79], projection='polar')
    ax.set_theta_zero_location('N'); ax.set_theta_direction(-1)
    ax.set_ylim(0,2.03); ax.set_axis_off()
    for sp in SPECIES:
        draw_tracks(ax, sp, sectors[sp], variant)
    # Compute the tree inset from the actual square polar axes position.
    fig.canvas.draw()
    bbox = ax.get_position()
    cx,cy = bbox.x0+bbox.width/2,bbox.y0+bbox.height/2
    tree_width = bbox.width*.235
    tree_height = bbox.height*.230
    draw_tree(fig, [cx-tree_width/2,cy-tree_height/2,tree_width,tree_height])
    handles = [Patch(facecolor=COLORS['chromosome'],label='Chromosome span (Mb)'),
               Patch(facecolor=COLORS['gene'],label='Gene density: 0–1'),
               Patch(facecolor=COLORS['TE'],label='TE density: 0–1')]
    first=fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.113),ncol=3,
               frameon=False,fontsize=11,handlelength=1.5,columnspacing=2.3)
    class_handles=[Patch(facecolor=color,label=cl if cl!='RC' else 'RC / Helitron') for cl,color in CLASS_COLORS.items()]
    if variant=='dominant':
        class_handles.append(Patch(facecolor=EMPTY_COLOR,edgecolor='#BBBBBB',label='No TE'))
        if any(w['dominant_class']=='Tie' for sp in SPECIES for s in sectors[sp] for w in s['windows']):
            class_handles.append(Patch(facecolor=TIE_COLOR,edgecolor='#BBBBBB',label='Equal top classes'))
    fig.legend(handles=class_handles,loc='lower center',bbox_to_anchor=(.5,.09),ncol=len(class_handles),
               frameon=False,fontsize=10,handlelength=1.5,columnspacing=1.3)
    notes = [
        'Common 0–1 scales: gene density = genes per 100 kb / 140; TE density = TE-covered bases / assembled window span.',
        explanation+' Mixed = overlapping classes.',
        '† Dmir Y: all 62 Y-associated scaffolds (110.50 Mb); length-ordered for display only. Ticks show cumulative span; small marks show boundaries.',
        'Each species spans the same arc independently: angular alignment does not imply synteny. Non-Y unplaced scaffolds and mitochondria are not shown.',
        '* Dmir XR includes A-derived material; Dpse X combines A + D. Gene/TE densities use assembled span; terminal windows use their actual length.',
        'Dated tree: Gao et al. (2007) model; plotted range midpoints 32.5 and 4.35 Ma, ranges 30–35 and 3.76–4.93 Ma (grey bars). Ages are model-dependent.'
    ]
    for i,note in enumerate(notes):
        fig.text(.5,.084-i*.0115,note,ha='center',va='top',fontsize=8.3,color='#58636A')
    out = ROOT/f'figures/{letter}_{variant}_circos'
    for ext in ['png','pdf','svg']:
        fig.savefig(out.with_suffix('.'+ext), dpi=220, facecolor='white')
    # Compact preview for inspection and inline delivery, full-resolution PNG above.
    fig.savefig(ROOT/f'figures/{letter}_{variant}_preview.png',dpi=110,facecolor='white')
    plt.close(fig)


def main():
    p=argparse.ArgumentParser()
    args=p.parse_args()
    receipt=ingest()
    sectors,summaries,windows=assemble()
    for variant in VARIANTS:
        render(sectors,variant)
    y=next(r for r in summaries if r['species']=='Dmir' and r['display_sector']=='Y scaffolds')
    stats = {'state':'NUMERICALLY_VALIDATED_RENDERED_AWAITING_VISUAL_REVIEW',
             'created_utc':datetime.now(timezone.utc).isoformat(),
             'source_state':receipt['state'],'tracks_per_species':4,'total_data_tracks':12,
             'ring_order_outer_to_inner':SPECIES,'data_window_count':len(windows),
             'max_gene_density':max(float(w['genes_per_100kb_span']) for w in windows),
             'common_gene_scale':[0,1],'common_TE_scale':[0,1],
             'gene_normalization':'genes_per_100kb_span / 140; same constant for all species',
             'TE_normalization':'TE_bp / span_bp',
             'variants':{k:{'label':v[0],'title':v[1],'description':v[2]} for k,v in VARIANTS.items()},
             'dominant_tie_window_count':sum(w['dominant_class']=='Tie' for w in windows),
             'no_TE_window_count':sum(w['dominant_class']=='No TE' for w in windows),
             'Dmir_Y':y,'scaffold_order':'descending assembled span; display only',
             'genome_based_not_expression':True,'new_mapping_or_annotation_performed':False,
             'phylogeny_newick':'(Dmel:32.5,(Dpse:4.345,Dmir:4.345):28.155);','phylogeny_source':PHYLOGENY_URL,
             'phylogeny_is_time_scaled':True,'phylogeny_ages':json.loads((ROOT/'source/divergence_ages.json').read_text()),
             'missing_class_species':[sp for sp in SPECIES if not (ROOT/f'tables/{sp}_TE_classes_100kb.tsv').exists()],
             'matplotlib':matplotlib.__version__,'numpy':np.__version__}
    if stats['missing_class_species']:
        stats['state']='INCOMPLETE_CLASS_DATA_PREVIEW'
    (ROOT/'provenance/figure_validation.json').write_text(json.dumps(stats,indent=2)+'\n')
    print(json.dumps(stats,indent=2))


if __name__=='__main__':
    main()
