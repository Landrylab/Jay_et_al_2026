#!/usr/bin/env python3
"""Render the complete within-background importance figure with an editable network."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import tempfile
os.environ.setdefault('MPLCONFIGDIR', str(Path(tempfile.gettempdir())/'mpl-within-groups'))
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.font_manager import FontProperties
from matplotlib.path import Path as MplPath
from matplotlib.patches import PathPatch
from matplotlib.textpath import TextPath
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd
from PIL import Image

OUT = SOURCE = IMPORTANCE_PATH = NETWORK_PATH = None
GROUPS = []
GROUP_COLOURS = {}
LEVEL_COLOURS = {'direct': '#8E145F', 'response': '#D55E00', 'pathway': '#000000', 'none': '#8B9298'}
WIDTH, HEIGHT = 180, 225
FONT_SIZE = 5.1


def style():
    for name in ['Arial.ttf', 'Arial Bold.ttf', 'Arial Italic.ttf', 'Arial Bold Italic.ttf']:
        path = Path('/System/Library/Fonts/Supplemental') / name
        if path.exists(): font_manager.fontManager.addfont(str(path))
    mpl.rcParams.update({'font.family': ['Arial', 'Liberation Sans', 'DejaVu Sans'], 'font.size': 5.5,
        'axes.labelsize': 6, 'axes.labelweight': 'bold', 'axes.linewidth': .6,
        'xtick.labelsize': 5.3, 'ytick.labelsize': 5.3, 'xtick.major.size': 2,
        'ytick.major.size': 2, 'xtick.major.width': .5, 'ytick.major.width': .5,
        'text.color': 'black', 'axes.labelcolor': 'black', 'xtick.color': 'black',
        'ytick.color': 'black', 'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none', 'text.hinting': 'none'})


def physical_axes(fig, x, y, w, h):
    return fig.add_axes([x/WIDTH, y/HEIGHT, w/WIDTH, h/HEIGHT])


def wrap_term(term, width_mm=53.5):
    if not term or pd.isna(term):
        return ''
    term = str(term).translate(str.maketrans({'₂':'2', '₆':'6'}))
    prop = FontProperties(family=['Arial', 'Liberation Sans', 'DejaVu Sans'], size=FONT_SIZE)
    max_width = width_mm / 25.4 * 72
    lines, line = [], ''
    for word in str(term).split():
        candidate = f'{line} {word}'.strip()
        if line and TextPath((0, 0), candidate, prop=prop).get_extents().width > max_width:
            lines.append(line)
            line = word
        else:
            line = candidate
    lines.append(line)
    return '\n'.join(lines)


LABEL_OFFSETS = {}

def evidence_labels(data, evidence, group_key, cluster):
    selected = evidence.loc[evidence.cluster.astype(str).eq(str(cluster)) &
                            evidence.evidence_level.isin(['direct', 'response'])]
    if selected.gene_id.duplicated().any():
        raise ValueError('A gene may have only one evidence category per group')
    return data.loc[data['rank'].le(20)].merge(
        selected[['gene_id', 'evidence_level']], on='gene_id', validate='one_to_one')


def draw_importance(fig, ax, importance, network, group_info, letter):
    cluster, group, n, r = group_info
    data = importance.loc[importance.group.eq(group)].copy()
    chroms = [str(c) for c in range(1, 8)] + ['R', 'MT']
    centres = []
    plotted = 0
    for i, chrom in enumerate(chroms):
        frame = data.loc[data.chrom_plot.astype(str).eq(chrom)]
        shared = importance.loc[importance.chrom_plot.astype(str).eq(chrom)]
        centres.append((shared.genome_position.min() + shared.genome_position.max()) / 2)
        plotted += len(frame)
        ax.scatter(frame.genome_position, frame.mean_importance_across_folds,
                   s=3.8, color=['#0072B2', '#999999'][i % 2], alpha=.8, linewidths=0)
    assert plotted == len(data)
    labels = evidence_labels(data, network, group, cluster)
    for row in labels.itertuples():
        color = LEVEL_COLOURS[row.evidence_level]
        ax.scatter([row.genome_position], [row.mean_importance_across_folds], s=13,
                   color=color, edgecolors='#222222', linewidths=.35, zorder=4)
        dx, dy, align = LABEL_OFFSETS.get(group, {}).get(row.gene_symbol, (0, 5, 'center'))
        ax.annotate(row.gene_symbol, (row.genome_position, row.mean_importance_across_folds),
                    xytext=(dx, dy), textcoords='offset points', ha=align,
                    va='bottom' if dy > 0 else 'top', fontsize=5.1, fontstyle='italic',
                    fontweight='bold', color=color, zorder=6,
                    arrowprops={'arrowstyle': '-', 'lw': .3, 'color': color,
                                'shrinkA': 1.0, 'shrinkB': 2.3} if dx or abs(dy)>7 else None)
    ax.set_xticks(centres, chroms)
    ax.set_xlabel('Chromosome', labelpad=2)
    ax.set_ylabel('Mean Gini importance', labelpad=3)
    ax.set_ylim(0, data.mean_importance_across_folds.max()*1.29)
    extent = importance.genome_position.max()
    ax.set_xlim(-extent*.018, extent*1.015)
    ax.yaxis.set_major_locator(MaxNLocator(4, min_n_ticks=3))
    ax.spines[['top', 'right']].set_visible(False)
    for tick in [*ax.get_xticklabels(), *ax.get_yticklabels()]:
        tick.set_fontweight('bold')
    ax.tick_params(pad=1.5)
    pos = ax.get_position()
    fig.text(pos.x0-.053, pos.y1+.010, letter, fontsize=8, fontweight='bold')
    fig.text((pos.x0+pos.x1)/2, pos.y1+.008, group,
             fontsize=6.5, fontweight='bold', ha='center', va='bottom')
    return labels[['group', 'gene_id', 'gene_symbol', 'genome_position',
                   'mean_importance_across_folds', 'evidence_level']]


def curved_link(ax, origin, target, color):
    x0, y0 = origin
    x1, y1 = target
    direction = np.sign(x1-x0)
    path = MplPath([(x0, y0), (x0+direction*4, y0), (x1-direction*4, y1), (x1, y1)],
                   [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4])
    ax.add_patch(PathPatch(path, fill=False, lw=.45, color=color, alpha=.62, zorder=1))


def draw_network(fig, network):
    unique = network.drop_duplicates('gene_symbol').set_index('gene_symbol')
    ordered = []
    for cluster, _, _, _ in GROUPS:
        for gene in network.loc[network.cluster.eq(cluster)].sort_values('rank').gene_symbol:
            if gene not in ordered: ordered.append(gene)
    left_groups = {group[0] for group in GROUPS[:2]}
    left_members = set(network.loc[network.cluster.isin(left_groups), 'gene_symbol'])
    banks = [[gene for gene in ordered if gene in left_members],
             [gene for gene in ordered if gene not in left_members]]
    # Shared genes appear once, with all links, in the first matching group bank.
    group_order = {group[0]: index for index, group in enumerate(GROUPS)}
    memberships = {gene: rows for gene, rows in network.groupby('gene_symbol')}
    def display_order(gene):
        rows = memberships[gene]
        first_group = min(group_order[cluster] for cluster in rows.cluster)
        return (-int(unique.loc[gene, 'n_cluster_rf_links']), first_group,
                int(rows['rank'].min()), str(gene))
    for bank in banks:
        bank.sort(key=display_order)
        single_links = [gene for gene in bank if int(unique.loc[gene, 'n_cluster_rf_links']) == 1]
        single_group_order = [group_order[int(memberships[gene].iloc[0].cluster)]
                              for gene in single_links]
        assert single_group_order == sorted(single_group_order)
    left, right = banks
    for bank in banks:
        assert all(unique.loc[a, 'n_cluster_rf_links'] >= unique.loc[b, 'n_cluster_rf_links']
                   for a, b in zip(bank, bank[1:]))
    assert set(left).isdisjoint(right) and len(left) + len(right) == network.gene_symbol.nunique()
    canvas = physical_axes(fig, 0, 0, WIDTH, HEIGHT)
    canvas.set(xlim=(0, WIDTH), ylim=(0, HEIGHT)); canvas.set_axis_off()
    canvas.text(3, 127.4, 'e', fontsize=8, fontweight='bold', va='bottom')
    terms = {gene: wrap_term(row.evidence_term) for gene, row in unique.iterrows()}
    coordinates, description_rows = {}, []
    for side, genes, gene_x, term_x, node_x in [('left',left,3,18,75), ('right',right,108,123,105)]:
        heights = np.array([2.3 + 1.9*terms[g].count('\n') for g in genes])
        assert heights.sum() <= 114, (side, heights.sum())
        heights += (114 - heights.sum()) / len(heights)
        top = 123
        for gene, height in zip(genes, heights):
            y = top-height/2; top -= height
            row = unique.loc[gene]
            coordinates[gene] = (node_x, y)
            color = LEVEL_COLOURS[row.evidence_level]
            weight = 'bold' if row.evidence_level in ['direct','response'] or row.n_cluster_rf_links>1 else 'normal'
            canvas.text(gene_x, y, gene, fontsize=FONT_SIZE,
                        fontstyle='italic' if row.evidence_level!='none' else 'normal',
                        fontweight=weight, color=color, ha='left', va='center', zorder=5)
            if terms[gene]:
                canvas.text(term_x, y, terms[gene], fontsize=FONT_SIZE, color=color,
                            ha='left', va='center', linespacing=1.12, zorder=5)
            canvas.scatter([node_x], [y], s=9*row.n_cluster_rf_links, color=color,
                           edgecolors='#111111' if row.n_cluster_rf_links>1 else 'white',
                           linewidths=.5 if row.n_cluster_rf_links>1 else .25, zorder=4)
            description_rows.append({**row.to_dict(), 'gene_symbol':gene, 'bank':side,
                                     'display_term':terms[gene], 'x_mm':node_x, 'y_mm':y})
        canvas.text(gene_x,126,'Gene',fontsize=5.4,fontweight='bold')
        canvas.text(term_x,126,'Functional or response evidence',fontsize=5.4,fontweight='bold')
    hubs = dict(zip([row[0] for row in GROUPS], [(85,96),(85,37),(95,101),(95,42)]))
    edge_rows = []
    for row in network.itertuples():
        curved_link(canvas,hubs[row.cluster],coordinates[row.gene_symbol],GROUP_COLOURS[row.cluster])
        edge_rows.append({'cluster':row.cluster,'gene_symbol':row.gene_symbol,
                         'rank':row.rank,'group_x_mm':hubs[row.cluster][0],
                         'group_y_mm':hubs[row.cluster][1],
                         'gene_x_mm':coordinates[row.gene_symbol][0],
                         'gene_y_mm':coordinates[row.gene_symbol][1]})
    for cluster,label,_,_ in GROUPS:
        x,y=hubs[cluster]
        canvas.scatter([x],[y],s=34,color=GROUP_COLOURS[cluster],edgecolors='white',linewidths=.45,zorder=6)
        text = label.replace('PCA cluster', 'PCA\ncluster')
        canvas.text(x,y+2.5,text,ha='center',va='bottom',fontsize=5.5,fontweight='bold',
                    color='black',linespacing=1.0,bbox=dict(facecolor='white',edgecolor='none',pad=.35),zorder=7)
    key=[(3,'direct','Azole response / target'),(50,'response','Other drug / wall evidence'),
         (111,'pathway','Other function'),(148,'none','Uninformative')]
    for x,level,label in key:
        canvas.scatter([x],[3.5],s=9,color=LEVEL_COLOURS[level],linewidths=0)
        canvas.text(x+2,3.5,label,fontsize=5.0,va='center')
    return pd.DataFrame(description_rows),pd.DataFrame(edge_rows)


def main():
    global OUT, SOURCE, IMPORTANCE_PATH, NETWORK_PATH, GROUPS, GROUP_COLOURS, LABEL_OFFSETS
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--importance-data',type=Path,required=True)
    parser.add_argument('--network-data',type=Path,required=True)
    parser.add_argument('--outdir',type=Path,required=True)
    parser.add_argument('--label-offsets',type=Path,help='Optional JSON mapping group and gene to x/y offsets and alignment')
    args=parser.parse_args()
    OUT=args.outdir; SOURCE=OUT/'source_data'/'plot_data'
    OUT.mkdir(parents=True,exist_ok=True); SOURCE.mkdir(parents=True,exist_ok=True)
    (OUT/'provenance').mkdir(exist_ok=True)
    IMPORTANCE_PATH=args.importance_data; NETWORK_PATH=args.network_data
    if args.label_offsets: LABEL_OFFSETS=json.loads(args.label_offsets.read_text())
    style()
    importance=pd.read_csv(IMPORTANCE_PATH,sep='\t')
    network=pd.read_csv(NETWORK_PATH,sep='\t')
    groups=network[['cluster','group']].drop_duplicates()
    assert len(groups)==4 and len(network)==80
    GROUPS=[(row.cluster,row.group,None,None) for row in groups.itertuples()]
    GROUP_COLOURS=dict(zip(groups.cluster,['#CE78A7','#08A985','#1683C5','#E8A000']))
    assert importance.groupby('group').size().min()>0
    assert network.groupby('cluster').size().eq(20).all()
    fig=plt.figure(figsize=(WIDTH/25.4,HEIGHT/25.4))
    positions=[(15,183,68,31),(105,183,68,31),(15,139,68,31),(105,139,68,31)]
    labels=[]
    for letter,position,group in zip('abcd',positions,GROUPS):
        labels.append(draw_importance(fig,physical_axes(fig,*position),importance,network,group,letter))
    description,edges=draw_network(fig,network)
    for ext in ['pdf','png','svg']:
        fig.savefig(OUT/f'figure.{ext}',dpi=600,facecolor='white')
    with Image.open(OUT/'figure.png') as im:
        im.convert('RGB').save(OUT/'figure.png',dpi=(600,600))
    plt.close(fig)
    import shutil
    shutil.copyfile(OUT/'figure.pdf', OUT/'figure_nature_size.pdf')
    description.to_csv(SOURCE/'network_gene_descriptions.tsv',sep='\t',index=False)
    edges.to_csv(SOURCE/'network_memberships.tsv',sep='\t',index=False)
    pd.concat(labels).to_csv(SOURCE/'highlighted_gene_positions.tsv',sep='\t',index=False)
    importance.to_csv(SOURCE/'importance_plot_data.tsv',sep='\t',index=False)
    report={'width_mm':WIDTH,'height_mm':HEIGHT,'importance_rows':len(importance),
            'highlighted_gene_group_points':len(pd.concat(labels)),'unique_network_genes':network.gene_symbol.nunique(),'network_edges':80,
            'smallest_font_pt':5.0,'network_font_pt':FONT_SIZE,
            'network_banks':{'left':[group[1] for group in GROUPS[:2]],
                             'right':[group[1] for group in GROUPS[2:]]},
            'network_order':'Decreasing link count within each bank, then group order, within-group rank and gene symbol; single-link genes follow clade 1 then clade 3 on the left and clade 4 then PCA cluster 6 on the right. Shared genes appear once with all links.',
            'wrapped_terms':int(description.display_term.fillna('').str.contains('\n').sum()),
            'shortened_terms':0,'input_files':[{ 'path':str(p.resolve()),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
             for p in [IMPORTANCE_PATH,NETWORK_PATH]]}
    (OUT/'provenance/source_audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__': main()
