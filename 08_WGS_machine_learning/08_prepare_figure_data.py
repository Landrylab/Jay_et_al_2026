#!/usr/bin/env python3
"""Prepare within-background gene ranks, genomic coordinates and evidence labels."""
from pathlib import Path
import argparse
import json
import re
import shutil
import numpy as np
import pandas as pd



def short_description(description):
    description=re.sub(r'^\s*\((?:orf|Caalf)[^)]*\)\s*','',str(description or '')).strip()
    if not description or re.match(r'^(Protein of unknown function|Uncharacterized|Hypothetical protein)',description,re.I):
        return 'none','Unknown function'
    return 'pathway',description.split(';')[0].strip().rstrip('.')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rf-root',type=Path,required=True)
    parser.add_argument('--groups',type=Path,required=True)
    parser.add_argument('--gene-map',type=Path,required=True)
    parser.add_argument('--reference-fai',type=Path,required=True)
    parser.add_argument('--evidence',type=Path,required=True)
    parser.add_argument('--outdir',type=Path,required=True)
    args=parser.parse_args()
    genes=pd.read_csv(args.gene_map,sep='\t').fillna('').set_index('gene_id')
    evidence=pd.read_csv(args.evidence,sep='\t').fillna('').set_index('gene_id')
    reference=pd.read_csv(args.reference_fai,sep='\t',header=None,names=['ref','length','offset','bases','bytes'])
    lengths={re.search(r'chr([1-7RM])',r.ref).group(1):int(r.length) for r in reference.itertuples()}
    lengths['MT']=lengths.pop('M')
    offsets={};total=0
    for chrom in ['1','2','3','4','5','6','7','R','MT']:
        offsets[chrom]=total;total+=lengths[chrom]+250000
    out=args.outdir/'source_data';plot=out/'plot_data';plot.mkdir(parents=True,exist_ok=True)
    tables,tops=[] , []
    group_table=pd.read_csv(args.groups,sep='\t')
    for group in group_table.itertuples():
        key,label,cluster=group.group_key,group.group_label,group.cluster
        source=args.rf_root/'within_group'/key/'gene'
        info=json.loads((source/'run_info.json').read_text())
        table=pd.read_csv(source/'feature_importance_summary.tsv',sep='\t').fillna('')
        assert np.isclose(table.mean_importance_across_folds.sum(),1)
        table['gene_symbol']=[str(genes.loc[g,'gene_symbol'] or g).removesuffix('_A') for g in table.gene_id]
        table['reference_description']=table.gene_id.map(genes.description)
        table['group']=label;table['cluster']=cluster;table['method_key']='tuned'
        table['chrom_plot']=table.chrom.astype(str).replace({'9':'R','8':'MT','M':'MT','Mt':'MT'})
        assert set(table.chrom_plot).issubset(offsets)
        table['genome_position']=table.chrom_plot.map(offsets)+table.pos
        tables.append(table)
        top=table.loc[table['rank']<=20].copy()
        assert len(top)==20
        levels,terms,sources,source_urls,contexts=[] , [] , [] , [], []
        for row in top.itertuples():
            if row.gene_id in evidence.index:
                e=evidence.loc[row.gene_id];level,term=e.evidence_level,e.evidence_term
                sources.append('Reviewed gene-response evidence catalogue')
                source_urls.append(e.get('evidence_sources',''))
                contexts.append(e.get('evidence_context',''))
            else:
                level,term=short_description(row.reference_description)
                sources.append('CGD reference gene annotation')
                source_urls.append('')
                contexts.append('')
            levels.append(level);terms.append(term)
        top['evidence_level']=levels;top['evidence_term']=terms;top['evidence_source']=sources
        top['evidence_source_urls']=source_urls
        top['evidence_context']=contexts
        top['mean_importance_all_folds']=top.mean_importance_across_folds
        top['mean_rf_gini_importance']=top.mean_importance_across_folds
        tops.append(top)
        target=out/'upstream'/key;target.mkdir(parents=True,exist_ok=True)
        for file in ['run_info.json','feature_importance_summary.tsv']:
            shutil.copyfile(source/file,target/file)
    table=pd.concat(tables,ignore_index=True);top=pd.concat(tops,ignore_index=True)
    memberships=top.groupby('gene_id').cluster.agg(lambda values:sorted(set(values)))
    top['n_cluster_rf_links']=top.gene_id.map(memberships.map(len))
    top['cluster_rf_memberships']=top.gene_id.map(memberships.map(lambda values:','.join(map(str,values))))
    assert top.groupby('gene_id').gene_symbol.nunique().max()==1
    assert top.gene_symbol.nunique()==top.gene_id.nunique()
    table.to_csv(plot/'ML_within_clade_gene_feature_importance_multipanel_plot_data.tsv',sep='\t',index=False)
    top.to_csv(plot/'Within_cluster_gene_burden_top20_evidence_network_associated_genes.tsv',sep='\t',index=False)
    print(json.dumps({'importance_rows':len(table),'network_genes':top.gene_id.nunique(),'connections':len(top),
                      'new_reference_only_genes':sorted(top.loc[top.evidence_source.eq('CGD reference gene annotation'),'gene_symbol'].unique())},indent=2))


if __name__=='__main__':main()
