#!/usr/bin/env python3
"""Collapse rare annotated target-ALT alleles into gene-level carrier indicators."""
import argparse,json,sys
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'06_Variant_GWAS'/'helpers'))
from master_inputs import load_master,phenotype_cohort,write_bed,DRUGS

def collapse(d):
    present=(d>0).any(axis=0)
    missing=(d<0).any(axis=0)&~present
    c=present.astype(np.int8);c[missing]=-1
    return c

def build(master,catalog,gene_map,phenotype_dir,outdir,drugs=DRUGS,cap=.05):
    metadata,samples,x=load_master(master);outdir=Path(outdir);outdir.mkdir(parents=True,exist_ok=True)
    ann=pd.read_csv(catalog,sep='\t').fillna('');genes=pd.read_csv(gene_map,sep='\t').fillna('').set_index('gene_id')
    if ann.variant_idx.duplicated().any() or not ann.impact.isin(['HIGH','MODERATE']).all() or not ann.gene_id.isin(genes.index).all():
        raise ValueError('Selected annotations must contain one mapped HIGH/MODERATE consequence per allele')
    if metadata.iloc[ann.variant_idx].id.tolist()!=ann.variant_id.tolist():
        raise ValueError('Annotation alleles do not match master row indices')
    manifest=[]
    for drug in drugs:
        sample_idx,ordered,p=phenotype_cohort(phenotype_dir,drug,samples)
        dest=outdir/drug;dest.mkdir(exist_ok=True)
        carrier=np.empty(len(ann),np.int64)
        for lo in range(0,len(ann),4096):
            d=x[ann.variant_idx.iloc[lo:lo+4096].to_numpy(int)][:,sample_idx]
            carrier[lo:lo+len(d)]=(d>0).sum(axis=1)
        membership=ann.copy();membership['n_carriers']=carrier
        membership['carrier_frequency_total_cohort']=carrier/len(ordered)
        membership=membership.loc[(carrier>0)&(carrier/len(ordered)<=cap)]
        rows=[];codes=[]
        for gid,part in membership.groupby('gene_id',sort=True):
            code=collapse(x[part.variant_idx.to_numpy(int)][:,sample_idx]);g=genes.loc[gid]
            nc=int((code==1).sum());nn=int((code==0).sum());nm=int((code<0).sum())
            rows.append({'gene_id':gid,'gene_symbol':g.gene_symbol,'id':'GB:'+gid,'gene_marker_id':'GB:'+gid,
                         'chromosome':int(g.chromosome),'pos':(int(g.start)+int(g.end))//2,'bp':(int(g.start)+int(g.end))//2,
                         'gene_start':int(g.start),'gene_end':int(g.end),'n_collapsed_variants':len(part),
                         'n_carriers':nc,'n_noncarriers':nn,'n_missing':nm,'n_called':nc+nn,
                         'carrier_frequency_called':nc/(nc+nn) if nc+nn else np.nan,
                         'is_polymorphic_gene_marker':nc>0 and nn>0,'passes_gene_marker_maf_filter':nc>0 and nn>0})
            codes.append(code)
        qc=pd.DataFrame(rows);all_codes=np.vstack(codes);keep=qc.is_polymorphic_gene_marker.to_numpy(bool)
        tested=qc.loc[keep].copy();c=all_codes[keep];np.save(dest/'codes.npy',c)
        tested.to_csv(dest/'gene_metadata.tsv',sep='\t',index=False)
        pd.DataFrame({'sample':ordered}).to_csv(dest/'samples.tsv',sep='\t',index=False)
        tested['cm']=0;tested['alt']='2';tested['ref']='1'
        dose=np.where(c<0,-1,2*c).astype(np.int8);prefix=dest/'genotypes'
        write_bed(prefix,tested,ordered,[dose])
        p[['FID','sample','log2_mic']].to_csv(dest/'phenotype.log2_mic.txt',sep='\t',index=False,header=False)
        qc.to_csv(dest/'gene_qc.tsv',sep='\t',index=False)
        membership['gene_marker_tested']=membership.gene_id.isin(tested.gene_id)
        membership.to_csv(dest/'variant_membership.tsv.gz',sep='\t',index=False)
        row={'drug':drug,'n':len(ordered),'tested_gene_markers':len(tested),'mitochondrial_gene_markers':int((tested.chromosome==8).sum()),
             'qualifying_variants':len(membership),'genotype_prefix':str(prefix.resolve()),'phenotype_path':str((dest/'phenotype.log2_mic.txt').resolve()),
             'qc_path':str((dest/'gene_qc.tsv').resolve()),'maximum_carrier_frequency':cap}
        manifest.append(row);print(json.dumps(row),flush=True)
    pd.DataFrame(manifest).to_csv(outdir/'prepared_scan_manifest.tsv',sep='\t',index=False)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--master-prefix',required=True);p.add_argument('--catalog',required=True);p.add_argument('--gene-map',required=True)
    p.add_argument('--phenotype-dir',required=True);p.add_argument('--outdir',required=True);p.add_argument('--drug',action='append',choices=DRUGS)
    a=p.parse_args();build(a.master_prefix,a.catalog,a.gene_map,a.phenotype_dir,a.outdir,a.drug or DRUGS)
