#!/usr/bin/env python3
"""Filter target-ALT dosages within each quantitative antifungal MIC cohort."""
import argparse,json,sys
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parent/'helpers'))
from master_inputs import load_master,counts,write_bed,phenotype_cohort,DRUGS

def eligible(q):
    return (q['maf']>=.005)&(q['missing_fraction']<=.10)&(q['minor_allele_carriers']>=5)

def prepare(master,phenotype_dir,outdir,drugs=DRUGS,chunk_size=4096):
    metadata,samples,x=load_master(master);outdir=Path(outdir);outdir.mkdir(parents=True,exist_ok=True)
    manifest=[]
    for drug in drugs:
        sample_idx,ordered,p=phenotype_cohort(phenotype_dir,drug,samples)
        dest=outdir/drug;dest.mkdir(exist_ok=True)
        pieces=[]
        for lo in range(0,len(metadata),chunk_size):
            d=x[lo:lo+chunk_size][:,sample_idx];q=counts(d);mask=eligible(q)
            part=metadata.iloc[lo:lo+chunk_size].loc[mask].copy()
            for key,val in q.items():part[key]=val[mask]
            pieces.append(part)
        selected=pd.concat(pieces,ignore_index=True);selected['cm']=0
        selected['SNP']=selected.id;selected['Chr']=selected.chromosome;selected['bp']=selected.pos
        selected['A1']=selected.alt;selected['A2']=selected.ref;selected['marker_type']=selected.type
        selected['effect_estimable']=selected.n_distinct_called_dosages>1
        # A constant dosage has no estimable additive contrast, even when allele MAF is positive.
        selected.to_csv(dest/'eligible_variants.tsv.gz',sep='\t',index=False)
        selected.loc[~selected.effect_estimable].to_csv(dest/'constant_dosage_variants.tsv.gz',sep='\t',index=False)
        tested=selected.loc[selected.effect_estimable].copy()
        ix=tested['index'].to_numpy(int);prefix=dest/'genotypes'
        write_bed(prefix,tested,ordered,(x[ix[i:i+chunk_size]][:,sample_idx] for i in range(0,len(ix),chunk_size)))
        p[['FID','sample','log2_mic']].to_csv(dest/'phenotype.log2_mic.txt',sep='\t',index=False,header=False)
        tested.to_csv(dest/'variant_qc.tsv.gz',sep='\t',index=False)
        record={'drug':drug,'n':len(ordered),'eligible_variants':len(selected),'constant_dosage_variants':int((~selected.effect_estimable).sum()),
                'tested_variants':len(tested),'tested_SNP':int((tested.type=='SNP').sum()),'tested_INDEL':int((tested.type=='INDEL').sum()),
                'tested_mitochondrial':int((tested.chromosome==8).sum()),'genotype_prefix':str(prefix.resolve()),
                'phenotype_path':str((dest/'phenotype.log2_mic.txt').resolve()),'qc_path':str((dest/'variant_qc.tsv.gz').resolve())}
        manifest.append(record);print(json.dumps(record),flush=True)
    pd.DataFrame(manifest).to_csv(outdir/'prepared_scan_manifest.tsv',sep='\t',index=False)
    return manifest
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--master-prefix',required=True);p.add_argument('--phenotype-dir',required=True);p.add_argument('--outdir',required=True)
    p.add_argument('--drug',action='append',choices=DRUGS);a=p.parse_args()
    prepare(a.master_prefix,a.phenotype_dir,a.outdir,a.drug or DRUGS)
