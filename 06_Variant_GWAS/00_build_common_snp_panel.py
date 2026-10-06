#!/usr/bin/env python3
"""Select one nuclear SNP panel for relatedness estimation and PC-based prediction."""
import argparse,json,subprocess,sys
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parent/'helpers'))
from master_inputs import load_master, counts, write_bed, NUCLEAR, sha256

def build(master, out_prefix, plink, chunk_size=4096):
    metadata,samples,x=load_master(master)
    indices=[]
    for lo in range(0,len(metadata),chunk_size):
        part=metadata.iloc[lo:lo+chunk_size]
        candidates=part.index[(part.type=='SNP')&part.chromosome.isin(NUCLEAR)].to_numpy()
        if not len(candidates): continue
        q=counts(x[candidates])
        indices.extend(candidates[(q['maf']>=.05)&(q['missing_fraction']<=.10)])
    indices=np.asarray(indices,dtype=np.int64)
    target=Path(out_prefix);target.parent.mkdir(parents=True,exist_ok=True)
    pre=target.parent/'common_nuclear_before_ld'
    selected=metadata.iloc[indices].copy();selected['cm']=0
    write_bed(pre,selected,samples,(x[indices[i:i+chunk_size]] for i in range(0,len(indices),chunk_size)))
    command=[str(plink),'--bfile',str(pre),'--chr-set','9','no-xy','--indep-pairwise','50','5','0.2','--out',str(target)]
    with open(str(target)+'.pruning_console.log','w') as log:
        subprocess.run(command,check=True,stdout=log,stderr=subprocess.STDOUT)
    retained=set(Path(str(target)+'.prune.in').read_text().splitlines())
    final=selected.loc[selected.id.isin(retained)].copy()
    if len(final)!=len(retained) or not len(final):raise ValueError('Invalid retained SNP identifiers')
    ix=final['index'].to_numpy(int)
    write_bed(target,final,samples,(x[ix[i:i+chunk_size]] for i in range(0,len(ix),chunk_size)))
    Path(str(target)+'.ids').write_text('\n'.join(final.id)+'\n')
    np.save(str(target)+'.variant_indices.npy',ix)
    final.to_csv(str(target)+'.metadata.tsv.gz',sep='\t',index=False)
    audit={'samples':len(samples),'before_ld':len(selected),'retained_snps':len(final),
           'maf_min':.05,'missing_max':.10,'ld_window_variants':50,'ld_step_variants':5,'ld_r2':.2,
           'frequency_cohort':'all master samples','nuclear_chromosomes':NUCLEAR,
           'master_metadata_sha256':sha256(str(master)+'.variants.tsv.gz'),
           'retained_ids_sha256':sha256(str(target)+'.ids'),'command':command,
           'PLINK_bed_effect_allele':'A1 is target ALT; A2 labels reference but codes all called non-target alleles',
           'pruning_order':'master normalized genomic order'}
    Path(str(target)+'.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(json.dumps(audit,indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--master-prefix',required=True);p.add_argument('--out-prefix',required=True)
    p.add_argument('--plink',default='plink');a=p.parse_args();build(a.master_prefix,a.out_prefix,a.plink)
