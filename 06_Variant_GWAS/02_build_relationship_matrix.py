#!/usr/bin/env python3
"""Build full-nuclear and chromosome-excluded GRMs from the identical shared SNP set."""
import argparse,json,sys
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parent/'helpers'))
from master_inputs import load_master,NUCLEAR,sha256

def write_grm(prefix,kernel,nmarkers,samples):
    lower=np.tril_indices(len(samples))
    kernel[lower].astype('<f4').tofile(str(prefix)+'.grm.bin')
    np.full(len(lower[0]),nmarkers,dtype='<f4').tofile(str(prefix)+'.grm.N.bin')
    Path(str(prefix)+'.grm.id').write_text(''.join(f'{s}\t{s}\n' for s in samples))

def build(master,panel,outdir,chunk_size=2048):
    metadata,samples,x=load_master(master);panel=Path(panel);outdir=Path(outdir);outdir.mkdir(parents=True,exist_ok=True)
    indices=np.load(str(panel)+'.variant_indices.npy')
    chosen=metadata.iloc[indices]
    if chosen.id.tolist()!=Path(str(panel)+'.ids').read_text().splitlines():
        raise ValueError('Shared panel indices and identifiers differ')
    if not chosen.type.eq('SNP').all() or not chosen.chromosome.isin(NUCLEAR).all():
        raise ValueError('The common panel must contain only nuclear SNPs')
    sums={c:np.zeros((len(samples),len(samples)),np.float64) for c in NUCLEAR};num={c:0 for c in NUCLEAR}
    for chromosome,part in chosen.groupby('chromosome'):
        ix=part['index'].to_numpy(int)
        for lo in range(0,len(ix),chunk_size):
            d=x[ix[lo:lo+chunk_size]].astype(float);missing=d<0;d[missing]=np.nan
            p=np.nanmean(d,axis=1)/2
            if not (np.minimum(p,1-p)>=.05-1e-12).all() or not (missing.mean(axis=1)<=.10+1e-12).all():
                raise ValueError('Shared SNP set violates baseline criteria; no additional SNP removal is performed')
            z=(d-2*p[:,None])/np.sqrt(2*p*(1-p))[:,None];z[missing]=0
            sums[chromosome]+=z.T@z;num[chromosome]+=len(d)
    total=sum(sums.values());m=sum(num.values())
    full=outdir/'common_nuclear';write_grm(full,total/m,m,samples)
    manifest=[]
    for chromosome in NUCLEAR:
        prefix=outdir/f'loco_chr{chromosome}';n=m-num[chromosome]
        write_grm(prefix,(total-sums[chromosome])/n,n,samples)
        manifest.append({'excluded_chr':chromosome,'grm_prefix':str(prefix.resolve()),'retained_snps':n,'adjustment':'LOCO'})
    # Mitochondrial markers use all nuclear SNPs, since none are on their chromosome.
    manifest.append({'excluded_chr':8,'grm_prefix':str(full.resolve()),'retained_snps':m,'adjustment':'full nuclear'})
    pd.DataFrame(manifest).to_csv(outdir/'loco_manifest.tsv',sep='\t',index=False)
    summary={'samples':len(samples),'retained_snps':m,'snps_by_chromosome':num,'shared_panel':str(panel.resolve()),
             'shared_ids_sha256':sha256(str(panel)+'.ids'),'formula':'K=Z.T@Z/M; Z=(ALT dosage-2p)/sqrt(2p(1-p)); missing Z=0',
             'frequency_cohort':'all master samples','additional_filtering':False,'minimum_eigenvalue':float(np.linalg.eigvalsh(total/m)[0]),
             'minimum_diagonal':float(np.diag(total/m).min()),'maximum_diagonal':float(np.diag(total/m).max())}
    (outdir/'grm_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--master-prefix',required=True)
    p.add_argument('--shared-panel-prefix',required=True);p.add_argument('--outdir',required=True);a=p.parse_args()
    build(a.master_prefix,a.shared_panel_prefix,a.outdir)
