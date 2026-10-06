"""Read the shared sequence-variant catalogue and preserve target-ALT coding."""
from pathlib import Path
import hashlib
import re
import numpy as np
import pandas as pd

NUCLEAR = [1, 2, 3, 4, 5, 6, 7, 9]
DRUGS = ['fluconazole','voriconazole','posaconazole','caspofungin','micafungin','anidulafungin']

def chromosome_number(value):
    s = str(value)
    if s.isdigit() and 1 <= int(s) <= 9:
        return int(s)
    m = re.search(r'chr([1-7RM])', s, re.I)
    if not m:
        raise ValueError('Unrecognized reference chromosome: ' + s)
    return {'R':9, 'M':8}.get(m[1].upper(), int(m[1]) if m[1].isdigit() else None)

def sha256(path):
    h = hashlib.sha256()
    with open(path,'rb') as f:
        for chunk in iter(lambda:f.read(8*1024**2), b''):
            h.update(chunk)
    return h.hexdigest()

def load_master(prefix):
    prefix = str(prefix)
    # Attribute access is deliberately avoided for pandas method names.
    samples = pd.read_csv(prefix+'.samples.tsv',sep='\t')['sample'].astype(str).tolist()
    metadata = pd.read_csv(prefix+'.variants.tsv.gz',sep='\t',usecols=['index','id','chrom','pos','ref','alt','type'])
    metadata['chromosome'] = metadata.chrom.map(chromosome_number)
    x = np.load(prefix+'.dosage.npy', mmap_mode='r')
    if x.shape != (len(metadata),len(samples)) or x.dtype != np.int8:
        raise ValueError('Master dosage dimensions or type disagree with its metadata')
    if len(samples) != len(set(samples)) or metadata.id.duplicated().any():
        raise ValueError('Duplicate sample or canonical allele identifiers')
    if not np.array_equal(metadata['index'],np.arange(len(metadata))):
        raise ValueError('Master allele indices must be consecutive and zero-based')
    return metadata, samples, x

def counts(x):
    n_called = (x >= 0).sum(axis=1)
    ac = np.maximum(x,0).sum(axis=1)
    af = np.divide(ac, 2*n_called, out=np.full(len(x),np.nan), where=n_called>0)
    return {'n_called':n_called, 'n_missing':x.shape[1]-n_called,
            'alt_frequency':af, 'maf':np.minimum(af,1-af),
            'missing_fraction':1-n_called/x.shape[1],
            'n_alt_carriers':(x>0).sum(axis=1),
            'minor_allele_carriers':np.where(af<=.5,(x>0).sum(axis=1),((x>=0)&(x<2)).sum(axis=1)),
            'n_distinct_called_dosages':sum((x==v).any(axis=1).astype(int) for v in range(3))}

def pack_dosages(x):
    if not np.isin(x,[-1,0,1,2]).all():
        raise ValueError('Dosages must be -1, 0, 1 or 2')
    lookup = np.array([1,3,2,0],dtype=np.uint8)
    codes = lookup[x+1]
    padded = np.zeros((len(x),4*((x.shape[1]+3)//4)),dtype=np.uint8)
    padded[:,:x.shape[1]]=codes
    return ((padded.reshape(len(x),-1,4).astype(np.uint16) << np.array([0,2,4,6]))
            .sum(axis=2).astype(np.uint8)).tobytes()

def write_bed(prefix, metadata, samples, blocks):
    prefix=Path(prefix); prefix.parent.mkdir(parents=True,exist_ok=True)
    with open(str(prefix)+'.bed','wb') as f:
        f.write(bytes([108,27,1]))
        for x in blocks:
            f.write(pack_dosages(x))
    metadata[['chromosome','id','cm','pos','alt','ref']].to_csv(str(prefix)+'.bim',sep='\t',header=False,index=False)
    with open(str(prefix)+'.fam','w') as f:
        for sample in samples:
            f.write(f'{sample}\t{sample}\t0\t0\t0\t-9\n')
    expected=3+len(metadata)*((len(samples)+3)//4)
    if Path(str(prefix)+'.bed').stat().st_size!=expected:
        raise ValueError('Written BED dimensions disagree with metadata')

def phenotype_cohort(directory, drug, samples):
    file=Path(directory)/f'pheno.mic.{drug}.txt'
    p=pd.read_csv(file,sep=r'\s+',header=None,names=['FID','sample','log2_mic'],dtype={'FID':str,'sample':str})
    if p['sample'].duplicated().any() or not p.FID.equals(p['sample']):
        raise ValueError('Phenotype identifiers are not unique matching FID/IID pairs')
    p=p.loc[np.isfinite(p.log2_mic)&p.log2_mic.ne(-9)].set_index('sample')
    unknown=set(p.index)-set(samples)
    if unknown:
        raise ValueError('Phenotype samples absent from master: '+str(sorted(unknown)))
    idx=np.array([i for i,s in enumerate(samples) if s in p.index])
    ordered=[samples[i] for i in idx]
    p=p.loc[ordered].reset_index()
    return idx,ordered,p
