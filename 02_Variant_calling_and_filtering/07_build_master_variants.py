#!/usr/bin/env python3
"""Create a normalized sequence-SNP/indel catalogue and target-ALT dosage matrix."""
import argparse,csv,fcntl,gzip,hashlib,json,re,subprocess,sys
from collections import Counter
from pathlib import Path
import numpy as np
from cyvcf2 import VCF,Writer

def run(cmd, **kw):
    print(' '.join(map(str,cmd)),file=sys.stderr,flush=True)
    subprocess.run(list(map(str,cmd)),check=True,**kw)

def dosage(gt, alt_index=1):
    """Count one target allele only for fully called diploid genotypes."""
    a=np.asarray(gt)
    out=np.full(a.shape[0],-1,dtype=np.int8)
    alleles=a[:,:-1]
    # cyvcf2 pads shorter genotypes with -2 when records contain mixed ploidies.
    valid=((alleles!=-2).sum(axis=1)==2)&np.all((alleles>=0)|(alleles==-2),axis=1)
    out[valid]=(alleles[valid]==alt_index).sum(axis=1)
    return out

def eligible(ref,alt):
    if not re.fullmatch('[ACGTN]+',ref.upper()) or not re.fullmatch('[ACGTN]+',alt.upper()):
        return False
    return (len(ref)==len(alt)==1 and ref.upper()!=alt.upper()) or len(ref)!=len(alt)

def hash_file(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(8*1024**2),b''): h.update(b)
    return h.hexdigest()

def read_samples(path):
    rows=[line.split() for line in Path(path).read_text().splitlines() if line.strip()]
    if not rows or any(len(row)!=1 for row in rows):
        raise ValueError('Sample list must contain one identifier per nonempty line')
    samples=[row[0] for row in rows]
    if len(samples)!=len(set(samples)):
        raise ValueError('Sample identifiers must be unique')
    return samples

def check_sample_order(path, expected):
    reader=VCF(str(path))
    try:
        if reader.samples!=expected:
            raise ValueError('Cached genotype sample order differs from the requested sample list: '+str(path))
    finally:
        reader.close()

def artifact_hashes(paths, parent):
    return {str(path.relative_to(parent)):hash_file(path) for path in paths}

def validate_artifacts(artifacts, parent):
    if not artifacts:
        raise ValueError('Cache manifest does not contain artifact checksums')
    for name,expected in artifacts.items():
        path=Path(name)
        if path.is_absolute() or '..' in path.parts:
            raise ValueError('Cache artifact names must be relative to the output directory')
        path=parent/path
        if not path.is_file() or hash_file(path)!=expected:
            raise ValueError('Cache artifact is missing or has a different checksum: '+str(path))

def write_manifest(path, manifest):
    temporary=path.with_suffix('.tmp')
    temporary.write_text(json.dumps(manifest,indent=2)+'\n')
    temporary.replace(path)

def prepare_cache(stem, resume, binding):
    work=Path(str(stem)+'.master_work')
    manifest_path=work/'manifest.json'
    if resume:
        if not manifest_path.is_file():
            raise ValueError('Resume requires a completed preparation manifest for this output prefix')
        manifest=json.loads(manifest_path.read_text())
        if manifest.get('format_version')!=1 or manifest.get('output_prefix')!=stem.name or manifest.get('inputs')!=binding:
            raise ValueError('Resume inputs, sample order or tool version differ from the cache manifest')
        if manifest.get('stage') not in {'prepared','complete'}:
            raise ValueError('Cache manifest does not describe a completed preparation stage')
        validate_artifacts(manifest.get('prepared_artifacts'),stem.parent)
        for path in [Path(str(stem)+'.original_genotypes.bcf'),work/'eligible.normalized.bcf',work/'eligible.sorted.bcf']:
            check_sample_order(path,binding['sample_order'])
        stored=(Path(str(stem)+'.samples.tsv').read_text().splitlines())
        if stored!=['sample']+binding['sample_order']:
            raise ValueError('Cached dosage sample table differs from the requested sample order')
        if manifest['stage']=='complete':
            validate_artifacts(manifest.get('output_artifacts'),stem.parent)
            audit=json.loads(Path(str(stem)+'.audit.json').read_text())
            if any(audit.get(key)!=binding[key] for key in ['source_sha256','reference_sha256','samples_sha256']):
                raise ValueError('Catalogue audit input hashes differ from the cache manifest')
            if audit.get('counts')!=manifest.get('counts'):
                raise ValueError('Catalogue audit counts differ from the cache manifest')
        return work,manifest
    existing=[path for path in stem.parent.iterdir() if path.name.startswith(stem.name+'.')]
    if existing:
        raise ValueError('Output prefix already has files; use --resume or choose a new output prefix')
    work.mkdir(parents=True,exist_ok=False)
    return work,{'format_version':1,'output_prefix':stem.name,'inputs':binding,'stage':'preparing'}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--vcf',required=True);p.add_argument('--samples',required=True)
    p.add_argument('--reference',required=True);p.add_argument('--out-prefix',required=True)
    p.add_argument('--resume',action='store_true',help='Validate and reuse completed prefix-specific preparation or catalogue outputs');p.add_argument('--bcftools',default='bcftools');p.add_argument('--threads',type=int,default=2)
    a=p.parse_args()
    if a.threads<1:p.error('--threads must be positive')
    stem=Path(a.out_prefix).expanduser().resolve();stem.parent.mkdir(parents=True,exist_ok=True)
    orig=Path(str(stem)+'.original_genotypes.bcf');b=a.bcftools
    wanted=read_samples(a.samples)
    binding={label+'_sha256':hash_file(path) for label,path in [('source',a.vcf),('reference',a.reference),('samples',a.samples)]}
    binding['sample_order']=wanted
    binding['bcftools_version']=subprocess.check_output([b,'--version']).decode().splitlines()[:2]
    if a.resume:
        work=Path(str(stem)+'.master_work')
        if not (work/'manifest.json').is_file():
            raise ValueError('Resume requires a completed preparation manifest for this output prefix')
    else:
        work,manifest=prepare_cache(stem,False,binding)
    lock=(work/'run.lock').open('a')
    try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError:raise ValueError('Another process is using this output prefix')
    if a.resume:
        # Validate the cached artifacts while holding the prefix lock.
        work,manifest=prepare_cache(stem,True,binding)
        if manifest['stage']=='complete':
            print(json.dumps(manifest['counts'],indent=2))
            return
    else:
        with open(work/'subset.log','w') as log:
            p1=subprocess.Popen([b,'view','-S',a.samples,'-Ou',a.vcf],stdout=subprocess.PIPE,stderr=log)
            p2=subprocess.Popen([b,'annotate','-x','^INFO/MQ,^FORMAT/GT','-Ob','--threads',str(a.threads),'-o',str(orig)],stdin=p1.stdout,stderr=log)
            p1.stdout.close();second_status=p2.wait();first_status=p1.wait()
            if first_status or second_status:raise RuntimeError('Sample selection failed; inspect subset.log')
        run([b,'index','-f',orig])
    source=VCF(str(orig));samples=source.samples
    if samples!=wanted:raise ValueError('Selected VCF sample order differs from the requested sample list')
    if not a.resume:
        with open(str(stem)+'.samples.tsv','w') as f:f.write('sample\n'+'\n'.join(samples)+'\n')
    eligible_bcf=work/'eligible.normalized.bcf';sorted_bcf=work/'eligible.sorted.bcf'
    counts=Counter(); origin_path=Path(str(stem)+'.original_records.tsv.gz')
    if not a.resume:
        hdr=['##fileformat=VCFv4.2']
        hdr.extend(x for x in source.raw_header.splitlines() if x.startswith('##contig=') or x.startswith('##FILTER='))
        hdr.extend(['##FORMAT=<ID=GT,Number=1,Type=String,Description="Target ALT copies: zero represents any called non-target allele">',
          '##INFO=<ID=SOURCE_REC,Number=1,Type=String,Description="One-based original cohort record ordinal">',
          '##INFO=<ID=SOURCE_ALT_IDX,Number=1,Type=Integer,Description="Original one-based ALT index">',
          '##INFO=<ID=MQ,Number=1,Type=Float,Description="Mapping quality retained from original joint call">',
          '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+'\t'.join(samples)])
        with open(work/'normalization.log','w') as log,gzip.open(origin_path,'wt',compresslevel=2) as provenance:
            norm=subprocess.Popen([b,'norm','-f',a.reference,'-c','e','--old-rec-tag','NORM_ORIG','-Ob','--threads',str(a.threads),'-o',str(eligible_bcf)],stdin=subprocess.PIPE,text=True,stderr=log)
            norm.stdin.write('\n'.join(hdr)+'\n');wr=csv.writer(provenance,delimiter='\t');wr.writerow(['source_record','chrom','pos','original_id','ref','alts','alt_index','target_alt'])
            lookup=np.array(['./.','0/0','0/1','1/1'])
            for n,v in enumerate(source,1):
                counts['original_records']+=1;gt=v.genotype.array(); src=f'src{n:07d}'
                for j,alt in enumerate(v.ALT,1):
                    if not eligible(v.REF,alt):counts['ineligible_alt_alleles']+=1;continue
                    counts['eligible_source_alleles']+=1
                    d=dosage(gt,j); encoded='\t'.join(lookup[d+1])
                    wr.writerow([src,v.CHROM,v.POS,v.ID or '.',v.REF,','.join(v.ALT),j,alt])
                    info=f'SOURCE_REC={src};SOURCE_ALT_IDX={j}';mq=v.INFO.get('MQ')
                    if mq is not None:info+=f';MQ={mq}'
                    norm.stdin.write(f'{v.CHROM}\t{v.POS}\t{src}:{j}\t{v.REF.upper()}\t{alt.upper()}\t.\t.\t{info}\tGT\t{encoded}\n')
                if n%100000==0:print(f'Extracted {n:,} source records',file=sys.stderr,flush=True)
            norm.stdin.close()
            if norm.wait():raise RuntimeError('Reference normalization failed; inspect normalization.log')
        (work/'extraction_counts.json').write_text(json.dumps(counts,indent=2))
    else:counts.update(json.loads((work/'extraction_counts.json').read_text()))
    source.close()
    if not a.resume:
        run([b,'sort','-m','1G','-Ob','-o',sorted_bcf,eligible_bcf])
        run([b,'index','-f',sorted_bcf])
        for path in [eligible_bcf,sorted_bcf]:check_sample_order(path,wanted)
        prepared=[orig,Path(str(orig)+'.csi'),Path(str(stem)+'.samples.tsv'),origin_path,
                  eligible_bcf,sorted_bcf,Path(str(sorted_bcf)+'.csi'),work/'extraction_counts.json',work/'normalization.log']
        manifest['prepared_artifacts']=artifact_hashes(prepared,stem.parent)
        manifest['stage']='prepared';write_manifest(work/'manifest.json',manifest)
    norm=VCF(str(sorted_bcf))
    nmax=int(subprocess.check_output([b,'index','-n',str(sorted_bcf)]).strip());ns=len(samples)
    raw=work/'dosage.raw';dm=np.memmap(raw,dtype=np.int8,mode='w+',shape=(nmax,ns))
    canonical=Path(str(stem)+'.unannotated.bcf');writer=Writer(str(canonical),norm,mode='wb')
    metadata=Path(str(stem)+'.variants.tsv.gz');mapping=Path(str(stem)+'.normalization_map.tsv.gz')
    with gzip.open(metadata,'wt',compresslevel=2) as mt,gzip.open(mapping,'wt',compresslevel=2) as mp:
        mw=csv.writer(mt,delimiter='\t');mw.writerow(['index','id','chrom','pos','ref','alt','type','is_mito','ac','an','af','missing_fraction','carrier_count'])
        pw=csv.writer(mp,delimiter='\t');pw.writerow(['index','id','source_record','source_alt_index','normalization_original'])
        position=None; buffer={};n=0
        def flush():
            nonlocal n
            for key,(v,d,origins) in sorted(buffer.items()):
                for other in origins:
                    pw.writerow([n,key,*other])
                v.ID=key;writer.write_record(v);dm[n]=d
                ac=int(d[d>=0].sum());an=int((d>=0).sum()*2);kind='SNP' if len(v.REF)==len(v.ALT[0])==1 else 'INDEL'
                mito=int('chrM_' in v.CHROM);counts[kind]+=1;counts['mitochondrial_'+kind if mito else 'nuclear_'+kind]+=1
                mw.writerow([n,key,v.CHROM,v.POS,v.REF,v.ALT[0],kind,mito,ac,an,ac/an if an else '.',float((d<0).mean()),int((d>0).sum())]);n+=1
        for v in norm:
            pos=(v.CHROM,v.POS)
            if position is not None and pos!=position:flush();buffer={}
            position=pos;key=f'{v.CHROM}:{v.POS}:{v.REF}:{v.ALT[0]}';d=dosage(v.genotype.array())
            origin=(v.INFO.get('SOURCE_REC'),v.INFO.get('SOURCE_ALT_IDX'),v.INFO.get('NORM_ORIG') or '.')
            if key in buffer:
                if not np.array_equal(buffer[key][1],d):raise RuntimeError(f'Conflicting normalized duplicate genotypes: {key}')
                buffer[key][2].append(origin);counts['collapsed_identical_duplicate_alleles']+=1
            else:buffer[key]=(v,d,[origin])
        flush();writer.close();norm.close();dm.flush()
    out=np.lib.format.open_memmap(str(stem)+'.dosage.npy',mode='w+',dtype=np.int8,shape=(n,ns))
    for lo in range(0,n,10000):out[lo:lo+10000]=dm[lo:min(lo+10000,n)]
    out.flush();del out,dm
    counts.update({'samples':ns,'canonical_alleles':n,'dosage_missing_value':-1})
    run([b,'index','-f',canonical]);run([b,'view','-G','-Oz','--threads',a.threads,'-o',str(stem)+'.sites.vcf.gz',canonical]);run([b,'index','-f','-t',str(stem)+'.sites.vcf.gz'])
    audit={'counts':dict(counts),'source_vcf':str(Path(a.vcf).resolve()),'source_sha256':hash_file(a.vcf),'reference':str(Path(a.reference).resolve()),'reference_sha256':hash_file(a.reference),'samples_sha256':hash_file(a.samples),'bcftools':subprocess.check_output([b,'--version']).decode().splitlines()[:2],'genotype_rule':'Number of target ALT copies in complete diploid calls; -1 for partial/full missing or non-diploid; other called ALTs contribute zero','eligibility':'Original sequence SNP/unequal-length indel alleles; no global frequency, missingness, MQ or FILTER restriction; no MNP atomization','annotation_status':'pending'}
    numbers=(work/'normalization.log').read_text().splitlines()[-1].split('\t')[-1].split('/')
    audit['normalization']={'realigned_allele_records':int(numbers[3]),'reference_mismatches':int(numbers[4]),'duplicate_conflicts':0}
    audit['counts'].setdefault('collapsed_identical_duplicate_alleles',0)
    Path(str(stem)+'.audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    completed=[canonical,Path(str(canonical)+'.csi'),metadata,mapping,Path(str(stem)+'.dosage.npy'),
               Path(str(stem)+'.sites.vcf.gz'),Path(str(stem)+'.sites.vcf.gz.tbi')]
    manifest['output_artifacts']=artifact_hashes(completed,stem.parent)
    manifest['counts']=audit['counts'];manifest['stage']='complete';write_manifest(work/'manifest.json',manifest)
    print(json.dumps(audit['counts'],indent=2))
if __name__=='__main__':main()
