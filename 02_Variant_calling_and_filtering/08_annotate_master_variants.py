#!/usr/bin/env python3
"""Annotate normalized alleles once and export shared allele-to-gene assignments."""
import argparse,csv,gzip,hashlib,json,re,shutil,subprocess
from collections import Counter
from pathlib import Path
from urllib.parse import unquote
from cyvcf2 import VCF,Writer

def run(cmd,**kwargs):subprocess.run(list(map(str,cmd)),check=True,**kwargs)
def chrom_number(chrom):
    match=re.search(r'chr([1-7MR])',chrom)
    if not match:raise ValueError('Unrecognized reference contig: '+chrom)
    return {'M':8,'R':9}.get(match[1],int(match[1]) if match[1].isdigit() else None)
def attributes(text):return {k:unquote(v) for k,v in (x.split('=',1) for x in text.split(';') if '=' in x)}
def publication_symbol(gene_id,symbol):
    return 'FKS1' if gene_id=='C1_02420C_A' else symbol
def main():
    p=argparse.ArgumentParser(description=__doc__)
    for x in ['out-prefix','reference','gff','snpeff-jar','base-config']:p.add_argument('--'+x,required=True)
    p.add_argument('--java',default='java');p.add_argument('--bcftools',default='bcftools');p.add_argument('--genome',default='SC5314_A22_s07_m01_r183')
    p.add_argument('--reuse-raw-annotation',help='Use an existing raw SnpEff sites VCF with identical canonical row order');p.add_argument('--reference-cds');p.add_argument('--reference-proteins');p.add_argument('--threads',type=int,default=2)
    a=p.parse_args();stem=Path(a.out_prefix);ann=stem.parent/'annotation';db=ann/'data'/a.genome;db.mkdir(parents=True,exist_ok=True)
    contigs={l[1:].split()[0] for l in open(a.reference) if l.startswith('>')}
    gene_map={};gff_lines=[];records=[];transcript_types={};pseudo=set()
    for line in open(a.gff):
        if line.startswith('#'):gff_lines.append(line);continue
        f=line.rstrip('\n').split('\t')
        if len(f)!=9 or f[0] not in contigs:continue
        records.append(f);at=attributes(f[8])
        if f[2]=='pseudogene':pseudo.add(at['ID'])
        if f[2] in ['mRNA','tRNA','rRNA','snoRNA','snRNA','ncRNA']:
            for gid in at.get('Parent','').split(','):transcript_types[gid]=f[2]
    def biotype(gid,chrom):
        if gid in pseudo:return 'pseudogene'
        typ=transcript_types.get(gid,'ncRNA')
        if typ=='mRNA':return 'protein_coding'
        if typ=='tRNA':return 'Mt_tRNA' if 'chrM_' in chrom else 'non_coding'
        if typ=='rRNA':return 'Mt_rRNA' if 'chrM_' in chrom else 'rRNA'
        return typ if typ in ['snoRNA','snRNA'] else 'non_coding'
    for f in records:
        at=attributes(f[8]);original_type=f[2]
        if original_type in ['gene','pseudogene']:
            gid=at['ID'];bt=biotype(gid,f[0]);symbol=publication_symbol(gid,at.get('Gene',at.get('Name',gid)))
            description=at.get('Note','')
            if gid=='C1_02420C_A':description=re.sub(r'\bgsc1\b','FKS1',description,flags=re.IGNORECASE)
            gene_map[gid]={'gene_id':gid,'gene_symbol':symbol,'chromosome':chrom_number(f[0]),'chrom':f[0],'start':int(f[3]),'end':int(f[4]),'strand':f[6],'description':description,'source_biotype':'pseudogene' if gid in pseudo else transcript_types.get(gid,'ncRNA'),'snpeff_biotype':bt}
            # Explicit gene and transcript biotypes preserve the reference feature classes.
            f[2]='gene';f[8]+=f';gene_biotype={bt};biotype={bt}'
        elif original_type in ['mRNA','tRNA','rRNA','snoRNA','snRNA','ncRNA']:
            gid=at.get('Parent','').split(',')[0];bt=biotype(gid,f[0]);f[8]+=f';transcript_biotype={bt};gene_biotype={bt};biotype={bt}'
        gff_lines.append('\t'.join(f)+'\n')
    with open(stem.parent/'gene_map.tsv','w') as out:
        w=csv.DictWriter(out,delimiter='\t',fieldnames=list(next(iter(gene_map.values()))));w.writeheader();w.writerows(gene_map.values())
    (db/'genes.gff').write_text(''.join(gff_lines));shutil.copyfile(a.reference,db/'sequences.fa')
    config=Path(a.base_config).read_text()+f'\n{a.genome}.genome : Candida albicans SC5314 Assembly22 haplotype A\n'
    for chrom in sorted(contigs):config+=f'{a.genome}.{chrom}.codonTable : '+('Mold_Mitochondrial' if 'chrM_' in chrom else 'Alternative_Yeast_Nuclear')+'\n'
    config+='\ndata.dir = '+str((ann/'data').resolve())+'\n';cfg=ann/'snpEff.config';cfg.write_text(config)
    cmd=[a.java,'-Xmx4g','-jar',a.snpeff_jar,'build','-gff3','-treatAllAsProteinCoding','false','-noLog','-v','-c',cfg]
    if a.reference_cds:shutil.copyfile(a.reference_cds,db/'cds.fa')
    else:cmd+=['-noCheckCds']
    if a.reference_proteins:shutil.copyfile(a.reference_proteins,db/'protein.fa')
    else:cmd+=['-noCheckProtein']
    with open(ann/'build.log','w') as log:run(cmd+[a.genome],stdout=log,stderr=log)
    annotated=Path(str(stem)+'.annotated_sites.vcf.gz')
    raw=Path(a.reuse_raw_annotation) if a.reuse_raw_annotation else ann/'raw_snpeff_sites.vcf.gz'
    if not a.reuse_raw_annotation:
        with open(ann/'annotation.log','w') as log:
            cmd=[a.java,'-Xmx4g','-jar',a.snpeff_jar,'ann','-noLog','-noStats','-nodownload','-c',str(cfg),a.genome,str(stem)+'.sites.vcf.gz']
            process=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=log)
            compressed=subprocess.Popen([a.bcftools,'view','-Oz','--threads',str(a.threads),'-o',str(raw)],stdin=process.stdout,stderr=log)
            process.stdout.close()
            if compressed.wait() or process.wait():raise RuntimeError('SnpEff annotation failed')
    if raw.resolve()==annotated.resolve():raise ValueError('Raw and standardized annotation paths must differ')
    run([a.bcftools,'index','-f','-t',raw])
    raw_reader=VCF(str(raw));standardized=Writer(str(annotated),raw_reader,mode='wz')
    coding_effects={'missense_variant','synonymous_variant','frameshift_variant','stop_gained','stop_lost','start_lost','inframe_insertion','inframe_deletion','disruptive_inframe_insertion','disruptive_inframe_deletion','conservative_inframe_insertion','conservative_inframe_deletion','protein_altering_variant','stop_retained_variant','start_retained_variant','initiator_codon_variant','incomplete_terminal_codon_variant','coding_sequence_variant'}
    adjustments=Counter()
    with gzip.open(ann/'noncoding_consequence_adjustments.tsv.gz','wt',compresslevel=2) as audit_file:
        audit_writer=csv.writer(audit_file,delimiter='\t');audit_writer.writerow(['variant_id','annotation_index','gene_id','reference_biotype','original_ann','standardized_ann'])
        for v in raw_reader:
            entries=(v.INFO.get('ANN') or '').split(',');changed=False;noncoding_changed=False
            for j,entry in enumerate(entries):
                f=entry.split('|');f+=['']*max(0,16-len(f));gene=gene_map.get(f[4])
                symbol=publication_symbol(f[4],f[3])
                if f[3]!=symbol:
                    f[3]=symbol;entries[j]='|'.join(f);changed=True
                if not gene or gene['source_biotype']=='mRNA':continue
                original=list(f);effects=f[1].split('&');protein_effects=set(effects)&coding_effects
                if protein_effects:
                    valid=[effect for effect in effects if effect not in coding_effects]
                    if valid:
                        # Mixed functional consequences require an explicit curated impact mapping.
                        raise ValueError('Mixed coding/noncoding consequence requires review: '+v.ID+' '+f[1])
                    f[1]='non_coding_transcript_exon_variant';f[2]='MODIFIER';adjustments['protein_effect_entries']+=1
                if f[10]:f[10]='';adjustments['protein_hgvs_entries']+=1
                f[12]='';f[13]=''
                if f!=original:
                    entries[j]='|'.join(f);changed=True;noncoding_changed=True;adjustments['annotation_entries']+=1
                    audit_writer.writerow([v.ID,j,f[4],gene['source_biotype'],entry,entries[j]])
            if changed:v.INFO['ANN']=','.join(entries)
            if noncoding_changed:adjustments['variant_records']+=1
            standardized.write_record(v)
    standardized.close();raw_reader.close();run([a.bcftools,'index','-f','-t',annotated])
    full_fields=['index','id','gene_id','gene_symbol','impact','effect','feature_id','transcript_biotype','hgvs_c','hgvs_p','annotation_allele','warnings']
    selected_fields=['index','id','variant_idx','variant_id','gene_id','gene_symbol','impact','effect','hgvs_c','hgvs_p','protein_change','chromosome','bp']
    counts=Counter();selected_path=stem.parent/'selected_gene_variants.tsv.gz'
    metadata_handle=gzip.open(str(stem)+'.variants.tsv.gz','rt');metadata=csv.DictReader(metadata_handle,delimiter='\t')
    with gzip.open(str(stem)+'.gene_annotations.tsv.gz','wt',compresslevel=2) as full,gzip.open(selected_path,'wt',compresslevel=2) as selected:
        w=csv.DictWriter(full,delimiter='\t',fieldnames=full_fields);w.writeheader();s=csv.DictWriter(selected,delimiter='\t',fieldnames=selected_fields);s.writeheader()
        for i,v in enumerate(VCF(str(annotated))):
            expected=next(metadata)
            if expected['id']!=v.ID or int(expected['index'])!=i:raise ValueError('Annotation and dosage row order differ')
            counts['annotated_alleles']+=1; choices=[]
            if len(v.ALT)!=1:raise ValueError('Annotation input must be biallelic')
            annotations=v.INFO.get('ANN') or ''
            for entry in annotations.split(','):
                if not entry:continue
                f=entry.split('|');f+=['']*max(0,16-len(f));gene=gene_map.get(f[4]);symbol=gene['gene_symbol'] if gene else f[3]
                row=dict(zip(full_fields,[i,v.ID,f[4],symbol,f[2],f[1],f[6],f[7],f[9],f[10],f[0],f[15]]));w.writerow(row);counts['annotation_records']+=1
                if f[2] in ['HIGH','MODERATE'] and gene:
                    # Each standardized record has one ALT; all ANN entries are specific to that allele.
                    choices.append(row)
            if choices:
                r=min(choices,key=lambda r:(0 if r['impact']=='HIGH' else 1,0 if r['hgvs_p'] else 1,r['gene_id'],r['effect'],r['hgvs_p'],r['feature_id']))
                g=gene_map[r['gene_id']]
                s.writerow({'index':i,'id':v.ID,'variant_idx':i,'variant_id':v.ID,'gene_id':r['gene_id'],'gene_symbol':r['gene_symbol'],'impact':r['impact'],'effect':r['effect'],'hgvs_c':r['hgvs_c'],'hgvs_p':r['hgvs_p'],'protein_change':r['hgvs_p'],'chromosome':g['chromosome'],'bp':v.POS})
                counts['selected_HIGH_MODERATE_alleles']+=1;counts['selected_mitochondrial_alleles' if g['chromosome']==8 else 'selected_nuclear_alleles']+=1
    assert next(metadata,None) is None;metadata_handle.close()
    final=Path(str(stem)+'.vcf.gz')
    annotate=subprocess.Popen([a.bcftools,'annotate','-a',str(annotated),'-c','INFO/ANN','--pair-logic','exact','-Ou',str(stem)+'.unannotated.bcf'],stdout=subprocess.PIPE)
    fill=subprocess.Popen([a.bcftools,'+fill-tags','-Oz','--threads',str(a.threads),'-o',str(final),'--','-t','AC,AN,AF,NS,F_MISSING,MAF'],stdin=annotate.stdout)
    annotate.stdout.close()
    if fill.wait() or annotate.wait():raise RuntimeError('Master annotation/tag export failed')
    run([a.bcftools,'index','-f','-t',final])
    audit_path=Path(str(stem)+'.audit.json');audit=json.loads(audit_path.read_text());audit['annotation_status']='complete';audit['annotation']={'counts':dict(counts),'snpeff_version':subprocess.check_output([a.java,'-jar',a.snpeff_jar,'-version']).decode().strip(),'genome':a.genome,'gene_map_count':len(gene_map),'nuclear_codon_table':'Alternative_Yeast_Nuclear (NCBI 12)','mitochondrial_codon_table':'Mold_Mitochondrial (NCBI 4)','external_CDS_check':bool(a.reference_cds),'external_protein_check':bool(a.reference_proteins),'noncoding_consequence_adjustments':dict(adjustments),'gff_sha256':hashlib.sha256(Path(a.gff).read_bytes()).hexdigest(),'gene_name_convention':{'C1_02420C_A':'FKS1'},'reference_biotypes_explicit':True,'nuclear_tRNA_snpeff_biotype':'non_coding','selected_annotation_priority':'HIGH before MODERATE; mapped reference gene; nonempty protein HGVS; gene ID/effect/protein HGVS/transcript ID lexicographic tie-break'};audit_path.write_text(json.dumps(audit,indent=2)+'\n');print(json.dumps(audit['annotation'],indent=2))
if __name__=='__main__':main()
