"""Exploratory gene-specific joint co-carriage adjustment, with PCA clusters only."""
from pathlib import Path
import collections,csv,gzip,json,math,subprocess,argparse,shutil

HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--gene-results', type=Path, required=True)
parser.add_argument('--gene-calls', type=Path, required=True)
parser.add_argument('--outdir', type=Path, required=True)
parser.add_argument('--rscript', default='Rscript')
parser.add_argument('--thresholds', type=float, nargs='+', default=[.20,.15,.25,.30])
args = parser.parse_args()
if not shutil.which(args.rscript): parser.error('Rscript is required; install R or pass --rscript')
if not args.thresholds or any(not math.isfinite(t) or t < 0 or t > 1 for t in args.thresholds): parser.error('Thresholds must be finite numbers between zero and one')
if len(set(args.thresholds)) != len(args.thresholds): parser.error('Thresholds must be distinct')
THRESHOLDS = tuple(args.thresholds)
if .20 not in THRESHOLDS: parser.error('Include the primary threshold 0.20')
OUT = args.outdir.resolve()
OUT.mkdir(parents=True, exist_ok=True)

def read(p):
    with (gzip.open(p,'rt',newline='') if p.suffix=='.gz' else p.open(newline='')) as f:
        return list(csv.DictReader(f))

def write(name,rows,fields=None):
    with (OUT/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields or list(rows[0]));w.writeheader();w.writerows(rows)

def threshold_label(value):
    """Preserve study IDs while retaining other thresholds without rounding."""
    short = f'{value:.2f}'
    return short if float(short) == value else repr(value)

original=read(args.gene_results)
assert len({(r['drug'],r['gene_id']) for r in original}) == len(original), 'Duplicate gene–drug hypotheses'
qualified=[r for r in original if r['pca_cluster_q_primary'] and float(r['pca_cluster_q_primary'])<.05 and float(r['pca_cluster_or'])>1]
if not qualified: raise ValueError('No qualifying gene associations')
DRUGS=tuple(dict.fromkeys(r['drug'] for r in original))
qual_by_drug={d:[r for r in qualified if r['drug']==d] for d in DRUGS}
keyset={(r['drug'],r['gene_id']) for r in qualified}
calls=[r for r in read(args.gene_calls) if (r['drug'],r['gene_id']) in keyset]
by_drug=collections.defaultdict(dict)
for r in calls:
    d=by_drug[r['drug']].setdefault(r['sample'],{'sample':r['sample'],'phenotype':r['phenotype'],'cluster':r['cluster_label'],'calls':{}})
    assert d['phenotype']==r['phenotype'] and d['cluster']==r['cluster_label']
    assert r['gene_id'] not in d['calls']
    d['calls'][r['gene_id']]=int(r['carrier'])

# Within-PCA-cluster residual Pearson correlation of binary carriage calls.
# Selection pools both phenotypes and does not consult target-outcome p values.
pairs=[];pairmap={}
for drug,genes in qual_by_drug.items():
    samples=list(by_drug[drug].values())
    for i,a in enumerate(genes):
        for b in genes[i+1:]:
            groups=collections.defaultdict(list)
            for s in samples:
                x,y=s['calls'][a['gene_id']],s['calls'][b['gene_id']]
                if x>=0 and y>=0:groups[s['cluster']].append((x,y,s['phenotype']))
            xx=yy=xy=0.
            for values in groups.values():
                mx=sum(x for x,y,p in values)/len(values);my=sum(y for x,y,p in values)/len(values)
                xx+=sum((x-mx)**2 for x,y,p in values);yy+=sum((y-my)**2 for x,y,p in values)
                xy+=sum((x-mx)*(y-my) for x,y,p in values)
            corr=xy/math.sqrt(xx*yy) if xx>0 and yy>0 else None
            values=[v for group in groups.values() for v in group]
            row=dict(drug=drug,gene_1=a['gene'],gene_id_1=a['gene_id'],gene_2=b['gene'],gene_id_2=b['gene_id'],
                     jointly_called=len(values),R_called=sum(p=='R' for x,y,p in values),S_called=sum(p=='S' for x,y,p in values),
                     both_carriers=sum(x==y==1 for x,y,p in values),pca_residual_correlation=corr,
                     selected_primary=corr is not None and abs(corr)>=0.20)
            pairs.append(row);pairmap[(drug,a['gene_id'],b['gene_id'])]=corr;pairmap[(drug,b['gene_id'],a['gene_id'])]=corr
write('pairwise_pca_carriage_correlations.csv',pairs,fields=['drug','gene_1','gene_id_1','gene_2','gene_id_2','jointly_called','R_called','S_called','both_carriers','pca_residual_correlation','selected_primary'])

specs=[];modelrows=[]
for threshold in THRESHOLDS:
    for drug,genes in qual_by_drug.items():
        for target in genes:
            gid=target['gene_id']
            neighbors=[g for g in genes if g['gene_id']!=gid and pairmap.get((drug,gid,g['gene_id'])) is not None and abs(pairmap[(drug,gid,g['gene_id'])])>=threshold]
            neighbors.sort(key=lambda g:(-abs(pairmap[(drug,gid,g['gene_id'])]),g['gene_id']))
            mid=f'{drug}|{gid}|{threshold_label(threshold)}'
            eligible=list(by_drug[drug].values())
            complete=[s for s in eligible if all(s['calls'][g]>=0 for g in [gid]+[g['gene_id'] for g in neighbors])]
            row=dict(model_id=mid,threshold=threshold,drug=drug,gene=target['gene'],gene_id=gid,
                     original_pca_q=target['pca_cluster_q_primary'],original_pca_or=target['pca_cluster_or'],
                     n_partners=len(neighbors),partners=';'.join(g['gene'] for g in neighbors),partner_ids=';'.join(g['gene_id'] for g in neighbors),
                     partner_correlations=';'.join(f"{pairmap[(drug,gid,g['gene_id'])]:.10g}" for g in neighbors),
                     R_eligible=sum(s['phenotype']=='R' for s in eligible),S_eligible=sum(s['phenotype']=='S' for s in eligible),
                     R_called=sum(s['phenotype']=='R' for s in complete),S_called=sum(s['phenotype']=='S' for s in complete))
            row['R_omitted']=row['R_eligible']-row['R_called'];row['S_omitted']=row['S_eligible']-row['S_called']
            specs.append(row)
            for s in complete:
                pattern=''.join(str(s['calls'][g['gene_id']]) for g in neighbors) or 'none'
                modelrows.append(dict(model_id=mid,sample=s['sample'],phenotype=s['phenotype'],cluster=s['cluster'],carrier=s['calls'][gid],partner_pattern=pattern))
assert len({r['model_id'] for r in specs}) == len(specs), 'Model identifiers must be unique'
write('model_specifications.csv',specs);write('model_sample_calls.csv',modelrows,fields=['model_id','sample','phenotype','cluster','carrier','partner_pattern'])
subprocess.run([args.rscript,'--vanilla',str(HERE/'helpers/conditional_tests.R'),str(OUT)],check=True)

results=read(OUT/'conditional_results.csv')
metadata={
 'analysis':'Exploratory exact joint co-carriage conditioning',
 'population':'Original drug-specific analysis cohorts, restricted to jointly known target and selected partner calls for each model.',
 'qualifying_set':f'{len(qualified)} gene–drug associations from the supplied PCA-only screen: BH q < 0.05 and PCA-only odds ratio > 1.',
 'partner_candidates':'Other PCA-only qualifying genes for the same drug.',
 'partner_selection':'Absolute within-PCA-cluster residual Pearson correlation >= 0.20 on jointly called binary carriage across pooled R phenotype and S phenotype isolates. The threshold is an operational exploratory choice, not a universal or validated cutoff. Negative correlations, if present, are also eligible.',
 'primary_threshold':0.20,
 'sensitivity_thresholds':[t for t in THRESHOLDS if t != .20],
 'target_test':'Exact two-sided conditional common-odds-ratio test, stratified jointly by PCA cluster and the complete selected-partner carriage pattern. With one informative stratum, Fisher exact test. This conditions on partners jointly rather than one at a time.',
 'baseline':'Same complete-call isolates, stratified by PCA cluster alone; reported alongside adjusted estimates.',
 'study_adjustment':False,
 'missingness':'Unknown calls are excluded per target-specific model, never recoded as non-carriage. No minimum carrier/called filter is imposed.',
 'inference':'Strata must contain both phenotype classes and both target-carriage classes to contribute. Rows with no informative strata are retained as not estimable. Boundary odds ratios 0/Inf with valid exact inference are retained.',
 'multiple_testing':f'BH across the fixed {len(qualified)} target-gene tests separately for each threshold and model type. Untestable tests enter as p=1 but have reported p/q missing. The legacy q_36_tests field records this follow-up family correction.',
 'interpretation':'Association conditional on the chosen correlated genes and PCA clusters. Same discovery cohort and data-informed gene/partner selection make these exploratory results, not independent replication, causal proof, or independence from all unmodelled genes.',
 'phenotype_display_mapping':{'R':'R phenotype','S':'S phenotype'},
 'sources':[args.gene_results.name,args.gene_calls.name],
}
(OUT/'analysis_metadata.json').write_text(json.dumps(metadata,indent=2))
for t in THRESHOLDS:
    rr=[r for r in results if float(r['threshold'])==t and r['method']=='pca_joint']
    significant=[r for r in rr if r['q_36_tests'] and float(r['q_36_tests'])<.05 and float(r['or'])>1]
    print(f"threshold={threshold_label(t)}: models={len(rr)}, estimable={sum(bool(r['p']) for r in rr)}, positive_q<.05={len(significant)}, positive_q<.05_with_partners={sum(int(r['n_partners'])>0 for r in significant)}")
for r in results:
    if float(r['threshold'])==.2 and r['method']=='pca_joint':print(r['drug'],r['gene'],'partners='+r['partners'],'OR='+r['or'],'q='+r['q_36_tests'],'informative_R='+r['informative_R'])
