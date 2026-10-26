#!/usr/bin/env python3
"""Compare all six Figure 4 models after exclusion of the fluconazole SDD MIC range."""
from __future__ import annotations

import argparse
import hashlib
import html
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

def project_root():
    for start in [Path.cwd(),Path(__file__).resolve().parent]:
        for candidate in [start,*start.parents]:
            if (candidate/'f/Nature_submission/Analysis/RF').is_dir():
                return candidate
    return Path.cwd()


ROOT = project_root()
BASE = ROOT/'f/Nature_submission/Analysis/RF'
STRICT = ROOT/'f/Nature_submission/Analysis/RF_SDD_sensitivity_20261007'
DESIGNS = ['stratified','clade_grouped']
INPUTS = ['variant','gene','pc']
METRICS = ['roc_auc','balanced_accuracy','precision','recall']
METRIC_LABELS = {'roc_auc':'ROC AUC','balanced_accuracy':'Balanced accuracy','precision':'Precision','recall':'Resistant-class recall'}
INPUT_LABELS = {'variant':'Variant','gene':'Gene','pc':'PC'}
DESIGN_LABELS = {'stratified':'Stratified CV','clade_grouped':'Clade-grouped CV'}
SCENARIOS = ['baseline_full','baseline_restricted','strict_retrained']
SCENARIO_LABELS = {'baseline_full':'Original training / original test cohort',
                   'baseline_restricted':'Original training / strict test cohort',
                   'strict_retrained':'Strict training / strict test cohort'}
KEYS = ['design_key','input_key','repeat','fold','fold_number']
PRED_KEYS = KEYS+['sample']


def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def write_tsv(df,path):
    df.to_csv(path,sep='\t',index=False,na_rep='',float_format='%.17g')


def metrics(frame,background=False):
    y=frame.true_label.eq('R').to_numpy()
    pred=frame.predicted_label.eq('R').to_numpy()
    score=frame.probability_R.to_numpy(float)
    tp=int((y&pred).sum());tn=int((~y&~pred).sum())
    fp=int((~y&pred).sum());fn=int((y&~pred).sum())
    recall=tp/(tp+fn) if tp+fn else (np.nan if background else 0.)
    specificity=tn/(tn+fp) if tn+fp else np.nan
    precision=tp/(tp+fp) if tp+fp else (np.nan if background else 0.)
    ba=(recall+specificity)/2 if tp+fn and tn+fp else np.nan
    auc=float(roc_auc_score(y,score)) if np.unique(y).size==2 else np.nan
    return {'n_test':len(frame),'n_R':tp+fn,'n_S':tn+fp,'tp':tp,'fn':fn,'tn':tn,'fp':fp,
            'roc_auc':auc,'balanced_accuracy':ba,'precision':precision,'recall':recall,'specificity':specificity}


def fold_metrics(predictions):
    records=[]
    for key,frame in predictions.groupby(['scenario']+KEYS,sort=False):
        record=dict(zip(['scenario']+KEYS,key))
        assert frame.threshold.nunique()==1
        record.update(metrics(frame));record['threshold']=float(frame.threshold.iloc[0])
        records.append(record)
    return pd.DataFrame(records)


def check_predictions(predictions,cohort_ids,expected_n,expected_r):
    assert not predictions.duplicated(PRED_KEYS).any(), 'Duplicate sample-fold prediction'
    assert set(zip(predictions.design_key,predictions.input_key))=={(d,m) for d in DESIGNS for m in INPUTS}
    assert set(predictions.true_label)=={'R','S'} and set(predictions.predicted_label)<= {'R','S'}
    assert np.isfinite(predictions.probability_R).all() and predictions.probability_R.between(0,1).all()
    assert np.isfinite(predictions.threshold).all() and predictions.threshold.between(0,1).all()
    expected=np.where(predictions.probability_R.ge(predictions.threshold),'R','S')
    assert np.array_equal(expected,predictions.predicted_label), 'Stored labels differ from stored decision thresholds'
    assert ((predictions.fold_number-1)//5+1).eq(predictions['repeat']).all()
    assert ((predictions.fold_number-1)%5+1).eq(predictions.fold).all()
    for (design,model),frame in predictions.groupby(['design_key','input_key']):
        assert set(frame.fold_number)==set(range(1,26))
        for rep,part in frame.groupby('repeat'):
            assert len(part)==expected_n and set(part['sample'])==cohort_ids
            assert part['sample'].is_unique and int(part.true_label.eq('R').sum())==expected_r
    for design in DESIGNS:
        ref=predictions[predictions.design_key.eq(design)&predictions.input_key.eq('variant')][['repeat','fold','sample','true_label']]
        for model in ['gene','pc']:
            other=predictions[predictions.design_key.eq(design)&predictions.input_key.eq(model)][ref.columns]
            pd.testing.assert_frame_equal(ref.sort_values(['repeat','fold','sample']).reset_index(drop=True),other.sort_values(['repeat','fold','sample']).reset_index(drop=True))


def load_baseline(base):
    pred_path=base/'results/combined_cv_fold_predictions.tsv'
    metadata_path=base/'inputs/labelled_metadata.tsv'
    pred=pd.read_csv(pred_path,sep='\t')
    meta=pd.read_csv(metadata_path,sep='\t')
    samples=pred[['sample','true_label']].drop_duplicates()
    assert samples['sample'].is_unique and len(samples)==594
    cohort=samples.merge(meta[['sample','Fluconazole','cluster']],on='sample',validate='one_to_one')
    cohort['MIC']=pd.to_numeric(cohort.Fluconazole,errors='raise')
    assert cohort.MIC.gt(0).all()
    assert cohort.true_label.eq(np.where(cohort.MIC.gt(2),'R','S')).all()
    cohort['excluded_SDD']=cohort.MIC.gt(2)&cohort.MIC.le(4)
    assert int(cohort.excluded_SDD.sum())==31
    assert set(cohort.loc[cohort.excluded_SDD,'MIC'])=={3.,4.}
    strict=cohort.loc[~cohort.excluded_SDD].copy()
    assert len(strict)==563 and int(strict.true_label.eq('R').sum())==204
    assert strict.true_label.eq(np.where(strict.MIC.gt(4),'R','S')).all()
    groups=pd.read_csv(base/'source_data/grouped_group_summary.tsv',sep='\t')
    groups['cluster']=groups.cluster.astype(int)
    cohort['cluster']=cohort.cluster.astype(int)
    cohort=cohort.merge(groups[['cluster','group']],on='cluster',validate='many_to_one')
    assert len(cohort)==594 and cohort['group'].nunique()==7
    check_predictions(pred,set(cohort['sample']),594,235)
    pred['scenario']='baseline_full'
    strict_ids=set(strict['sample'])
    restricted=pred[pred['sample'].isin(strict_ids)].copy()
    restricted['scenario']='baseline_restricted'
    check_predictions(restricted,strict_ids,563,204)
    return pred,restricted,cohort,[pred_path,metadata_path,base/'source_data/grouped_group_summary.tsv']


def load_strict(strict_root,base,strict_ids,restricted):
    frames=[];sources=[];run_audits=[]
    for design in DESIGNS:
        for model in INPUTS:
            for number in range(1,26):
                folder=strict_root/'results'/f'{design}_{model}'/'folds'/f'{number:02d}'
                pp=folder/'cv_fold_predictions.tsv'
                mp=folder/'cv_fold_metrics.tsv'
                rp=folder/'run_info.json'
                if not all(p.exists() for p in [pp,mp,rp]):
                    raise RuntimeError(f'Complete fold outputs required: {folder}')
                frame=pd.read_csv(pp,sep='\t')
                assert frame.input_key.eq(model).all() and frame.design_key.eq(design).all()
                assert frame.fold_number.eq(number).all()
                frame['scenario']='strict_retrained';frames.append(frame)
                report=pd.read_csv(mp,sep='\t')
                assert len(report)==1
                calc=metrics(frame)
                for metric in METRICS+['tp','tn','fp','fn','n_test']:
                    if metric in report:
                        assert np.isclose(float(report.iloc[0][metric]),calc[metric],rtol=1e-10,atol=1e-12,equal_nan=True),(pp,metric)
                info=json.loads(rp.read_text())
                original_path=base/'results'/f'{design}_{model}'/'folds'/f'{number:02d}'/'run_info.json'
                original=json.loads(original_path.read_text())
                fields=['input_type','design','outer_folds','outer_repeats','inner_folds','seed','trees','top_k','maximum_training_missingness','rf_grid','pc_grid','feature_selection_scope','versions','n_input_features']
                for field in fields:
                    assert info.get(field)==original.get(field),(rp,field,info.get(field),original.get(field))
                assert info['n_samples']==563 and info['n_resistant']==204 and info['n_susceptible']==359
                assert info.get('outer_partition_source') and info.get('inner_partition_source'),f'Missing partition provenance: {rp}'
                run_audits.append({'design_key':design,'input_key':model,'fold_number':number,'seed':info['seed'],
                                   'outer_partition_source':info['outer_partition_source'],'inner_partition_source':info['inner_partition_source'],
                                   'core_settings_match_baseline':True,'run_info':str(rp)})
                sources.extend([pp,mp,rp,original_path])
    pred=pd.concat(frames,ignore_index=True)
    check_predictions(pred,strict_ids,563,204)
    cols=PRED_KEYS+['true_label']
    left=restricted[cols].sort_values(PRED_KEYS).reset_index(drop=True)
    right=pred[cols].sort_values(PRED_KEYS).reset_index(drop=True)
    pd.testing.assert_frame_equal(left,right)
    return pred,sources,run_audits


def metric_summaries(folds):
    long=folds.melt(id_vars=['scenario']+KEYS,value_vars=METRICS,var_name='metric',value_name='value')
    summary=long.groupby(['scenario','design_key','input_key','metric'],sort=False)['value'].agg(mean='mean',sd='std',n_folds='count',minimum='min',maximum='max').reset_index()
    assert summary.n_folds.eq(25).all()
    pairs=long.pivot(index=KEYS+['metric'],columns='scenario',values='value').reset_index()
    pairs['delta_retrained_minus_restricted']=pairs.strict_retrained-pairs.baseline_restricted
    pairs['delta_restricted_minus_full']=pairs.baseline_restricted-pairs.baseline_full
    pairs['delta_retrained_minus_full']=pairs.strict_retrained-pairs.baseline_full
    delta_summary=pairs.groupby(['design_key','input_key','metric'],sort=False)['delta_retrained_minus_restricted'].agg(mean_delta='mean',sd_delta='std',n_folds='count',minimum_delta='min',maximum_delta='max').reset_index()
    repeat_deltas=pairs.groupby(['design_key','input_key','repeat','metric'],sort=False)[['baseline_full','baseline_restricted','strict_retrained','delta_retrained_minus_restricted','delta_restricted_minus_full','delta_retrained_minus_full']].mean().reset_index()
    consistency=repeat_deltas.groupby(['design_key','input_key','metric'],sort=False)['delta_retrained_minus_restricted'].agg(n_repeats='count',positive_repeats=lambda s:int(s.gt(1e-12).sum()),negative_repeats=lambda s:int(s.lt(-1e-12).sum()),zero_repeats=lambda s:int(s.abs().le(1e-12).sum()),min_repeat_delta='min',max_repeat_delta='max').reset_index()
    delta_summary=delta_summary.merge(consistency,on=['design_key','input_key','metric'],validate='one_to_one')
    return long,summary,pairs,delta_summary,repeat_deltas


def background_summaries(predictions,cohort):
    grouped=predictions[predictions.design_key.eq('clade_grouped')].merge(cohort[['sample','cluster','group']],on='sample',validate='many_to_one')
    rows=[]
    for key,frame in grouped.groupby(['scenario','input_key','repeat','cluster','group'],sort=False):
        assert frame.fold_number.nunique()==1 and frame['sample'].is_unique
        row=dict(zip(['scenario','input_key','repeat','cluster','group'],key))
        row.update(metrics(frame,background=True));rows.append(row)
    detail=pd.DataFrame(rows)
    long=detail.melt(id_vars=['scenario','input_key','repeat','cluster','group','n_test','n_R','n_S'],value_vars=METRICS+['specificity'],var_name='metric',value_name='value')
    summary=long.groupby(['scenario','input_key','cluster','group','metric','n_test','n_R','n_S'],sort=False)['value'].agg(mean='mean',sd='std',n_defined_repeats='count',minimum='min',maximum='max').reset_index()
    paired=long.pivot(index=['input_key','repeat','cluster','group','metric'],columns='scenario',values='value').reset_index()
    paired['delta_retrained_minus_restricted']=paired.strict_retrained-paired.baseline_restricted
    return detail,summary,paired


def draw_figures(summary,pairs,delta_summary,repeat_deltas,out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':12,'axes.labelsize':10,
                         'pdf.fonttype':42,'ps.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
    colors={'baseline_full':'#7A7A7A','baseline_restricted':'#D98E04','strict_retrained':'#0072B2'}
    markers={'baseline_full':'o','baseline_restricted':'s','strict_retrained':'D'}
    model_order=[(d,m) for d in DESIGNS for m in INPUTS]
    ticklabels=[INPUT_LABELS[m]+'\n'+('Stratified' if d=='stratified' else 'Grouped') for d,m in model_order]
    fig,axes=plt.subplots(2,2,figsize=(12.5,8.5),sharex=True)
    fig.subplots_adjust(left=.08,right=.985,top=.82,bottom=.14,wspace=.20,hspace=.20)
    for j,(ax,metric) in enumerate(zip(axes.flat,METRICS)):
        upper=1.02;lower=0.
        for si,scenario in enumerate(SCENARIOS):
            frame=summary[summary.metric.eq(metric)&summary.scenario.eq(scenario)].set_index(['design_key','input_key']).loc[model_order]
            xx=np.arange(6)+(si-1)*.19
            ax.errorbar(xx,frame['mean'],yerr=frame.sd,fmt=markers[scenario],color=colors[scenario],mfc='white' if scenario=='baseline_full' else colors[scenario],markersize=5.5,capsize=2.5,lw=1.1,label=SCENARIO_LABELS[scenario])
            upper=max(upper,float((frame['mean']+frame.sd).max())+.035)
            lower=min(lower,float((frame['mean']-frame.sd).min())-.025)
        ax.axvline(2.5,color='#BBBBBB',lw=.8,ls=':')
        ax.set_title(chr(97+j)+'. '+METRIC_LABELS[metric],loc='left',fontweight='bold')
        ax.set_ylim(lower,upper);ax.set_xlim(-.6,5.6);ax.set_xticks(range(6),ticklabels)
        ax.grid(axis='y',color='#E5E5E5',lw=.6);ax.set_axisbelow(True)
        ax.set_ylabel('Mean across 25 outer folds')
    handles=[Line2D([],[],color=colors[s],marker=markers[s],mfc='white' if s=='baseline_full' else colors[s],ls='none',label=SCENARIO_LABELS[s]) for s in SCENARIOS]
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,.96),ncol=1,frameon=False)
    fig.suptitle('Fluconazole SDD-range exclusion: prediction performance',fontsize=15,fontweight='bold',y=.995)
    fig.supxlabel('Original cohort: 594 isolates (235 originally labelled R, 359 S). Strict cohort: 563 isolates (204 R, 359 S).\nPoints: fold means; error bars: fold SD.',fontsize=9,y=.02)
    for ext in ['png','pdf']:fig.savefig(out/f'performance_comparison.{ext}',dpi=300,bbox_inches='tight',facecolor='white')
    plt.close(fig)

    fig,axes=plt.subplots(2,2,figsize=(12.5,8.3),sharex=True,layout='constrained')
    for j,(ax,metric) in enumerate(zip(axes.flat,METRICS)):
        sub=delta_summary[delta_summary.metric.eq(metric)].set_index(['design_key','input_key']).loc[model_order]
        ax.axhline(0,color='#333333',lw=.9)
        ax.errorbar(np.arange(6),sub.mean_delta,yerr=sub.sd_delta,fmt='D',color='#0072B2',markersize=6,capsize=3,lw=1.2,zorder=4)
        for i,(design,model) in enumerate(model_order):
            vals=repeat_deltas[repeat_deltas.metric.eq(metric)&repeat_deltas.design_key.eq(design)&repeat_deltas.input_key.eq(model)].sort_values('repeat').delta_retrained_minus_restricted.to_numpy()
            ax.scatter(i+np.linspace(-.12,.12,len(vals)),vals,s=18,c='#333333',alpha=.65,zorder=5)
        extent=max(.03,float(np.max(np.abs(np.r_[sub.mean_delta-sub.sd_delta,sub.mean_delta+sub.sd_delta])))*1.15)
        ax.set_ylim(-extent,extent);ax.set_xlim(-.6,5.6);ax.set_xticks(range(6),ticklabels)
        ax.axvline(2.5,color='#BBBBBB',lw=.8,ls=':');ax.grid(axis='y',color='#E5E5E5',lw=.6);ax.set_axisbelow(True)
        ax.set_title(chr(97+j)+'. '+METRIC_LABELS[metric],loc='left',fontweight='bold')
        ax.set_ylabel('Retrained − original training\n(same strict test samples)')
    fig.suptitle('Effect of retraining on matched strict test folds',fontsize=15,fontweight='bold')
    fig.supxlabel('Blue diamonds: mean paired change ± SD across 25 outer folds.\nSmall grey circles: mean paired change within each of five repeats.',fontsize=9)
    for ext in ['png','pdf']:fig.savefig(out/f'paired_retraining_change.{ext}',dpi=300,bbox_inches='tight',facecolor='white')
    plt.close(fig)


def write_report(summary,deltas,cohort,out,audit):
    indexed=summary.set_index(['scenario','design_key','input_key','metric'])
    dindex=deltas.set_index(['design_key','input_key','metric'])
    lines=['FLUCONAZOLE SDD-RANGE EXCLUSION SENSITIVITY ANALYSIS','',
           'The strict cohort contains 563 isolates: 204 resistant and 359 susceptible. Thirty-one isolates with 2 < MIC ≤ 4 µg/ml were removed from the original 594-isolate cohort, which contained 235 originally labelled resistant and 359 susceptible isolates. The retained resistant isolates all have recorded MICs ≥8 µg/ml. The recorded excluded MICs are 3 and 4 µg/ml; 3 rounds up to the next twofold-dilution category. Twenty-nine excluded isolates belong to clade 1, one to clade 2 and one to clade 3.',
           'This retrospective analysis applies the CLSI susceptible-dose-dependent (SDD) MIC range to harmonized numeric measurements from published and laboratory sources. The exclusion is defined from the recorded MIC, with source susceptibility methods and labels retaining their original provenance.',
           'Threshold reference: CLSI M27M44S, 3rd edition, Table 1 and the next-higher-category rule for intermediate dilution values: https://clsi.org/media/osthhxax/m27m44sed3e_sample.pdf','',
           'Procedure',
           'All six Figure 4 model combinations were retrained: variant, gene and PC inputs under stratified and clade-grouped validation. The original outer and inner partitions were retained after removing the SDD-range isolates. Rare-variant eligibility remains defined in the original 871-isolate fluconazole cohort, and the common LD-pruned PC input panel remains fixed. Feature selection, PC fitting, hyperparameter retuning and threshold selection used strict training data. Each model has five repeats of five outer folds.','',
           'Three comparisons separate test-cohort composition from training changes: (1) original training and all 594 test isolates; (2) original out-of-fold predictions evaluated on the 563 strict isolates; and (3) strict retraining evaluated on the identical strict test folds. Primary estimates are unweighted means across 25 outer-fold metrics. SD describes repeated-validation variability. Retraining differences are paired by model and outer fold, and repeat-level summaries average the five paired fold differences within each repeat. These folds share training observations, so changes are interpreted descriptively.','',
           'Findings']
    for design in DESIGNS:
        parts=[]
        for model in INPUTS:
            x=indexed.loc[('strict_retrained',design,model,'roc_auc')]
            parts.append(f'{INPUT_LABELS[model].lower()} {x["mean"]:.3f} ± {x.sd:.3f}')
        lines.append(f'{DESIGN_LABELS[design]} strict-retrained ROC AUC (mean ± fold SD): '+', '.join(parts)+'.')
    for design in DESIGNS:
        lines.append('')
        lines.append(DESIGN_LABELS[design]+' matched retraining changes:')
        for model in INPUTS:
            parts=[]
            for metric in METRICS:
                x=dindex.loc[(design,model,metric)]
                parts.append(f'{METRIC_LABELS[metric]} {x.mean_delta:+.3f} (positive/negative/zero repeats: {int(x.positive_repeats)}/{int(x.negative_repeats)}/{int(x.zero_repeats)})')
            lines.append(INPUT_LABELS[model]+': '+'; '.join(parts)+'.')
    lines += ['', 'Model performance (mean ± fold SD)', 'Validation\tInput\tMetric\tOriginal 594\tOriginal evaluated on strict 563\tRetrained strict 563']
    for design in DESIGNS:
        for model in INPUTS:
            for metric in METRICS:
                vals=[]
                for scenario in SCENARIOS:
                    x=indexed.loc[(scenario,design,model,metric)]
                    vals.append(f'{x["mean"]:.4f} ± {x.sd:.4f}')
                lines.append('\t'.join([DESIGN_LABELS[design],INPUT_LABELS[model],METRIC_LABELS[metric]]+vals))
    if audit.get('original_reproduction_control'):
        control=audit['original_reproduction_control']['evidence']
        lines += ['',f'An original-cohort control reproduced the archived fold with {control["predictions"]["label_differences"]} label differences and a maximum absolute probability difference of {control["predictions"]["max_absolute_probability_difference"]:.2g}. Selected hyperparameters, feature indices and inner tuning results matched. Numerical library versions matched the original run.']
    lines += ['', 'Interpretation of denominators',
              'Fold-level precision is zero when no resistant isolate is predicted, matching the Figure 4 convention. Background-specific precision is blank when no positive predictions occur; recall requires resistant isolates, and balanced accuracy and ROC AUC require both phenotype classes. Background summaries include the number of defined repeats. Changes in resistant prevalence can affect precision, so paired strict-test comparisons hold the evaluation population fixed.','',
              'Verification',f'All {audit["strict_completed_folds"]} strict fold outputs were present. Test-sample identities and labels matched the restricted baseline in every corresponding model and fold. Each strict isolate appeared once per repeat in each model. Recalculated stored metrics and original baseline summaries agreed within numerical precision. Model settings and seed records matched their original fold records.','',
              'Files','performance_comparison.png / .pdf','paired_retraining_change.png / .pdf','performance_summary.tsv','paired_fold_changes.tsv','paired_change_summary.tsv','repeat_mean_changes.tsv','grouped_background_summary.tsv','grouped_background_paired_changes.tsv','source_audit.json']
    (out/'comparison_report.txt').write_text('\n'.join(lines)+'\n')
    table_start=lines.index('Model performance (mean ± fold SD)')
    table_end=table_start+2+len(DESIGNS)*len(INPUTS)*len(METRICS)
    table_rows=[line.split('\t') for line in lines[table_start+1:table_end] if '\t' in line]
    table='<table><thead><tr>'+''.join('<th>'+html.escape(cell)+'</th>' for cell in table_rows[0])+'</tr></thead><tbody>'
    table+=''.join('<tr>'+''.join('<td>'+html.escape(cell)+'</td>' for cell in row)+'</tr>' for row in table_rows[1:])+'</tbody></table>'
    paragraphs=[]
    headings={'Procedure','Findings','Interpretation of denominators','Verification','Files'}
    for i,line in enumerate(lines[1:],1):
        if not line or table_start<=i<table_end:continue
        if line in headings:paragraphs.append('<h2>'+html.escape(line)+'</h2>')
        else:paragraphs.append('<p>'+html.escape(line)+'</p>')
    body=''.join(paragraphs)
    body=body.replace('<h2>Interpretation of denominators</h2>','<h2>Model performance (mean ± fold SD)</h2>'+table+'<h2>Interpretation of denominators</h2>')
    content='<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>Fluconazole SDD-range sensitivity analysis</title><style>body{max-width:1120px;margin:36px auto;padding:0 22px;color:#192c3a;font:16px/1.6 system-ui,sans-serif}h1{line-height:1.2}h2{font-size:1.2rem;margin-top:2rem}img{width:100%;height:auto}table{border-collapse:collapse;width:100%;font-size:13px}th,td{text-align:left;padding:8px;border-bottom:1px solid #d7e0e6}th{background:#eaf2f7}tr:nth-child(even){background:#f8fafb}a{color:#0072b2}p{overflow-wrap:anywhere}</style><h1>Fluconazole SDD-range sensitivity analysis</h1>'+body+'<h2>Performance across the three conditions</h2><a href="performance_comparison.pdf"><img src="performance_comparison.png" alt="Mean ROC AUC, balanced accuracy, precision and resistant-class recall across all six model combinations for the original cohort, restricted original predictions and strict retraining."></a><h2>Paired effect of retraining</h2><a href="paired_retraining_change.pdf"><img src="paired_retraining_change.png" alt="Paired changes in performance on identical strict test folds. Blue diamonds show the mean plus or minus the fold standard deviation; grey circles show repeat means."></a></html>'
    (out/'comparison_report.html').write_text(content)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root','--baseline-dir',dest='baseline_dir',type=Path,default=BASE,help='Original Figure 4 RF analysis directory (inputs, results and source_data).')
    parser.add_argument('--analysis-root','--strict-dir',dest='strict_dir',type=Path,default=STRICT,help='Directory containing the strict retraining results.')
    parser.add_argument('--outdir',type=Path)
    parser.add_argument('--reproduction-control',type=Path,help='Optional original-cohort reproduction comparison JSON to include in the audit.')
    parser.add_argument('--validate-baseline-only',action='store_true')
    args=parser.parse_args()
    out=args.outdir or args.strict_dir/'comparison'
    baseline,restricted,cohort,sources=load_baseline(args.baseline_dir)
    base_metrics=fold_metrics(baseline)
    stored=pd.read_csv(args.baseline_dir/'results/combined_cv_fold_metrics.tsv',sep='\t')
    joined=base_metrics.merge(stored,on=KEYS,suffixes=('_recomputed','_stored'),validate='one_to_one')
    for metric in METRICS:
        assert np.allclose(joined[metric+'_recomputed'],joined[metric+'_stored'],rtol=1e-10,atol=1e-12,equal_nan=True),metric
    if args.validate_baseline_only:
        print(json.dumps({'baseline_samples':594,'baseline_R':235,'baseline_S':359,'strict_samples':563,'strict_R':204,'strict_S':359,'excluded':31,'baseline_fold_metrics_match':True,'groups_excluded':cohort[cohort.excluded_SDD].groupby('group').size().to_dict()},indent=2))
        return
    strict_ids=set(cohort.loc[~cohort.excluded_SDD,'sample'])
    strict,strict_sources,run_audits=load_strict(args.strict_dir,args.baseline_dir,strict_ids,restricted)
    predictions=pd.concat([baseline,restricted,strict],ignore_index=True)
    folds=fold_metrics(predictions)
    long,summary,pairs,delta_summary,repeat_deltas=metric_summaries(folds)
    background_detail,background_summary,background_pairs=background_summaries(predictions,cohort)
    out.mkdir(parents=True,exist_ok=True)
    sources += strict_sources + [args.baseline_dir/'results/combined_cv_fold_metrics.tsv']
    audit={'strict_completed_folds':150,'scenario_counts':{'baseline_full':594,'baseline_restricted':563,'strict_retrained':563},
           'strict_class_counts':{'R':204,'S':359},'exclusion_rule':'2 < fluconazole MIC <= 4 ug/ml',
           'exact_test_membership_and_labels_match':True,'baseline_recomputed_metrics_match':True,
           'paired_comparison':'strict_retrained minus baseline_restricted on identical sample-fold memberships',
           'primary_summary':'Unweighted mean and sample SD (ddof=1) of 25 outer-fold metrics',
           'repeat_summary':'Unweighted mean of the five matched fold differences in each repeat',
           'run_info_checks':run_audits,'sources':{str(p):sha(p) for p in dict.fromkeys(sources)}}
    if args.reproduction_control:
        evidence=json.loads(args.reproduction_control.read_text())
        assert evidence['predictions']['label_differences']==0
        assert evidence['predictions']['max_absolute_probability_difference']<1e-12
        assert evidence['predictions']['threshold_max_absolute_difference']==0
        assert evidence['features']['ordered_indices_equal']
        assert evidence['selected_hyperparameters_by_fold.tsv']['exact_equal']
        assert evidence['inner_cv_tuning_results.tsv']['exact_equal']
        assert evidence['versions']==evidence['original_versions']
        audit['original_reproduction_control']={'path':str(args.reproduction_control),'sha256':sha(args.reproduction_control),'evidence':evidence}
    write_tsv(cohort[['sample','MIC','true_label','cluster','group','excluded_SDD']],out/'cohort_comparison.tsv')
    write_tsv(predictions[['scenario']+PRED_KEYS+['true_label','predicted_label','probability_R','threshold']],out/'comparison_predictions.tsv.gz')
    write_tsv(folds,out/'fold_metrics.tsv')
    write_tsv(summary,out/'performance_summary.tsv')
    write_tsv(pairs,out/'paired_fold_changes.tsv')
    write_tsv(delta_summary,out/'paired_change_summary.tsv')
    write_tsv(repeat_deltas,out/'repeat_mean_changes.tsv')
    write_tsv(background_detail,out/'grouped_background_metrics.tsv')
    write_tsv(background_summary,out/'grouped_background_summary.tsv')
    write_tsv(background_pairs,out/'grouped_background_paired_changes.tsv')
    draw_figures(summary,pairs,delta_summary,repeat_deltas,out)
    write_report(summary,delta_summary,cohort,out,audit)
    (out/'source_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(json.dumps({'outdir':str(out),'prediction_rows':len(predictions),'fold_metric_rows':len(folds),'summary_rows':len(summary),'paired_metric_rows':len(pairs),'strict_completed_folds':150},indent=2))


if __name__=='__main__':
    main()
