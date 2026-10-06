#!/usr/bin/env python3
"""Validate held-out predictions and export performance and partition summaries."""
from __future__ import annotations
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

SPEC=importlib.util.spec_from_file_location('prediction_plot',Path(__file__).with_name('04_plot_prediction.py'))
PLOT=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PLOT)
METRICS=['balanced_accuracy','precision','recall','specificity','roc_auc']


def evaluate(frame,fold_level=False):
    y=frame.true_label.eq('R').to_numpy()
    pred=frame.predicted_label.eq('R').to_numpy()
    tp=int((y&pred).sum());tn=int((~y&~pred).sum())
    fp=int((~y&pred).sum());fn=int((y&~pred).sum())
    recall=tp/(tp+fn) if tp+fn else np.nan
    specificity=tn/(tn+fp) if tn+fp else np.nan
    precision=tp/(tp+fp) if tp+fp else (0.0 if fold_level else np.nan)
    return dict(n=len(y),n_R=int(y.sum()),n_S=int((~y).sum()),tp=tp,tn=tn,fp=fp,fn=fn,
                balanced_accuracy=(recall+specificity)/2,precision=precision,recall=recall,
                specificity=specificity,roc_auc=roc_auc_score(y,frame.probability_R) if y.any() and not y.all() else np.nan)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True)
    args=parser.parse_args()
    config=json.loads((args.root/'run_configuration.json').read_text())
    if config.get('status')!='complete':raise ValueError('All model tasks must be complete before summary')
    out=args.root/'source_data'
    out.mkdir(exist_ok=True)
    metadata=pd.read_csv(config['cohort_metadata'],sep='\t')
    metrics,predictions=PLOT.load_results(args.root/'results')
    expected_n=len(metadata)
    reference={}
    audit=[]
    for (design,kind),part in predictions.groupby(['design','input_key']):
        if len(part)!=expected_n*5 or part.groupby('repeat').size().ne(expected_n).any():
            raise ValueError('Incomplete repeated prediction coverage')
        assignments=part[['repeat','fold','fold_number','sample','true_label']].sort_values(['repeat','sample']).reset_index(drop=True)
        if design not in reference:reference[design]=assignments
        else:pd.testing.assert_frame_equal(reference[design],assignments)
        for number,fold in part.groupby('fold_number'):
            actual=evaluate(fold,fold_level=True)
            stored=metrics.loc[metrics.design.eq(design)&metrics.input_key.eq(kind)&metrics.fold_number.eq(number)].iloc[0]
            for key in ['balanced_accuracy','precision','recall','roc_auc','tp','tn','fp','fn']:
                if not np.isclose(actual[key],stored[key],atol=1e-12):raise ValueError(f'Metric mismatch {design}/{kind}/{number}/{key}')
            if not np.array_equal(fold.probability_R.ge(fold.threshold),fold.predicted_label.eq('R')):
                raise ValueError('Stored labels differ from thresholded probabilities')
            audit.append(dict(design=design,input_key=kind,fold_number=number,n_test=len(fold),metrics_recomputed=True))
    metrics.to_csv(out/'fold_metrics.tsv',sep='\t',index=False)
    PLOT.roc_summary(predictions).to_csv(out/'mean_roc_coordinates.tsv',sep='\t',index=False)
    summary=metrics.groupby(['design','input_key'])[['roc_auc','balanced_accuracy','precision','recall']].agg(['mean','std','median','min','max'])
    summary.columns=['_'.join(c) for c in summary.columns]
    summary.reset_index().to_csv(out/'performance_summary.tsv',sep='\t',index=False)
    boxrows=[]
    for (design,kind),part in metrics.groupby(['design','input_key']):
        for metric in ['balanced_accuracy','precision','recall']:
            a=part[metric].to_numpy();q1,median,q3=np.percentile(a,[25,50,75]);iqr=q3-q1
            boxrows.append(dict(design=design,input_key=kind,metric=metric,n=len(a),q1=q1,median=median,q3=q3,
                                whisker_low=a[a>=q1-1.5*iqr].min(),whisker_high=a[a<=q3+1.5*iqr].max()))
    pd.DataFrame(boxrows).to_csv(out/'boxplot_summary.tsv',sep='\t',index=False)
    cohort=metadata[['sample','cluster','rbin_flu','background']].rename(columns={'rbin_flu':'true_label','background':'group'})
    groups=cohort.groupby(['cluster','group'],as_index=False).agg(n=('sample','size'),n_R=('true_label',lambda x:x.eq('R').sum()),n_S=('true_label',lambda x:x.eq('S').sum()))
    groups.to_csv(out/'grouped_group_summary.tsv',sep='\t',index=False)
    cohort.to_csv(out/'grouped_cohort_metadata.tsv',sep='\t',index=False)
    assigned=reference['clade_grouped'].merge(cohort[['sample','cluster','group']],on='sample',validate='many_to_one')
    if assigned.groupby(['repeat','cluster']).fold.nunique().ne(1).any():raise ValueError('A genetic group crosses test folds')
    composition,membership,sample_membership=[],[],[]
    all_samples=set(cohort['sample'])
    for number,rows in assigned.groupby('fold_number'):
        repeat=int(rows.repeat.iloc[0]);fold=int(rows.fold.iloc[0]);test=set(rows['sample'])
        train=cohort.loc[~cohort['sample'].isin(test)]
        train_groups=set(train.cluster);test_groups=set(rows.cluster)
        if train_groups&test_groups:raise ValueError('Train/test genetic-group overlap')
        for g in groups.itertuples():
            membership.append(dict(repeat=repeat,fold=fold,fold_number=number,cluster=g.cluster,group=g.group,
                                   partition='Test' if g.cluster in test_groups else 'Training'))
        for s in cohort.itertuples():
            sample_membership.append(dict(repeat=repeat,fold=fold,fold_number=number,sample=s.sample,
                                          cluster=s.cluster,group=s.group,true_label=s.true_label,
                                          partition='Test' if s.sample in test else 'Training'))
        composition.append(dict(repeat=repeat,fold=fold,fold_number=number,n_train=len(train),n_test=len(rows),
                                train_R=int(train.true_label.eq('R').sum()),train_S=int(train.true_label.eq('S').sum()),
                                test_R=int(rows.true_label.eq('R').sum()),test_S=int(rows.true_label.eq('S').sum()),
                                test_R_fraction=rows.true_label.eq('R').mean(),
                                training_clusters=','.join(map(str,sorted(train_groups))),test_clusters=','.join(map(str,sorted(test_groups))),
                                training_groups='; '.join(sorted(train.group.unique())),test_groups='; '.join(sorted(rows.group.unique())),
                                clade1_in_training='Clade 1' in set(train.group)))
    composition=pd.DataFrame(composition)
    composition.to_csv(out/'grouped_fold_composition.tsv',sep='\t',index=False)
    pd.DataFrame(membership).to_csv(out/'grouped_group_partition_membership.tsv',sep='\t',index=False)
    pd.DataFrame(sample_membership).to_csv(out/'grouped_sample_partition_membership.tsv.gz',sep='\t',index=False)
    joined=predictions.loc[predictions.design.eq('clade_grouped')].merge(cohort[['sample','cluster','group']],on='sample',validate='many_to_one')
    joined.to_csv(out/'grouped_heldout_predictions_with_groups.tsv.gz',sep='\t',index=False)
    group_metrics=[]
    for (kind,repeat,cluster),rows in joined.groupby(['input_key','repeat','cluster']):
        group_metrics.append(dict(model=kind,repeat=repeat,fold=int(rows.fold.iloc[0]),fold_number=int(rows.fold_number.iloc[0]),
                                  cluster=cluster,group=rows.group.iloc[0],**evaluate(rows)))
    group_metrics=pd.DataFrame(group_metrics)
    group_metrics.to_csv(out/'grouped_heldout_group_metrics_by_repeat.tsv',sep='\t',index=False)
    bygroup=group_metrics.groupby(['cluster','group','model'])[METRICS].agg(['mean','min','max','count'])
    bygroup.columns=['_'.join(c) for c in bygroup.columns]
    bygroup.reset_index().merge(groups,on=['cluster','group'],validate='many_to_one').to_csv(out/'grouped_heldout_group_metric_summary.tsv',sep='\t',index=False)
    grouped=metrics.loc[metrics.design.eq('clade_grouped')].rename(columns={'input_key':'model'}).merge(composition,on=['repeat','fold','fold_number','n_train','n_test'],validate='many_to_one')
    grouped.to_csv(out/'grouped_fold_metrics_with_composition.tsv',sep='\t',index=False)
    conditional=grouped.groupby(['model','clade1_in_training'])[['balanced_accuracy','precision','recall','roc_auc','n_train','train_R','n_test','test_R']].agg(['mean','min','max','count'])
    conditional.columns=['_'.join(c) for c in conditional.columns]
    conditional.reset_index().to_csv(out/'grouped_clade1_training_presence_descriptive_summary.tsv',sep='\t',index=False)
    pd.DataFrame(audit).to_csv(out/'prediction_metric_audit.tsv',sep='\t',index=False)
    print(f'Validated {len(audit)} outer models and exported grouped performance summaries')


if __name__=='__main__':main()
