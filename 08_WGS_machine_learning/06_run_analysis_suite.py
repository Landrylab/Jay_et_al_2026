#!/usr/bin/env python3
"""Checkpoint repeated prediction and within-group random forests by outer fold."""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import time
import traceback

import numpy as np
import pandas as pd

from rf_common import aligned_metadata, load_cache, repeated_outer_splits, inner_splits

SPEC = importlib.util.spec_from_file_location("rf_driver", Path(__file__).with_name("02_run_random_forest.py"))
DRIVER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(DRIVER)
LABELS = {"variant": "Variant based", "gene": "Gene-level variant based", "pc": "Phylogeny based"}


def prepare(args):
    manifest = pd.read_csv(args.inputs, sep="\t")
    caches = {r.input_type: str(Path(r.cache_dir).resolve()) for r in manifest.itertuples()}
    metadata = pd.read_csv(args.metadata, sep="\t", dtype={"sample":str})
    cohort = pd.read_csv(args.cohort_metadata, sep="\t", dtype={"sample":str})
    groups = pd.read_csv(args.groups, sep="\t", dtype={"cluster":str})
    for kind, cache in caches.items():
        _, _, samples = load_cache(cache)
        if kind == "variant": full_samples = samples
        elif samples != full_samples: raise ValueError("All input sample orders must agree")
    sample_order = {s:i for i,s in enumerate(full_samples)}
    if not set(cohort['sample']).issubset(sample_order):
        raise ValueError("The analysis cohort includes samples absent from the input cache")
    cohort_samples = sorted(cohort['sample'].tolist(), key=sample_order.get)
    if cohort_samples != cohort['sample'].tolist():
        raise ValueError("Cohort metadata must preserve the canonical labelled-sample order")
    tasks = []
    for kind in ["gene", "pc", "variant"]:
        for design in ["stratified", "clade_grouped"]:
            tasks.append(dict(task=f"{design}_{kind}", scope="prediction", input_type=kind,
                              design=design, n_splits=5, n_repeats=5, samples=cohort_samples))
    for group in groups.to_dict("records"):
        gmeta=metadata.loc[metadata.cluster.astype(str).eq(str(group["cluster"]))]
        gsamples=[s for s in full_samples if s in set(gmeta['sample'])]
        for kind in ["gene", "variant"]:
            tasks.append(dict(task=f"{group['group_key']}/{kind}",scope="within_group", input_type=kind,
                              design="stratified",n_splits=int(group["n_splits"]), n_repeats=5,
                              samples=gsamples,group_key=group["group_key"],group=group["group_label"],
                              cluster=str(group["cluster"])))
    args.outdir.mkdir(parents=True, exist_ok=True)
    config = dict(input_caches=caches, metadata=str(args.metadata.resolve()),
                  cohort_metadata=str(args.cohort_metadata.resolve()), groups=str(args.groups.resolve()),
                  seed=42, trees=500, inner_folds=3, top_k=10000, max_missing=0.10,
                  tasks=tasks, created_utc=datetime.now(timezone.utc).isoformat(),
                  initial_gene_frequency_filter=None,
                  gene_variant_frequency_scope="full original binary-phenotyped cohort before sample subsetting",
                  model_settings="Seeded partitions, ordered candidate grids, training-based thresholds and preprocessing scopes")
    (args.outdir/"run_configuration.json").write_text(json.dumps(config, indent=2)+"\n")
    outer,inner=[] , []
    for task in tasks:
        y, labels=aligned_metadata(metadata,task["samples"])
        active=labels if task["design"]=="clade_grouped" else None
        for number,repeat,fold,train,test in repeated_outer_splits(y,task["n_splits"],5,42,task["design"],active):
            for partition, indices in [("Training",train),("Test",test)]:
                for i in indices:
                    outer.append(dict(task=task["task"],scope=task["scope"],design=task["design"],
                                      input_type=task["input_type"],fold_number=number,repeat=repeat,fold=fold,
                                      sample=task["samples"][i],partition=partition,cluster=labels[i],
                                      label="R" if y[i] else "S"))
            for inn,(tr,te) in enumerate(inner_splits(y[train],None if active is None else active[train],42+number),1):
                for partition, indices in [("Training",train[tr]),("Validation",train[te])]:
                    for i in indices:
                        inner.append(dict(task=task["task"],fold_number=number,inner_fold=inn,
                                          sample=task["samples"][i],partition=partition,seed=42+number))
    source=args.outdir/"source_data"
    source.mkdir(exist_ok=True)
    pd.DataFrame(outer).to_csv(source/"outer_sample_assignments.tsv.gz",sep="\t",index=False)
    pd.DataFrame(inner).to_csv(source/"inner_sample_assignments.tsv.gz",sep="\t",index=False)
    pd.DataFrame([{k:v for k,v in t.items() if k!="samples"}|{"n_samples":len(t["samples"]),
                   "n_outer_folds":t["n_splits"]*t["n_repeats"]} for t in tasks]).to_csv(source/"tasks.tsv",sep="\t",index=False)
    print(f"Prepared {len(tasks)} tasks, {sum(t['n_splits']*5 for t in tasks)} outer models",flush=True)


def destination(root,task):
    return root/("results" if task["scope"]=="prediction" else "within_group")/task["task"]


def work(args):
    config=json.loads((args.outdir/"run_configuration.json").read_text())
    metadata=pd.read_csv(config["metadata"],sep="\t",dtype={"sample":str})
    loaded={}
    current=None
    while True:
        claimed=None
        for task in config["tasks"]:
            if args.only_input and task["input_type"] != args.only_input: continue
            if args.only_scope and task["scope"] != args.only_scope: continue
            target=destination(args.outdir,task)
            (target/"folds").mkdir(parents=True,exist_ok=True)
            for number in range(1,task["n_splits"]*5+1):
                fold=target/"folds"/f"{number:02d}"
                if (fold/"run_info.json").exists(): continue
                lock=target/"folds"/f".{number:02d}.lock"
                try: lock.mkdir()
                except FileExistsError: continue
                if (fold/"run_info.json").exists():
                    lock.rmdir()
                    continue
                (lock/"owner.json").write_text(json.dumps({"pid":os.getpid(),"started":time.time()}))
                claimed=(task,number,fold,lock)
                break
            if claimed: break
        if claimed is None: break
        task,number,fold,lock=claimed
        kind=task["input_type"]
        try:
            if current!=task["task"]:
                if kind not in loaded:
                    loaded[kind]=load_cache(config["input_caches"][kind])
                matrix,features,all_samples=loaded[kind]
                lookup={sample:i for i,sample in enumerate(all_samples)}
                samples=task["samples"]
                rows=np.array([lookup[s] for s in samples])
                x=np.asarray(matrix[rows],dtype=np.int8)
                y,groups=aligned_metadata(metadata,samples)
                current=task["task"]
            if fold.exists():
                raise RuntimeError("Incomplete fold requires inspection: "+str(fold))
            print(datetime.now().isoformat(timespec="seconds"),task["task"],number,"start",flush=True)
            DRIVER.run_analysis(x,features,samples,y,groups,kind,task["design"],fold,
                                n_splits=task["n_splits"],n_repeats=5,seed=42,inner_folds=3,
                                n_estimators=500,n_jobs=args.jobs,top_k=10000,max_missing=.10,
                                selected_fold_numbers=[number])
            print(datetime.now().isoformat(timespec="seconds"),task["task"],number,"complete",flush=True)
        except BaseException:
            (lock/"error.txt").write_text(traceback.format_exc())
            raise
        else:
            (lock/"owner.json").unlink()
            lock.rmdir()
    print("Worker has no unclaimed outer folds",flush=True)


def aggregate_task(root,task):
    target=destination(root,task)
    folders=[target/"folds"/f"{number:02d}" for number in range(1,task["n_splits"]*5+1)]
    if not all((p/"run_info.json").exists() for p in folders): return False
    context={k:task[k] for k in ["group_key","group","cluster"] if k in task}
    if task["scope"]=="within_group": context["method_key"]="tuned"
    for name in sorted(p.name for p in folders[0].glob("*.tsv")):
        table=pd.concat([pd.read_csv(p/name,sep="\t",low_memory=False) for p in folders],ignore_index=True)
        for key,value in context.items(): table[key]=value
        if "gene_symbol" in table and "gene_id" in table:
            table.loc[table.gene_id.astype(str).str.removesuffix("_A").eq("C1_02420C"), "gene_symbol"] = "FKS1"
        table.to_csv(target/name,sep="\t",index=False)
    by_fold=pd.read_csv(target/"feature_importance_by_fold.tsv",sep="\t",low_memory=False)
    key="feature_idx" if task["input_type"]=="pc" else "variant_idx"
    other=[c for c in by_fold.columns if c not in {"gini_importance","repeat","fold","fold_number"}]
    metadata=by_fold[other].drop_duplicates(key)
    summary=by_fold.groupby(key,as_index=False).agg(importance_sum=("gini_importance","sum"),
                                                   n_folds_present=("fold_number","nunique"))
    summary["mean_importance_across_folds"]=summary.importance_sum/len(folders)
    summary["rank"]=summary.mean_importance_across_folds.rank(method="first",ascending=False).astype(int)
    summary=summary.merge(metadata,on=key,validate="one_to_one")
    summary.to_csv(target/"feature_importance_summary.tsv",sep="\t",index=False)
    info=json.loads((folders[0]/"run_info.json").read_text())
    info["completed_fold_numbers"]=list(range(1,len(folders)+1))
    info.update(context)
    if task["scope"] == "within_group":
        info.update(samples=info["n_samples"], R=info["n_resistant"], S=info["n_susceptible"],
                    pca_cluster=task["cluster"])
    info["input_cache"]=json.loads((root/"run_configuration.json").read_text())["input_caches"][task["input_type"]]
    (target/"run_info.json").write_text(json.dumps(info,indent=2)+"\n")
    return True


def summarize(args):
    config=json.loads((args.outdir/"run_configuration.json").read_text())
    finished=[]
    for task in config["tasks"]:
        if aggregate_task(args.outdir,task): finished.append(task)
    if len(finished)!=len(config["tasks"]):
        raise ValueError(f"Only {len(finished)}/{len(config['tasks'])} task results are complete")
    for scope in ["prediction","within_group"]:
        tasks=[t for t in finished if t["scope"]==scope]
        folder=args.outdir/("results" if scope=="prediction" else scope)
        for name in ["cv_fold_metrics.tsv","cv_fold_predictions.tsv","feature_importance_summary.tsv"]:
            frames=[]
            for task in tasks:
                frame=pd.read_csv(destination(args.outdir,task)/name,sep="\t",low_memory=False)
                frame["task"]=task["task"]
                frame["input_key"]=task["input_type"]
                frame["input_type"]=LABELS[task["input_type"]]
                frames.append(frame)
            result=pd.concat(frames,ignore_index=True)
            result.to_csv(folder/("combined_"+name),sep="\t",index=False)
            if name=="feature_importance_summary.tsv":
                result.loc[result['rank']<=20].to_csv(folder/"top20_features.tsv",sep="\t",index=False)
    config["completed_utc"]=datetime.now(timezone.utc).isoformat()
    config["status"]="complete"
    (args.outdir/"run_configuration.json").write_text(json.dumps(config,indent=2)+"\n")
    print("All repeated-validation tasks are complete",flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir",type=Path,required=True)
    parser.add_argument("--prepare",action="store_true")
    parser.add_argument("--summarize",action="store_true")
    parser.add_argument("--inputs",type=Path)
    parser.add_argument("--metadata",type=Path)
    parser.add_argument("--cohort-metadata",type=Path)
    parser.add_argument("--groups",type=Path)
    parser.add_argument("--jobs",type=int,default=2)
    parser.add_argument("--only-input",choices=["variant","gene","pc"])
    parser.add_argument("--only-scope",choices=["prediction","within_group"])
    args=parser.parse_args()
    if args.prepare: prepare(args)
    elif args.summarize: summarize(args)
    else: work(args)


if __name__=="__main__":main()
