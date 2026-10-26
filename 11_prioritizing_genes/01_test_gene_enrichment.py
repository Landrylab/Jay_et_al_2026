"""Test prepared network-gene burdens in drug-specific eligible cohorts."""
from pathlib import Path
import hashlib
import json
import subprocess
import argparse
import shutil

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, false_discovery_control
from scipy.stats.contingency import odds_ratio

HERE = Path(__file__).resolve().parent
OUT = None
DRUGS = ['fluconazole', 'voriconazole', 'posaconazole', 'caspofungin', 'anidulafungin', 'micafungin']
AZOLES = DRUGS[:3]


def bool_series(s):
    return s.astype(str).str.lower().isin(['true', '1'])


def collapse(a):
    if a.shape[0] == 0:
        raise ValueError('No genes in requested pool')
    return np.where((a == 1).any(axis=0), 1, np.where((a < 0).any(axis=0), -1, 0)).astype(np.int8)


def adjust(s):
    q = np.array(false_discovery_control(s.fillna(1).to_numpy(), method='bh'), copy=True)
    q[s.isna().to_numpy()] = np.nan
    return q


def calculate(frame):
    r, s, c = frame.phenotype.eq('R'), frame.phenotype.eq('S'), frame.carrier
    a, b = int((r & c.eq(1)).sum()), int((r & c.eq(0)).sum())
    cc, d = int((s & c.eq(1)).sum()), int((s & c.eq(0)).sum())
    tab = np.array([[a, b], [cc, d]])
    valid = bool(tab.sum(axis=0).min() > 0 and tab.sum(axis=1).min() > 0)
    fit = fisher_exact(tab, alternative='two-sided') if valid else None
    cond = odds_ratio(tab) if valid else None
    ci = cond.confidence_interval() if cond else None
    pct = lambda x, n: 100 * x / n if n else np.nan
    return dict(R_carrier=a, R_noncarrier=b, R_unknown=int((r & c.lt(0)).sum()), R_total=int(r.sum()),
                S_carrier=cc, S_noncarrier=d, S_unknown=int((s & c.lt(0)).sum()), S_total=int(s.sum()),
                R_carrier_pct_all=pct(a, int(r.sum())), S_carrier_pct_all=pct(cc, int(s.sum())),
                R_carrier_pct_called=pct(a, a+b), S_carrier_pct_called=pct(cc, cc+d),
                odds_ratio=float(fit.statistic) if fit else np.nan,
                conditional_odds_ratio=float(cond.statistic) if cond else np.nan,
                or_ci95_low=float(ci.low) if ci else np.nan, or_ci95_high=float(ci.high) if ci else np.nan,
                fisher_p=float(fit.pvalue) if fit else np.nan, testable=valid)


def export_postselection_coverage(results, calls, negative):
    """Descriptive coverage after same-cohort selection; deliberately no tests."""
    rows, sample_rows = [], []
    for selection, flag in [('pca_cluster_significant', 'pca_cluster_resistance_enriched'),
                            ('pca_cluster_study_significant', 'pca_cluster_study_resistance_enriched')]:
        for drug in DRUGS:
            chosen = results[results.drug.eq(drug) & results[flag]].copy()
            frame = negative[negative.drug.eq(drug)][['sample','drug','phenotype','cluster_label','Study_clean']].copy()
            ids = chosen.gene_id.tolist()
            if ids:
                selected = calls[calls.drug.eq(drug) & calls.gene_id.isin(ids)]
                matrix = selected.pivot(index='gene_id',columns='sample',values='carrier').reindex(index=ids,columns=frame['sample'])
                assert matrix.notna().all().all()
                frame['carrier'] = collapse(matrix.to_numpy())
            else:
                frame['carrier'] = 0
            frame['selection'] = selection
            sample_rows.append(frame)
            record = dict(selection=selection,drug=drug,n_genes=len(ids),genes=';'.join(chosen.gene),gene_ids=';'.join(ids),
                          interpretation='Descriptive carrier coverage after selection in the same cohort; not independent validation or causal explanation.')
            for ph in ['R','S']:
                a = frame.loc[frame.phenotype.eq(ph),'carrier']
                for label, value in [('total',len(a)),('carrier',int(a.eq(1).sum())),('noncarrier',int(a.eq(0).sum())),('unknown',int(a.lt(0).sum()))]:
                    record[f'{ph}_{label}'] = value
                record[f'{ph}_carrier_pct_all'] = 100*a.eq(1).sum()/len(a) if len(a) else np.nan
            rows.append(record)
    pd.DataFrame(rows).to_csv(OUT/'postselection_carrier_coverage.csv',index=False)
    pd.concat(sample_rows,ignore_index=True).to_csv(OUT/'postselection_sample_calls.csv.gz',index=False)



def main():
    global OUT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, required=True)
    parser.add_argument('--outdir', type=Path, required=True)
    parser.add_argument('--rscript', default='Rscript')
    args = parser.parse_args()
    OUT = args.outdir.resolve()
    source = args.input_dir.resolve()
    if OUT == source:
        parser.error('Input and output directories must differ')
    if not shutil.which(args.rscript):
        parser.error('Rscript is required; install R or pass --rscript')
    OUT.mkdir(parents=True, exist_ok=True)
    hc_path = source / 'high_confidence_classifications.csv'
    inv_path = source / 'network_gene_inventory.csv'
    hc = pd.read_csv(hc_path)
    inventory = pd.read_csv(inv_path).set_index('gene_id', drop=False)
    for flag in ['in_gwas_network', 'in_rf_network']:
        inventory[flag] = bool_series(inventory[flag])
    assert not hc.duplicated(['sample', 'antifungal']).any()
    assert len(hc) > 0 and set(hc.phenotype).issubset({'R', 'S'})
    hc['drug'] = hc.antifungal.str.lower()
    assert set(hc.high_confidence_status).issubset({'present', 'absent', 'unknown'})
    negative = hc[hc.high_confidence_status.eq('absent')].copy()
    assert negative.cluster_label.notna().all() and negative.Study_clean.notna().all()
    assert negative.cluster_label.astype(str).str.match(r'^C\d+$').all()
    negative.to_csv(OUT/'analysis_cohort.csv', index=False)
    cohort_counts = []
    for drug in DRUGS:
        for phenotype in ['R', 'S']:
            z = hc[hc.drug.eq(drug) & hc.phenotype.eq(phenotype)]
            cohort_counts.append(dict(drug=drug, phenotype=phenotype, total=len(z),
                                     hc_present=int(z.high_confidence_status.eq('present').sum()),
                                     hc_absent=int(z.high_confidence_status.eq('absent').sum()),
                                     hc_unknown=int(z.high_confidence_status.eq('unknown').sum())))
    pd.DataFrame(cohort_counts).to_csv(OUT/'cohort_counts.csv', index=False)

    gene_calls, pool_calls, gene_results, pool_results, input_paths = [], [], [], [], [hc_path, inv_path]
    core = ['sample', 'drug', 'phenotype', 'cluster_label', 'Study_clean']
    for drug in DRUGS:
        carrier_path = source/'carriers'/f'{drug}.npz'
        input_paths.append(carrier_path)
        with np.load(carrier_path, allow_pickle=False) as data:
            codes = data['codes']; genes = data['gene_ids'].astype(str).tolist(); samples = data['sample_ids'].astype(str).tolist()
        assert codes.shape == (len(genes),len(samples)) and np.isin(codes, [-1,0,1]).all()
        assert len(genes) == len(set(genes)) and len(samples) == len(set(samples))
        gi, si = {g:i for i,g in enumerate(genes)}, {s:i for i,s in enumerate(samples)}
        frame = negative[negative.drug.eq(drug)].copy().reset_index(drop=True)
        assert frame['sample'].isin(samples).all()
        cols = [si[s] for s in frame['sample']]
        active = inventory[inventory.in_gwas_network | (inventory.in_rf_network & (drug in AZOLES))]
        for spec in active.to_dict('records'):
            gene = spec['gene_id']
            row = frame[core].copy()
            row['gene_id'], row['gene'] = gene, spec['gene']
            row['carrier'] = codes[gi[gene],cols]
            gene_calls.append(row)
            original_drugs = str(spec.get('original_gwas_drugs', '')).split(';')
            edge = (drug in original_drugs) or (bool(spec['in_rf_network']) and drug == 'fluconazole')
            record = {**spec, 'drug':drug, 'original_network_edge':edge, **calculate(row)}
            gene_results.append(record)
        pools = {'gwas_all':inventory.index[inventory.in_gwas_network].tolist()}
        if drug in AZOLES:
            pools['rf_all'] = inventory.index[inventory.in_rf_network].tolist()
            pools['gwas_rf_union'] = active.index.tolist()
        for name, ids in pools.items():
            row = frame[core].copy()
            row['pool'] = name
            row['carrier'] = collapse(codes[[gi[g] for g in ids]][:,cols])
            pool_calls.append(row)
            pool_results.append(dict(drug=drug,pool=name,n_genes=len(ids),gene_ids=';'.join(ids),**calculate(row)))

    calls = pd.concat(gene_calls, ignore_index=True)
    pools = pd.concat(pool_calls, ignore_index=True)
    assert not calls.duplicated(['drug','sample','gene_id']).any()
    assert not pools.duplicated(['drug','sample','pool']).any()
    calls.to_csv(OUT/'sample_gene_calls.csv.gz',index=False)
    pools.to_csv(OUT/'sample_pool_calls.csv.gz',index=False)
    results = pd.DataFrame(gene_results)
    assert not results.duplicated(['drug','gene_id']).any()
    results['fisher_q_primary'] = adjust(results.fisher_p)
    results['fisher_q_original_edges'] = np.nan
    edge = results.original_network_edge
    results.loc[edge,'fisher_q_original_edges'] = adjust(results.loc[edge,'fisher_p'])
    pool_result = pd.DataFrame(pool_results)
    pool_result['fisher_q_primary'] = adjust(pool_result.fisher_p)
    results.to_csv(OUT/'gene_results.csv',index=False)
    pool_result.to_csv(OUT/'pool_results.csv',index=False)
    subprocess.run([args.rscript, '--vanilla', str(HERE/'helpers/gene_tests.R'), str(OUT)], check=True)

    results = pd.read_csv(OUT/'gene_results.csv')
    pool_result = pd.read_csv(OUT/'pool_results.csv')
    for method, or_col in [('fisher','odds_ratio'),('pca_cluster','pca_cluster_or'),('pca_cluster_study','pca_cluster_study_or')]:
        results[f'{method}_resistance_enriched'] = results[f'{method}_q_primary'].lt(.05) & results[or_col].gt(1)
        results[f'{method}_susceptibility_enriched'] = results[f'{method}_q_primary'].lt(.05) & results[or_col].lt(1)
    results.to_csv(OUT/'gene_results.csv',index=False)
    positive_direction = results.odds_ratio.gt(1) | results.pca_cluster_or.gt(1) | results.pca_cluster_study_or.gt(1)
    positive = results[positive_direction].sort_values(['pca_cluster_study_q_primary','fisher_q_primary'],na_position='last')
    positive.to_csv(OUT/'positive_enrichment_results.csv',index=False)
    significant = positive[positive[['fisher_resistance_enriched','pca_cluster_resistance_enriched','pca_cluster_study_resistance_enriched']].any(axis=1)]
    significant.to_csv(OUT/'positive_enrichment_signals.csv',index=False)
    results[results.pca_cluster_resistance_enriched].sort_values('pca_cluster_q_primary').to_csv(OUT/'primary_pca_enriched_genes.csv',index=False)
    results[results.pca_cluster_susceptibility_enriched].sort_values('pca_cluster_q_primary').to_csv(OUT/'primary_pca_susceptibility_enriched_genes.csv',index=False)
    export_postselection_coverage(results,calls,negative)
    results[results.original_network_edge.astype(str).str.lower().eq('true')].to_csv(OUT/'original_edge_sensitivity_results.csv',index=False)
    audit = dict(question='Are prepared gene burdens enriched among R phenotype versus S phenotype isolates in the eligible drug-specific cohorts?',
                 source_cohort_rows=len(hc),source_cohort_isolates=hc['sample'].nunique(),negative_rows=len(negative),negative_unique_samples=negative['sample'].nunique(),
                 gwas_genes=int(inventory.in_gwas_network.sum()),rf_genes=int(inventory.in_rf_network.sum()),union_genes=len(inventory),
                 primary_gene_drug_tests=len(results),original_edge_gene_drug_tests=int(edge.sum()),primary_pool_tests=len(pool_result),
                 high_confidence='Prepared classifications follow the Methods-defined drug-specific evidence rules.',
                 gene_carrier='Prepared fixed per-drug HIGH/MODERATE burden with original-cohort allele carrier frequency >0 and <=5%; Methods-defined variant selection. Present overrides missing; otherwise any missing qualifying allele gives unknown.',
                 test='Two-sided Fisher exact on known gene calls; sample odds ratio plus conditional exact95% CI. PCA-cluster and PCA-cluster×study exact conditional sensitivities.',
                 correction='BH over all primary gene×drug tests separately for each statistical method. Original-edge subset separate BH sensitivity. Twelve pooled comparisons form a separate family. Untestable entries count as p=1 internally and retain missing reported p/q.',
                 pca_cluster_source='cluster_label supplied in the prepared classification table; PCA cluster assignments.',
                 validation='Every Fisher p value and primary BH q value independently reproduced in R during execution. Unique gene-drug tests and sample-level joins asserted.',
                 postselection_coverage='Descriptive union of carriers among genes selected by each adjustment method; no additional p values or claims of independent validation.',
                 limitations=['Exploratory reanalysis of genes selected using the same cohort; no independent validation or causal attribution.','RF importance is unsigned. Association direction is estimated here.','PCA clusters and study stratification do not fully account for genetic relatedness.','Gene and pooled burdens can contain resistance-associated, neutral, and susceptibility-associated variants.','Unknown gene calls excluded from each corresponding test, shown in counts.'],
                 phenotype_display_mapping={'R':'R phenotype','S':'S phenotype'},
                 inputs={str(p.relative_to(source)):hashlib.sha256(p.read_bytes()).hexdigest() for p in input_paths})
    (OUT/'enrichment_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(json.dumps({k:audit[k] for k in ['negative_rows','negative_unique_samples','gwas_genes','rf_genes','union_genes','primary_gene_drug_tests','original_edge_gene_drug_tests','primary_pool_tests']},indent=2))
    print(significant[['drug','gene','R_carrier','R_total','R_unknown','S_carrier','S_total','S_unknown','odds_ratio','fisher_q_primary','pca_cluster_q_primary','pca_cluster_study_q_primary']].to_string(index=False))


if __name__ == '__main__':
    main()
