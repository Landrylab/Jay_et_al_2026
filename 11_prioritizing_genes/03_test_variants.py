"""Test prepared individual-variant calls and report conditional gene follow-up."""
from pathlib import Path
import argparse
import hashlib
import json
import shutil
import subprocess

import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, required=True)
    parser.add_argument('--gene-results-dir', type=Path, required=True)
    parser.add_argument('--conditional-results', type=Path, required=True)
    parser.add_argument('--outdir', type=Path, required=True)
    parser.add_argument('--rscript', default='Rscript')
    args = parser.parse_args()
    if not shutil.which(args.rscript):
        parser.error('Rscript is required; install R or pass --rscript')
    out = args.outdir.resolve()
    if out == args.input_dir.resolve():
        parser.error('Input and output directories must differ')
    out.mkdir(parents=True, exist_ok=True)
    paths = [args.input_dir/'variant_metadata.csv', args.input_dir/'sample_variant_calls.csv.gz']
    metadata = pd.read_csv(paths[0])
    calls = pd.read_csv(paths[1])
    screen = pd.read_csv(args.gene_results_dir/'gene_results.csv')
    genes = pd.read_csv(args.gene_results_dir/'sample_gene_calls.csv.gz')
    cohort = pd.read_csv(args.gene_results_dir/'analysis_cohort.csv')
    conditional = pd.read_csv(args.conditional_results, float_precision='round_trip')
    qualified = screen[screen.pca_cluster_q_primary.lt(.05) & screen.pca_cluster_or.gt(1)]
    qualified_keys = set(zip(qualified.drug, qualified.gene_id))
    assert set(zip(metadata.drug, metadata.gene_id)) == qualified_keys, 'Variant families must cover the original qualifying gene screen'
    assert not metadata.key.duplicated().any()
    assert not calls.duplicated(['key', 'sample']).any()
    assert set(calls.key) == set(metadata.key)
    assert calls.carrier.isin([-1, 0, 1]).all()
    assert set(calls.phenotype).issubset({'R', 'S'})
    assert metadata.family.eq(metadata.drug + '|' + metadata.gene_id).all()

    # Verify that the prepared individual calls reconstruct each fixed gene burden.
    for (drug, gene_id), family in metadata.groupby(['drug', 'gene_id'], sort=False):
        c = cohort[cohort.drug.eq(drug)].set_index('sample')
        v = calls[calls.key.isin(family.key)]
        matrix = v.pivot(index='key', columns='sample', values='carrier').reindex(index=family.key, columns=c.index)
        assert matrix.notna().all().all(), (drug, gene_id, 'Missing variant calls')
        assert set(v['sample']) == set(c.index), (drug, gene_id, 'Unexpected samples')
        for column in ['phenotype', 'cluster_label']:
            mapping = c[column].to_dict()
            assert v[column].eq(v['sample'].map(mapping)).all(), (drug, gene_id, column)
        codes = matrix.to_numpy()
        union = np.where((codes == 1).any(axis=0), 1, np.where((codes < 0).any(axis=0), -1, 0))
        g = genes[genes.drug.eq(drug) & genes.gene_id.eq(gene_id)].set_index('sample')
        assert np.array_equal(union, g.loc[c.index, 'carrier']), (drug, gene_id, 'Burden reconstruction')

    metadata.to_csv(out/'variant_metadata.csv', index=False)
    calls.to_csv(out/'sample_variant_calls.csv.gz', index=False)
    subprocess.run([args.rscript, '--vanilla', str(Path(__file__).parent/'helpers/variant_tests.R'), str(out)], check=True)
    results = pd.read_csv(out/'variant_results.csv')
    results['R_minus_S_percentage_points'] = 100 * (results.R_carrier / results.R_called - results.S_carrier / results.S_called)
    results['enriched'] = results.pca_q_within_gene_drug.lt(.05) & results.pca_or.gt(1)
    results.to_csv(out/'variant_results.csv', index=False)

    retained = conditional[
        conditional.threshold.eq(.20) & conditional.method.eq('pca_joint') &
        conditional.n_partners.gt(0) & conditional.q_36_tests.lt(.05) & conditional['or'].gt(1)
    ].copy()
    assert not retained.duplicated(['drug', 'gene_id']).any()
    retained.to_csv(out/'retained_gene_associations.csv', index=False)
    retained_keys = set(zip(retained.drug, retained.gene_id))
    selected = results[[key in retained_keys for key in zip(results.drug, results.gene_id)]].copy()
    selected.to_csv(out/'selected_gene_variant_results.csv', index=False)
    prioritized = selected[selected.enriched].copy()
    prioritized.to_csv(out/'prioritized_variants.csv', index=False)

    amino = dict(zip(['Ala','Arg','Asn','Asp','Cys','Gln','Glu','Gly','His','Ile','Leu','Lys','Met','Phe','Pro','Ser','Thr','Trp','Tyr','Val'], list('ARNDCQEGHILKMFPSTWYV')))
    def label(row):
        value = row.protein_change
        if pd.isna(value) or not value:
            return row.variant_id
        value = str(value).removeprefix('p.')
        for source, target in amino.items():
            value = value.replace(source, target)
        return value
    summaries = []
    for row in retained.itertuples():
        family = selected[selected.drug.eq(row.drug) & selected.gene_id.eq(row.gene_id)]
        enriched = family[family.enriched].sort_values(['R_minus_S_percentage_points','pca_q_within_gene_drug','variant_id'], ascending=[False,True,True])
        names = [label(v) for v in enriched.itertuples()]
        names = [n + ' (' + str(v.hgvs_c) + ')' if names.count(n) > 1 else n for n, v in zip(names, enriched.itertuples())]
        summaries.append(dict(drug=row.drug, gene_id=row.gene_id, gene=row.gene,
                              n_variants_tested=len(family), n_enriched=len(enriched),
                              enriched_variants=', '.join(names) if names else 'none',
                              variant_ids=';'.join(enriched.variant_id)))
    fields = ['drug','gene_id','gene','n_variants_tested','n_enriched','enriched_variants','variant_ids']
    pd.DataFrame(summaries, columns=fields).to_csv(out/'retained_gene_summary.csv', index=False)
    audit = dict(method='Two-sided exact conditional test stratified by original PCA clusters; BH within each gene–drug family. Untestable variants enter correction as p=1 and retain missing reported p/q.',
                 variant_selection=f'The original qualifying variants and {len(qualified)} screened gene–drug families are retained; gene selection does not change within-family multiplicity.',
                 displayed_gene_selection='Primary joint co-carriage threshold 0.20, at least one selected partner, positive conditional odds ratio and follow-up q < 0.05.',
                 scope='Prepared fixed rare HIGH/MODERATE variant calls in the original eligible drug-specific cohorts.',
                 limitation='Exploratory same-cohort follow-up. Individual-variant tests adjust for PCA clusters, not co-carried genes. Gene selection uses the co-carriage adjustment.',
                 gene_drug_rows=len(qualified), variant_tests=len(results), enriched_variants_all_screened_genes=int(results.enriched.sum()),
                 retained_gene_associations=len(retained), prioritized_variants=len(prioritized),
                 gene_call_reconstruction='All variant families exactly reproduce the supplied gene-carriage codes.',
                 phenotype_display_mapping={'R':'R phenotype','S':'S phenotype'},
                 input_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    (out/'analysis_metadata.json').write_text(json.dumps(audit, indent=2)+'\n')
    print(json.dumps({k:audit[k] for k in ['gene_drug_rows','variant_tests','retained_gene_associations','prioritized_variants']}, indent=2))


if __name__ == '__main__':
    main()
