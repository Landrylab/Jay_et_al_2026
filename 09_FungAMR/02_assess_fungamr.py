#!/usr/bin/env python3
"""Assess drug-matched FungAMR carriers in a study-selected phenotype cohort.

This stage reads the full carrier catalogue. Study exclusions
apply only to the phenotype cohort used for Figure 2b and its association tests.
"""
import argparse
import hashlib
import json
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy.stats import fisher_exact
from fungamr_evidence import SIGNED_EVIDENCE_RULE


TRAITS = {
    'Fluconazole': 'rbin_flu', 'Voriconazole': 'rbin_vor',
    'Posaconazole': 'rbin_pos', 'Caspofungin': 'rbin_cas',
    'Anidulafungin': 'rbin_ani', 'Micafungin': 'rbin_mic',
}
GENE_COLUMNS = [
    'antifungal', 'gene', 'count', 'denominator', 'percentage',
    'resistant_samples_with_gene_mutation', 'resistant_n',
    'resistant_samples_with_gene_mutation_pct',
]


def sha(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def percentage(numerator, denominator):
    if not denominator:
        return 'n/a'
    value = Decimal(100) * Decimal(int(numerator)) / Decimal(int(denominator))
    return str(value.quantize(Decimal('.1'), rounding=ROUND_HALF_UP))


def require_unique_samples(table, name):
    if 'sample' not in table or table['sample'].isna().any():
        raise ValueError(f'{name} requires nonmissing sample identifiers')
    if table['sample'].duplicated().any():
        raise ValueError(f'{name} contains duplicate sample identifiers')


def benjamini_hochberg(pvalues):
    """Adjust only evaluable tests; preserve the input ordering."""
    pvalues = np.asarray(pvalues, dtype=float)
    if not len(pvalues):
        return pvalues
    order = np.argsort(pvalues)
    scaled = pvalues[order] * len(pvalues) / np.arange(1, len(pvalues) + 1)
    adjusted = np.empty(len(pvalues))
    adjusted[order] = np.minimum(1, np.minimum.accumulate(scaled[::-1])[::-1])
    return adjusted


def build_assessment(matches_path, phenotypes_path, metadata_path,
                     exclude_studies=(), sample_ids=None):
    matches = pd.read_csv(matches_path, dtype=str)
    pheno = pd.read_csv(phenotypes_path, sep='\t', dtype=str)
    metadata = pd.read_csv(metadata_path, sep='\t', dtype=str)
    require_unique_samples(pheno, 'Phenotypes')
    require_unique_samples(metadata, 'Metadata')
    for name in TRAITS.values():
        if name not in pheno:
            raise ValueError(f'Missing phenotype column: {name}')
    if 'Study_clean' not in metadata:
        raise ValueError('Metadata requires a Study_clean column')
    for name in ['sample', 'gene', 'antifungal', 'confidence_score', 'mutation_found']:
        if name not in matches:
            raise ValueError(f'Missing carrier column: {name}')
    if matches[['sample', 'gene', 'antifungal']].isna().any().any():
        raise ValueError('Carrier sample, gene and antifungal values must be nonmissing')

    # The matching stage can constrain assessment to the master genotype samples.
    # A standalone assessment uses the supplied phenotype table as its cohort.
    pheno = pheno.set_index('sample')
    if sample_ids is not None:
        sample_ids = list(sample_ids)
        if len(sample_ids) != len(set(sample_ids)):
            raise ValueError('The master sample list contains duplicate identifiers')
        pheno = pheno.reindex(sample_ids)
        pheno.index.name = 'sample'
    cohort = pd.DataFrame(index=pheno.index)
    cohort['Study_clean'] = metadata.set_index('sample')['Study_clean'].reindex(pheno.index)
    if cohort['Study_clean'].isna().any():
        missing = cohort.index[cohort['Study_clean'].isna()].tolist()
        raise ValueError(f'No study assignment for {len(missing)} cohort samples: {missing[:5]}')
    if 'Study_clean' in pheno:
        conflict = pheno['Study_clean'].notna() & pheno['Study_clean'].ne(cohort['Study_clean'])
        if conflict.any():
            raise ValueError('Phenotype and metadata study assignments disagree')
    exclusions = sorted(set(exclude_studies))
    unknown = set(exclusions) - set(metadata['Study_clean'].dropna())
    if unknown:
        raise ValueError(f'Excluded study labels were not found in metadata: {sorted(unknown)}')
    binary = pheno[list(TRAITS.values())].replace({'R': '1', 'S': '0'})
    cohort['excluded_by_study'] = cohort['Study_clean'].isin(exclusions)
    cohort['included_in_assessment'] = ~cohort['excluded_by_study']
    cohort['n_nonmissing_drug_phenotypes'] = binary.isin(['0', '1']).sum(axis=1)
    cohort['has_any_nonmissing_phenotype'] = cohort['n_nonmissing_drug_phenotypes'].gt(0)
    cohort['included_in_panel_B'] = (
        cohort['included_in_assessment'] & cohort['has_any_nonmissing_phenotype']
    )

    # Input catalogue must already resolve signed evidence before +1/+2 filtering.
    # Carrier sets
    # collapse repeated database evidence, alleles and genes to one isolate/drug.
    raw_match_rows = len(matches)
    matches['gene'] = matches['gene'].replace({'GSC1': 'FKS1'})
    scores = pd.to_numeric(matches['confidence_score'], errors='coerce')
    matches = matches.loc[
        scores.isin([1, 2]) & matches['mutation_found'].ne('V488I')
    ].copy()
    matches['antifungal_norm'] = matches['antifungal'].str.casefold()
    cells, gene_rows, metric_rows, isolate_rows, fisher_rows = [], [], [], [], []
    for drug, trait in TRAITS.items():
        phenotype = binary.loc[cohort['included_in_assessment'], trait]
        valid = phenotype.isin(['0', '1'])
        drug_matches = matches.loc[matches['antifungal_norm'].eq(drug.casefold())]
        positive = set(drug_matches['sample'])
        resistant = set(phenotype.index[phenotype.eq('1')])
        susceptible = set(phenotype.index[phenotype.eq('0')])
        tp, fp, fn, tn = map(len, [resistant & positive, susceptible & positive,
                                  resistant - positive, susceptible - positive])
        present_n, absent_n = tp + fp, fn + tn
        total = int(valid.sum())
        metric = dict(
            antifungal=drug, trait=trait, n_nonmissing=total,
            resistant_n=len(resistant), susceptible_n=len(susceptible),
            mutation_present_n=present_n, mutation_absent_n=absent_n,
            tp_resistant_with_mutation=tp, fn_resistant_without_mutation=fn,
            fp_susceptible_with_mutation=fp, tn_susceptible_without_mutation=tn,
        )
        for name, n, d in [
            ('resistance_explained_pct', tp, len(resistant)),
            ('resistance_unexplained_pct', fn, len(resistant)),
            ('susceptible_with_mutation_pct', fp, len(susceptible)),
            ('susceptible_without_mutation_pct', tn, len(susceptible)),
            ('ppv_pct', tp, present_n), ('npv_pct', tn, absent_n),
            ('tp_total_pct', tp, total), ('fn_total_pct', fn, total),
            ('fp_total_pct', fp, total), ('tn_total_pct', tn, total),
            ('mutation_present_resistant_pct', tp, present_n),
            ('mutation_present_susceptible_pct', fp, present_n),
            ('mutation_absent_resistant_pct', fn, absent_n),
            ('mutation_absent_susceptible_pct', tn, absent_n),
        ]:
            metric[name] = 100 * n / d if d else np.nan
        metric_rows.append(metric)
        for state, phen, count, denominator in [
            ('Present', 'Resistant', tp, present_n),
            ('Present', 'Susceptible', fp, present_n),
            ('Absent', 'Resistant', fn, absent_n),
            ('Absent', 'Susceptible', tn, absent_n),
        ]:
            cells.append(dict(antifungal=drug, mutation=state, phenotype=phen,
                              count=count, denominator=denominator,
                              percentage=percentage(count, denominator)))
        for gene, group in drug_matches.groupby('gene'):
            count = len(resistant & set(group['sample']))
            if count:
                gene_rows.append(dict(
                    antifungal=drug, gene=gene, count=count, denominator=len(resistant),
                    percentage=percentage(count, len(resistant)),
                    resistant_samples_with_gene_mutation=count, resistant_n=len(resistant),
                    resistant_samples_with_gene_mutation_pct=100 * count / len(resistant),
                ))
        for sample in phenotype.index[valid]:
            isolate_rows.append(dict(sample=sample, antifungal=drug,
                                     phenotype='R' if sample in resistant else 'S',
                                     mutation_present=sample in positive))

        # Rows are mutation present/absent; columns are resistant/susceptible.
        # A zero row or column margin supplies no association information.
        margins = [present_n, absent_n, len(resistant), len(susceptible)]
        evaluable = min(margins) > 0
        odds_ratio, pvalue = (fisher_exact([[tp, fp], [fn, tn]], alternative='two-sided')
                              if evaluable else (np.nan, np.nan))
        fisher_rows.append(dict(
            antifungal=drug, trait=trait, n_nonmissing=total,
            excluded_nonmissing_n=int((binary[trait].isin(['0', '1']) & cohort['excluded_by_study']).sum()),
            mutation_present_resistant_n=tp, mutation_present_susceptible_n=fp,
            mutation_absent_resistant_n=fn, mutation_absent_susceptible_n=tn,
            mutation_present_n=present_n, mutation_absent_n=absent_n,
            resistant_n=len(resistant), susceptible_n=len(susceptible),
            mutation_present_resistant_pct=metric['mutation_present_resistant_pct'],
            mutation_absent_resistant_pct=metric['mutation_absent_resistant_pct'],
            odds_ratio=odds_ratio, fisher_p_two_sided=pvalue,
            fisher_q_bh=np.nan, test_evaluable=evaluable,
            test_status='evaluable' if evaluable else 'not_evaluable_zero_margin',
        ))
        assert tp + fp + fn + tn == total

    fisher = pd.DataFrame(fisher_rows)
    evaluable = fisher['test_evaluable']
    fisher.loc[evaluable, 'fisher_q_bh'] = benjamini_hochberg(fisher.loc[evaluable, 'fisher_p_two_sided'])
    fisher['bh_family_size'] = int(evaluable.sum())
    genes = pd.DataFrame(gene_rows, columns=GENE_COLUMNS)
    tables = {
        'antifungal_resistance_explained_confusion_metrics.csv': pd.DataFrame(metric_rows),
        'antifungal_resistance_explained_gene_contributions.csv': genes,
        'Figure_2_phylogeny_and_FungAMR_panel_B_verified_cells.csv': pd.DataFrame(cells),
        'Figure_2_phylogeny_and_FungAMR_panel_B_verified_genes.csv': genes[GENE_COLUMNS[:5]],
        'panel_B_isolate_classifications.csv': pd.DataFrame(isolate_rows, columns=['sample', 'antifungal', 'phenotype', 'mutation_present']),
        'panel_B_cohort_selection.csv': cohort.reset_index(),
        'panel_B_fisher_enrichment.csv': fisher,
    }
    audit = dict(
        excluded_studies=exclusions, phenotype_cohort_n=len(cohort),
        excluded_cohort_n=int(cohort['excluded_by_study'].sum()),
        retained_cohort_n=int(cohort['included_in_assessment'].sum()),
        excluded_by_study=cohort.loc[cohort['excluded_by_study'], 'Study_clean'].value_counts().sort_index().to_dict(),
        panel_B_isolates_before=int(cohort['has_any_nonmissing_phenotype'].sum()),
        panel_B_isolates_after=int(cohort['included_in_panel_B'].sum()),
        panel_B_observations_before=int(cohort['n_nonmissing_drug_phenotypes'].sum()),
        panel_B_observations_after=len(isolate_rows),
        input_carrier_rows=raw_match_rows, qualifying_carrier_rows=len(matches),
        study_selection='Exact Study_clean labels; exclusions affect assessment only',
        carrier_rule='At least one qualifying drug-matched substitution; one record per isolate and drug',
        signed_evidence_rule=SIGNED_EVIDENCE_RULE,
        phenotype_rule='R or 1 is resistant; S or 0 is susceptible; all other values omitted separately for each drug',
        absence_rule='No detected qualifying substitution; missing genotype calls do not count as carriers',
        heatmap_denominators='Mutation-present and mutation-absent isolates with an available phenotype',
        gene_bar_denominator='All resistant isolates with the corresponding phenotype in the selected cohort',
        fisher_alternative='two-sided',
        fisher_table_order={'rows': ['mutation_present', 'mutation_absent'], 'columns': ['resistant', 'susceptible']},
        odds_ratio_definition='Sample odds ratio (present_R * absent_S) / (present_S * absent_R), as returned by scipy.stats.fisher_exact',
        non_evaluable_rule='Any zero row or column margin: odds ratio, p and q are unavailable',
        multiple_testing='Benjamini-Hochberg across evaluable drugs only',
        bh_family_size=int(evaluable.sum()),
        interpretation='See Methods, FungAMR annotation and phenotypic concordance.',
        versions={'numpy': np.__version__, 'pandas': pd.__version__, 'scipy': scipy.__version__},
        inputs={str(path): sha(path) for path in [matches_path, phenotypes_path, metadata_path]},
    )
    return tables, audit


def run_assessment(matches_path, phenotypes_path, metadata_path, exclude_studies,
                   outdir, sample_ids=None):
    tables, audit = build_assessment(matches_path, phenotypes_path, metadata_path,
                                     exclude_studies, sample_ids)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    for name, table in tables.items():
        table.to_csv(outdir / name, index=False, na_rep='NA' if name == 'panel_B_fisher_enrichment.csv' else '')
    audit['output_sha256'] = {name: sha(outdir / name) for name in tables}
    (outdir / 'panel_B_assessment_audit.json').write_text(json.dumps(audit, indent=2) + '\n')
    return tables, audit


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--matches', type=Path, required=True)
    parser.add_argument('--phenotypes', type=Path, required=True)
    parser.add_argument('--metadata', type=Path, required=True)
    parser.add_argument('--exclude-studies', nargs='*', default=[])
    parser.add_argument('--outdir', type=Path, required=True)
    args = parser.parse_args()
    tables, audit = run_assessment(args.matches, args.phenotypes, args.metadata,
                                   args.exclude_studies, args.outdir)
    print(tables['panel_B_fisher_enrichment.csv'].to_string(index=False))
    print(json.dumps({key: audit[key] for key in [
        'excluded_studies', 'excluded_cohort_n', 'retained_cohort_n',
        'panel_B_isolates_after', 'panel_B_observations_after', 'bh_family_size',
    ]}, indent=2))


if __name__ == '__main__':
    main()
