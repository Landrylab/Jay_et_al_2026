#!/usr/bin/env python3
"""Replay portable Figure 2b inputs and optionally compare supplied source tables."""
import argparse
import hashlib
import importlib.util
import json
import math
import tempfile
import zipfile
from pathlib import Path

import pandas as pd
from fungamr_evidence import best_signed_score, is_excluded_substitution


def load_inputs(source, destination):
    """Copy only declared input members, without extracting arbitrary paths."""
    names = ('input_manifest.json', 'matches.csv', 'phenotypes.tsv')
    if source.is_dir():
        contents = {name: (source / name).read_bytes() for name in names}
    else:
        with zipfile.ZipFile(source) as archive:
            if len(archive.namelist()) != len(set(archive.namelist())):
                raise ValueError('Input archive contains duplicate member names')
            contents = {name: archive.read(name) for name in names}
    manifest = json.loads(contents['input_manifest.json'])
    for name in names[1:]:
        expected = manifest.get('members', {}).get(name)
        if not expected or hashlib.sha256(contents[name]).hexdigest() != expected:
            raise ValueError(f'Input checksum does not match manifest: {name}')
    if manifest.get('confidence_scores') != [1, 2, 3, 4]:
        raise ValueError('This validator requires the current Figure 2 confidence configuration')
    for name, data in contents.items():
        (destination / name).write_bytes(data)
    return manifest


def compare_tables(actual, reference):
    if list(actual.columns) != list(reference.columns):
        raise AssertionError('Reference table columns or column order differ')
    keys = list(actual.columns)
    pd.testing.assert_frame_equal(
        actual.sort_values(keys, na_position='last').reset_index(drop=True),
        reference.sort_values(keys, na_position='last').reset_index(drop=True),
        check_dtype=False, check_exact=False, rtol=1e-12, atol=1e-14)


def validate_fisher(table):
    """Enumerate the hypergeometric probabilities independently of SciPy."""
    evaluated = []
    for row in table.itertuples(index=False):
        a, b, c, d = (int(getattr(row, name)) for name in (
            'mutation_present_resistant_n', 'mutation_present_susceptible_n',
            'mutation_absent_resistant_n', 'mutation_absent_susceptible_n'))
        first, second, total_r = a + b, c + d, a + c
        if min(first, second, total_r, b + d) == 0:
            if not (pd.isna(row.fisher_p_two_sided) and pd.isna(row.fisher_q_bh)):
                raise AssertionError('A non-evaluable test has reported probabilities')
            continue
        denominator = math.comb(first + second, total_r)
        observed = math.comb(first, a) * math.comb(second, total_r - a) / denominator
        probabilities = [math.comb(first, x) * math.comb(second, total_r - x) / denominator
                         for x in range(max(0, total_r - second), min(first, total_r) + 1)]
        pvalue = sum(p for p in probabilities if p <= observed * (1 + 1e-12))
        if not math.isclose(pvalue, row.fisher_p_two_sided, rel_tol=1e-10, abs_tol=1e-15):
            raise AssertionError(f'Fisher probability differs for {row.antifungal}')
        evaluated.append((row.antifungal, pvalue, row.fisher_q_bh))
    family_size = len(evaluated)
    if not table.bh_family_size.eq(family_size).all():
        raise AssertionError('Reported BH family size differs from evaluable tests')
    running = 1.0
    for rank, (drug, pvalue, qvalue) in reversed(list(enumerate(sorted(evaluated, key=lambda x: x[1]), 1))):
        running = min(running, pvalue * family_size / rank)
        if not math.isclose(running, qvalue, rel_tol=1e-10, abs_tol=1e-15):
            raise AssertionError(f'BH-adjusted probability differs for {drug}')
    return family_size


def validate(args, work):
    manifest = load_inputs(args.inputs, work)
    spec = importlib.util.spec_from_file_location('figure2_assessment', Path(__file__).with_name('02_assess_fungamr.py'))
    assessment = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(assessment)
    tables, audit = assessment.run_assessment(
        work / 'matches.csv', work / 'phenotypes.tsv', work / 'phenotypes.tsv',
        manifest.get('cohort_exclusions', []), work / 'results',
        confidence_scores=manifest['confidence_scores'])
    carriers = pd.read_csv(work / 'matches.csv')
    if len(carriers) != manifest['full_cohort_carrier_rows']:
        raise AssertionError('Carrier row count differs from input manifest')
    if not carriers.confidence_score.isin(manifest['confidence_scores']).all():
        raise AssertionError('Carrier input does not match its confidence configuration')
    if any(is_excluded_substitution(g, m) for g, m in carriers[['gene', 'mutation_found']].itertuples(index=False, name=None)):
        raise AssertionError('Carrier input does not match the Methods-defined marker configuration')
    assert best_signed_score([8, -2]) == -2 and best_signed_score([-2, 2]) == 2
    assert is_excluded_substitution('ERG11', 'D116E') and is_excluded_substitution('ERG11', 'V488I')
    assert not is_excluded_substitution('FKS1', 'D116E')
    classifications = tables['panel_B_isolate_classifications.csv']
    cells = tables['Figure_2_phylogeny_and_FungAMR_panel_B_verified_cells.csv']
    if set(cells.phenotype) != {'R phenotype', 'S phenotype'}:
        raise AssertionError('Unexpected displayed phenotype classes')
    if classifications.duplicated(['sample', 'antifungal']).any():
        raise AssertionError('Duplicate isolate–drug classifications')
    for drug, group in cells.groupby('antifungal'):
        if int(group['count'].sum()) != int(classifications.antifungal.eq(drug).sum()):
            raise AssertionError(f'Heatmap counts do not reconcile for {drug}')
    family_size = validate_fisher(tables['panel_B_fisher_enrichment.csv'])
    compared = []
    if args.reference_dir is not None:
        for name in tables:
            target = args.reference_dir / name
            if not target.is_file():
                raise FileNotFoundError(f'Missing reference assessment table: {target}')
            compare_tables(pd.read_csv(work / 'results' / name), pd.read_csv(target))
            compared.append(name)
        catalogue = args.reference_dir / 'fungamr_high_confidence_carriers.csv'
        if catalogue.exists():
            compare_tables(carriers, pd.read_csv(catalogue))
            compared.append(catalogue.name)
    return {
        'status': 'passed', 'input_manifest_checksums_verified': True,
        'full_cohort_carrier_rows': len(carriers),
        'assessment_isolates': int(classifications['sample'].nunique()),
        'isolate_drug_observations': len(classifications), 'heatmap_cells': len(cells),
        'gene_bars': len(tables['Figure_2_phylogeny_and_FungAMR_panel_B_verified_genes.csv']),
        'fisher_probabilities_independently_enumerated': True, 'bh_family_size': family_size,
        'reference_tables_compared': compared, 'full_master_matching_rerun': False,
        'assessment_input_sha256': manifest['members'], 'output_sha256': audit['output_sha256']}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, default=Path(__file__).with_name('input_data.zip'),
                        help='Portable input ZIP or directory with its manifest and two input tables')
    parser.add_argument('--reference-dir', type=Path,
                        help='Optional directory of current Figure 2b source CSVs for table comparisons')
    parser.add_argument('--report', type=Path, help='Optional JSON report; stdout is always written')
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='figure2-validation-') as temporary:
        result = validate(args, Path(temporary))
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
