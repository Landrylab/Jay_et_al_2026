"""Compare a study-input rerun with the deposited statistical reference tables."""
from pathlib import Path
import argparse
import hashlib
import json

import pandas as pd
from pandas.testing import assert_frame_equal


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, required=True)
    parser.add_argument('--results-dir', type=Path, required=True)
    args = parser.parse_args()
    manifest = json.loads((args.input_dir/'manifest.json').read_text())
    for relative, checksum in manifest['sha256'].items():
        assert hashlib.sha256((args.input_dir/relative).read_bytes()).hexdigest() == checksum, relative
    reference = args.input_dir/'reference_results'
    comparisons = []
    for stage, name, keys in [
        ('gene_screen','gene_results.csv',['drug','gene_id']),
        ('gene_screen','pool_results.csv',['drug','pool']),
        ('cocarriage','conditional_results.csv',['model_id','method']),
        ('variants','variant_results.csv',['key']),
    ]:
        previous = pd.read_csv(reference/name).sort_values(keys).reset_index(drop=True)
        current = pd.read_csv(args.results_dir/stage/name).sort_values(keys).reset_index(drop=True)
        assert not current.duplicated(keys).any()
        assert_frame_equal(previous, current, check_dtype=False, check_exact=False, rtol=1e-9, atol=1e-12)
        comparisons.append(dict(stage=stage, file=name, rows=len(current), columns=len(current.columns), matched=True))
    genes = pd.read_csv(args.results_dir/'gene_screen/gene_results.csv')
    retained = pd.read_csv(args.results_dir/'variants/retained_gene_associations.csv')
    variants = pd.read_csv(args.results_dir/'variants/variant_results.csv')
    prioritized = pd.read_csv(args.results_dir/'variants/prioritized_variants.csv')
    counts = dict(primary_gene_drug_tests=len(genes),
                  qualifying_gene_drug_associations=int((genes.pca_cluster_q_primary.lt(.05) & genes.pca_cluster_or.gt(1)).sum()),
                  retained_gene_drug_associations=len(retained), all_variant_tests=len(variants), prioritized_variants=len(prioritized))
    assert counts == manifest['expected_results'], (counts, manifest['expected_results'])
    report = dict(status='passed', input_hashes_verified=len(manifest['sha256']), results=counts,
                  comparison_tolerance=dict(relative=1e-9, absolute=1e-12), tables=comparisons)
    (args.results_dir/'validation.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
