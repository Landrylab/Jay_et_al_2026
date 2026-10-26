#!/usr/bin/env python3
"""Compare sampled master target-ALT dosages with original VCF/BCF genotypes."""
import argparse
import csv
import gzip
import hashlib
import json
import math
import re
import subprocess
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np


def stream_query(command):
    """Stream query output and propagate HTSlib failures without a stderr pipe deadlock."""
    with tempfile.TemporaryFile(mode='w+t') as log:
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=log, text=True)
        try:
            for line in process.stdout:
                yield line.rstrip('\n').split('\t')
            if process.wait():
                log.seek(0)
                raise RuntimeError('BCFtools query failed: ' + log.read().strip())
        finally:
            process.stdout.close()
            if process.poll() is None:
                process.terminate()
                process.wait()


def target_dosage(genotype, alt_index, n_alts):
    alleles = re.split(r'[/|]', genotype)
    if len(alleles) != 2 or any(a == '.' for a in alleles):
        return -1
    if any(not a.isdigit() for a in alleles):
        raise ValueError(f'Invalid genotype: {genotype}')
    numbers = [int(a) for a in alleles]
    if any(a > n_alts for a in numbers):
        raise ValueError(f'Genotype refers to an absent ALT: {genotype}')
    return sum(a == alt_index for a in numbers)


def read_tsv(path):
    with gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path) as stream:
        yield from csv.DictReader(stream, delimiter='\t')


def source_number(name):
    match = re.fullmatch(r'src(\d+)', name)
    if not match or int(match.group(1)) < 1:
        raise ValueError(f'Invalid original-record identifier: {name}')
    return int(match.group(1))


def file_hash(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def run(args):
    paths = {name: Path(str(args.prefix) + suffix) for name, suffix in {
        'dosage': '.dosage.npy', 'samples': '.samples.tsv', 'variants': '.variants.tsv.gz',
        'mapping': '.normalization_map.tsv.gz', 'original_records': '.original_records.tsv.gz',
        'original_genotypes': '.original_genotypes.bcf'}.items()}
    if args.original_vcf:
        paths['original_genotypes'] = args.original_vcf
    for path in paths.values():
        if not path.is_file():
            raise FileNotFoundError(path)
    if args.out.resolve() in {p.resolve() for p in paths.values()}:
        raise ValueError('Report path must not overwrite a source file')
    samples = [r['sample'] for r in read_tsv(paths['samples'])]
    if not samples or len(samples) != len(set(samples)) or any(not s for s in samples):
        raise ValueError('Sample IDs must be present and unique')
    original_samples = subprocess.check_output(
        [args.bcftools, 'query', '-l', str(paths['original_genotypes'])], text=True).splitlines()
    if samples != original_samples:
        raise ValueError('Original genotype sample order differs from dosage columns')
    matrix = np.load(paths['dosage'], mmap_mode='r', allow_pickle=False)
    if matrix.ndim != 2 or matrix.shape[1] != len(samples) or matrix.dtype != np.int8:
        raise ValueError('Dosage must be a two-dimensional int8 matrix matching the sample table')
    if not matrix.shape[0]:
        raise ValueError('Cannot audit an empty catalogue')
    n = min(args.rows, matrix.shape[0])
    selected = set(map(int, np.random.default_rng(args.seed).choice(matrix.shape[0], n, replace=False)))
    metadata = {}
    row_count = 0
    for row_count, row in enumerate(read_tsv(paths['variants']), 1):
        index = int(row['index'])
        if index != row_count - 1:
            raise ValueError('Variant indices must be contiguous and match matrix row order')
        if index not in selected:
            continue
        values = matrix[index]
        if not np.isin(values, [-1, 0, 1, 2]).all():
            raise ValueError(f'Invalid dosage value at row {index}')
        called = values >= 0
        observed = {'ac': int(values[called].sum()), 'an': int(called.sum() * 2),
                    'carrier_count': int((values > 0).sum()), 'missing_fraction': float((~called).mean())}
        for key, value in observed.items():
            if not math.isclose(float(row[key]), value, rel_tol=1e-12, abs_tol=1e-15):
                raise ValueError(f'Metadata {key} disagrees with dosage at row {index}')
        if observed['an'] and not math.isclose(float(row['af']), observed['ac']/observed['an'], rel_tol=1e-12, abs_tol=1e-15):
            raise ValueError(f'Metadata allele frequency disagrees at row {index}')
        metadata[index] = row
    if row_count != matrix.shape[0] or len(metadata) != n:
        raise ValueError('Variant metadata dimensions differ from the dosage matrix')
    mappings = defaultdict(list)
    for row in read_tsv(paths['mapping']):
        index = int(row['index'])
        if index in selected:
            if row['id'] != metadata[index]['id']:
                raise ValueError(f'Normalization-map identifier disagrees at row {index}')
            pair = (source_number(row['source_record']), int(row['source_alt_index']))
            if pair in mappings[index]:
                raise ValueError(f'Duplicate source mapping for row {index}')
            mappings[index].append(pair)
    if set(mappings) != selected:
        raise ValueError('Selected variants are missing normalization mappings')
    wanted = {pair for pairs in mappings.values() for pair in pairs}
    originals = {}
    for row in read_tsv(paths['original_records']):
        pair = (source_number(row['source_record']), int(row['alt_index']))
        if pair in wanted:
            if pair in originals:
                raise ValueError(f'Duplicate original-record mapping: {pair}')
            alts = row['alts'].split(',')
            if pair[1] < 1 or pair[1] > len(alts) or row['target_alt'] != alts[pair[1]-1]:
                raise ValueError(f'Invalid original ALT index: {pair}')
            originals[pair] = (row['chrom'], row['pos'], row['original_id'], row['ref'], row['alts'])
    if set(originals) != wanted:
        raise ValueError('Original-record table is missing selected allele mappings')
    # First verify the one-based source ordinals using compact site-only output.
    expected_sites = {}
    for (ordinal, _), key in originals.items():
        if ordinal in expected_sites and expected_sites[ordinal] != key:
            raise ValueError(f'Conflicting site metadata for source ordinal {ordinal}')
        expected_sites[ordinal] = key
    seen_ordinals = set()
    fmt = '%CHROM\t%POS\t%ID\t%REF\t%ALT\n'
    for ordinal, fields in enumerate(stream_query([args.bcftools, 'query', '-f', fmt, str(paths['original_genotypes'])]), 1):
        if ordinal in expected_sites:
            if tuple(fields) != expected_sites[ordinal]:
                raise ValueError(f'Source ordinal {ordinal} does not match the retained original record')
            seen_ordinals.add(ordinal)
    if seen_ordinals != set(expected_sites):
        raise ValueError('Selected source ordinals extend beyond the original genotype file')
    # Fetch genotypes only at sampled original positions through the existing index.
    targets = defaultdict(list)
    for index, pairs in mappings.items():
        for pair in pairs:
            targets[originals[pair]].append((index, pair))
    counts = Counter()
    compared_pairs = set()
    mismatches = []
    with tempfile.TemporaryDirectory(prefix='master-dosage-audit-') as temporary:
        regions = Path(temporary) / 'positions.tsv'
        coordinates = sorted({(key[0], int(key[1])) for key in targets})
        regions.write_text(''.join(f'{chrom}\t{position}\t{position}\n' for chrom, position in coordinates))
        command = [args.bcftools, 'query', '-R', str(regions), '-f', '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n', str(paths['original_genotypes'])]
        for fields in stream_query(command):
            key = tuple(fields[:5])
            if key not in targets:
                continue
            if len(fields[5:]) != len(samples):
                raise ValueError('Source genotype count does not match the sample table')
            for index, pair in targets[key]:
                expected = np.array([target_dosage(gt, pair[1], len(key[4].split(','))) for gt in fields[5:]], dtype=np.int8)
                differing = np.flatnonzero(expected != matrix[index])
                counts['dosage_disagreements'] += len(differing)
                for column in differing[:max(0, 10-len(mismatches))]:
                    mismatches.append({'index':index, 'sample':samples[int(column)], 'observed':int(matrix[index,column]), 'expected':int(expected[column]), 'source_record':pair[0], 'source_alt_index':pair[1]})
                # Repeated identical source records are checked, but not counted twice.
                if (index, pair) not in compared_pairs:
                    counts['genotype_comparisons'] += len(samples)
                    counts['missing_comparisons'] += int((expected < 0).sum())
                    counts['multiallelic_source_alleles'] += int(len(key[4].split(',')) > 1)
                    compared_pairs.add((index, pair))
    if compared_pairs != {(i,p) for i, pairs in mappings.items() for p in pairs}:
        raise ValueError('Indexed source lookup did not recover every selected allele mapping')
    report = {'status':'passed' if not counts['dosage_disagreements'] else 'failed',
              'sampled_rows':n, 'seed':args.seed, 'matrix_shape':list(matrix.shape),
              'matrix_dtype':str(matrix.dtype), 'sample_order_verified':True,
              'sampled_metadata_verified':True, 'source_record_ordinals_verified':True,
              'all_selected_source_mappings_checked':True, 'counts':dict(counts),
              'dosage_disagreements':counts['dosage_disagreements'], 'mismatch_examples':mismatches,
              'inputs':{name:{'file':path.name,'sha256':file_hash(path)} for name,path in paths.items()},
              'bcftools_version':subprocess.check_output([args.bcftools,'--version'],text=True).splitlines()[0]}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report,indent=2)+'\n')
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prefix', type=Path, required=True, help='Shared catalogue prefix')
    parser.add_argument('--out', type=Path, required=True, help='JSON validation report')
    parser.add_argument('--rows', type=int, default=4000, help='Number of deterministic sampled rows (default: 4000; capped at catalogue size)')
    parser.add_argument('--seed', type=int, default=20260923)
    parser.add_argument('--original-vcf', type=Path, help='Indexed original cohort VCF/BCF; defaults to <prefix>.original_genotypes.bcf')
    parser.add_argument('--bcftools', default='bcftools', help='BCFtools executable')
    args = parser.parse_args()
    if args.rows < 1 or args.seed < 0:
        parser.error('--rows must be positive and --seed must be nonnegative')
    report = run(args)
    print(json.dumps({k:v for k,v in report.items() if k != 'inputs'},indent=2))
    if report['status'] != 'passed':
        raise SystemExit(1)


if __name__ == '__main__':
    main()
