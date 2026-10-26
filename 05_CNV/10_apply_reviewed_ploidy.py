#!/usr/bin/env python3
"""Apply documented ploidy reviews and summarize accepted cohort assignments.

This changes the baseline metadata only. Refit affected sample calls with
04_segment_and_assign_copy_number.R before collating copy-number tracks.
"""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path


ALLOWED_PLOIDIES = {2, 3, 4, 6}


def read_rows(path, required):
    with Path(path).open(newline='', encoding='utf-8-sig') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        if not set(required).issubset(reader.fieldnames or []):
            raise ValueError(f'{path}: required columns are {", ".join(required)}')
        rows = list(reader)
    if not rows:
        raise ValueError(f'{path}: empty table')
    seen = set()
    for row in rows:
        if any(not (row.get(key) or '').strip() for key in required):
            raise ValueError(f'{path}: blank required value')
        sample = row['sample_id']
        if sample in seen:
            raise ValueError(f'{path}: duplicate sample {sample}')
        seen.add(sample)
        try:
            numeric = float(row['baseline_ploidy'])
            ploidy = int(numeric)
        except (ValueError, OverflowError):
            raise ValueError(f'{sample}: invalid baseline ploidy') from None
        if numeric != ploidy or ploidy not in ALLOWED_PLOIDIES:
            raise ValueError(f'{sample}: baseline must be one of {sorted(ALLOWED_PLOIDIES)}')
        row['baseline_ploidy'] = ploidy
    return rows


def apply_reviews(calls, reviews):
    accepted = {row['sample_id']: row['baseline_ploidy'] for row in calls}
    unknown = {row['sample_id'] for row in reviews} - set(accepted)
    if unknown:
        raise ValueError(f'Review contains unknown samples: {sorted(unknown)}')
    changes = []
    for row in reviews:
        sample = row['sample_id']
        original = accepted[sample]
        accepted[sample] = row['baseline_ploidy']
        changes.append(dict(row, input_baseline_ploidy=original,
                            changed=original != row['baseline_ploidy']))
    return accepted, changes


def write_table(path, fields, rows, delimiter='\t'):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter=delimiter,
                                lineterminator='\n', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ploidy', required=True, type=Path,
                        help='TSV with sample_id and baseline_ploidy')
    parser.add_argument('--review', required=True, type=Path,
                        help='TSV with sample_id, baseline_ploidy, reviewer and evidence')
    parser.add_argument('--output-dir', required=True, type=Path)
    args = parser.parse_args()
    calls = read_rows(args.ploidy, ['sample_id', 'baseline_ploidy'])
    reviews = read_rows(args.review, ['sample_id', 'baseline_ploidy', 'reviewer', 'evidence'])
    accepted, changes = apply_reviews(calls, reviews)
    out = args.output_dir
    output_names = ['baseline_ploidy.tsv', 'applied_reviews.tsv', 'ploidy_counts.csv', 'ploidy_review_audit.json']
    inputs = {args.ploidy.resolve(), args.review.resolve()}
    if any((out / name).resolve() in inputs for name in output_names):
        raise ValueError('Choose an output directory that does not overwrite an input table.')
    out.mkdir(parents=True, exist_ok=True)
    write_table(out/'baseline_ploidy.tsv', ['sample_id', 'baseline_ploidy'],
                [{'sample_id': sample, 'baseline_ploidy': value} for sample, value in accepted.items()])
    write_table(out/'applied_reviews.tsv', ['sample_id', 'input_baseline_ploidy', 'baseline_ploidy',
                'changed', 'reviewer', 'evidence'], changes)
    counts = Counter(accepted.values())
    total = len(accepted)
    write_table(out/'ploidy_counts.csv', ['ploidy', 'isolates', 'cohort_total', 'percentage'],
                [{'ploidy': f'{ploidy}N', 'isolates': count, 'cohort_total': total,
                  'percentage': 100 * count / total} for ploidy, count in sorted(counts.items())], ',')
    audit = {'cohort_size': total, 'ploidy_counts': dict(sorted(counts.items())),
             'reviewed_samples': len(reviews),
             'changed_samples': [r['sample_id'] for r in changes if r['changed']],
             'requires_copy_number_reassignment': [r['sample_id'] for r in changes if r['changed']],
             'input_sha256': {key: hashlib.sha256(path.read_bytes()).hexdigest()
                              for key, path in [('ploidy', args.ploidy), ('review', args.review)]}}
    (out/'ploidy_review_audit.json').write_text(json.dumps(audit, indent=2)+'\n')
    print(json.dumps(audit, indent=2))


if __name__ == '__main__':
    main()
