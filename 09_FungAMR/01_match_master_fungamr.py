#!/usr/bin/env python3
"""Match allele-specific protein substitutions to high-confidence FungAMR records."""
import argparse
import csv
import gzip
import hashlib
import importlib.util
import json
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

AA3 = dict(zip('Ala Arg Asn Asp Cys Gln Glu Gly His Ile Leu Lys Met Phe Pro Ser Thr Trp Tyr Val Ter Sec Pyl'.split(), 'ARNDCQEGHILKMFPSTWYV*UO'))
PROTEIN_EFFECTS = {'missense_variant', 'stop_gained', 'stop_lost', 'start_lost'}
GENE_ALIASES = {'CYP51': 'ERG11', 'CYP51A': 'ERG11', 'FKS': 'FKS1', 'GSC1': 'FKS1', 'SQUALENE EPOXIDASE': 'ERG1'}
EXCLUDED_MUTATIONS = {'V488I'}
EXCLUDED_PAIRS = {('FKS1', 'Fluconazole')}


def short_protein_change(value):
    m = re.fullmatch(r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})', value or '')
    return AA3[m[1]] + m[2] + AA3[m[3]] if m and m[1] in AA3 and m[3] in AA3 else ''


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--master-prefix', type=Path, required=True)
    parser.add_argument('--gene-map', type=Path, required=True)
    parser.add_argument('--fungamr', type=Path, required=True)
    parser.add_argument('--phenotypes', type=Path, required=True)
    parser.add_argument('--metadata', type=Path, help='Study metadata; defaults to the phenotype table')
    parser.add_argument('--exclude-studies', nargs='*', default=[], help='Exact Study_clean labels excluded from assessment only')
    parser.add_argument('--outdir', type=Path, required=True)
    args = parser.parse_args()
    out = args.outdir
    out.mkdir(parents=True, exist_ok=True)
    prefix = str(args.master_prefix)
    samples = pd.read_csv(prefix + '.samples.tsv', sep='\t')['sample'].astype(str).tolist()
    dosage = np.load(prefix + '.dosage.npy', mmap_mode='r')
    assert dosage.shape[1] == len(samples) and len(samples) == len(set(samples))
    genes = pd.read_csv(args.gene_map, sep='\t').fillna('')
    by_symbol = {row.gene_symbol.upper(): row.gene_id for row in genes.itertuples()}
    gene_symbol = genes.set_index('gene_id')['gene_symbol'].to_dict()
    database = pd.read_csv(args.fungamr, sep=None, engine='python', dtype=str).fillna('')
    lookup = defaultdict(list)
    db_rows = []
    for row in database.to_dict('records'):
        score = pd.to_numeric(row['confidence score'], errors='coerce')
        mutation = row['mutation'].strip().upper()
        if row['species'].strip().upper() != 'CANDIDA ALBICANS' or score not in {1, 2}:
            continue
        if row['mutation_type'].upper() not in {'SNP', 'SNPS'} or not re.fullmatch(r'[A-Z*]\d+[A-Z*]', mutation):
            continue
        symbols = {GENE_ALIASES.get(v.strip().upper(), v.strip().upper()) for field in ['gene or protein name', 'ortho_homolog'] for v in re.split(r'[,;/&]+', row[field])}
        gene_ids = {by_symbol[v] for v in symbols if v in by_symbol}
        if not gene_ids:
            raise ValueError(f"Unmapped FungAMR gene: {row['gene or protein name']} / {row['ortho_homolog']}")
        for gene_id in gene_ids:
            symbol = gene_symbol[gene_id]
            if mutation in EXCLUDED_MUTATIONS or (symbol, row['drug']) in EXCLUDED_PAIRS:
                continue
            record = dict(gene_id=gene_id, gene=symbol, mutation_found=mutation, confidence_score=int(score), antifungal=row['drug'], fungamr_gene=row['gene or protein name'], fungamr_id=row['ID'], pubmedid=row['pubmedid'], first_author=row['first author name'], year=row['year'], strongest_resistance_evidence=row['Strongest resistance evidence reported'])
            lookup[(gene_id, mutation)].append(record)
            db_rows.append(record)
    pd.DataFrame(db_rows).to_csv(out / 'fungamr_high_confidence_reference_records.csv', index=False)
    candidates = []
    seen = set()
    with gzip.open(prefix + '.gene_annotations.tsv.gz', 'rt') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            if row['gene_id'] not in gene_symbol:
                continue
            mutation = short_protein_change(row['hgvs_p'])
            matches = lookup.get((row['gene_id'], mutation))
            if not matches or not set(row['effect'].split('&')) & PROTEIN_EFFECTS:
                continue
            chrom, pos, ref, alt = row['id'].split(':')
            assert row['annotation_allele'] == alt
            for entry in matches:
                key = (row['id'], row['gene_id'], row['feature_id'], entry['fungamr_id'])
                if key in seen:
                    continue
                seen.add(key)
                candidates.append({**entry, 'variant_index': int(row['index']), 'variant_id': row['id'], 'chrom': chrom, 'pos': int(pos), 'ref': ref, 'alt': alt, 'effect': row['effect'], 'transcript_id': row['feature_id'], 'hgvs_c': row['hgvs_c'], 'hgvs_p': row['hgvs_p']})
    candidate_df = pd.DataFrame(candidates)
    candidate_df.to_csv(out / 'master_fungamr_allele_matches.csv', index=False)
    carriers = []
    callability = []
    for entry in candidates:
        values = dosage[entry['variant_index']]
        callability.append({**entry, 'n_called': int((values >= 0).sum()), 'n_missing': int((values < 0).sum()), 'n_carriers': int((values > 0).sum())})
        for s in np.flatnonzero(values > 0):
            carriers.append({**entry, 'sample': samples[s], 'alternate_allele_dosage': int(values[s])})
    pd.DataFrame(callability).to_csv(out / 'master_fungamr_allele_callability.csv', index=False)
    matches = pd.DataFrame(carriers).sort_values(['sample', 'gene', 'mutation_found', 'antifungal', 'confidence_score', 'fungamr_id']).drop_duplicates(['sample', 'variant_id', 'gene', 'mutation_found', 'antifungal'])
    matches.to_csv(out / 'fungamr_score_1_or_2_nonsynonymous_matches.csv', index=False)
    # Full catalogue matching above is independent of study selection. Only the
    # phenotype assessment below uses the requested cohort exclusions.
    assessment_path = Path(__file__).with_name('02_assess_fungamr.py')
    spec = importlib.util.spec_from_file_location('fungamr_assessment', assessment_path)
    assessment = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(assessment)
    tables, assessment_audit = assessment.run_assessment(
        out / 'fungamr_score_1_or_2_nonsynonymous_matches.csv', args.phenotypes,
        args.metadata or args.phenotypes, args.exclude_studies, out,
        sample_ids=samples,
    )
    audit = dict(master_shape=list(dosage.shape), high_confidence_database_records=len(db_rows), matched_alleles=int(candidate_df.variant_id.nunique()), matched_protein_changes=int(candidate_df[['gene','mutation_found']].drop_duplicates().shape[0]), carrier_drug_rows=len(matches), carrier_rule='Called dosage > 0 for the specified alternate allele', confidence_scores=[1,2], mutation_scope='Exact single amino-acid substitutions in protein-changing allele-specific annotations', excluded_mutations=sorted(EXCLUDED_MUTATIONS), excluded_gene_drug_pairs=[list(v) for v in sorted(EXCLUDED_PAIRS)], absence_rule='No detected drug-matched high-confidence substitution; missing calls do not count as carriers', inputs={str(v):sha(v) for v in [args.fungamr,args.gene_map,args.phenotypes,Path(prefix+'.gene_annotations.tsv.gz'),Path(prefix+'.samples.tsv')]})
    manifest = args.master_prefix.parent / 'master_SHA256SUMS'
    if manifest.exists():
        checksums = {line.split(None, 1)[1].strip(): line.split(None, 1)[0] for line in manifest.read_text().splitlines() if line.strip()}
        names = [args.master_prefix.name + suffix for suffix in ['.dosage.npy', '.vcf.gz']]
        audit['master_checksums'] = {name: checksums[name] for name in names if name in checksums}
        audit['master_checksums_source'] = str(manifest)
    (out / 'master_fungamr_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(tables['panel_B_fisher_enrichment.csv'].to_string(index=False))
    print(json.dumps({k:v for k,v in audit.items() if k!='inputs'},indent=2))

if __name__ == '__main__':
    main()
