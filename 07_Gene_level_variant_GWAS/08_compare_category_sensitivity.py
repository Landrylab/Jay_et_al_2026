#!/usr/bin/env python3
"""Compare quantitative gene GWAS after MIC-defined I/SDD-range exclusions."""
from __future__ import annotations

import argparse
import hashlib
import html
import json
from pathlib import Path

import numpy as np
import pandas as pd

DRUGS = ['fluconazole', 'voriconazole', 'caspofungin', 'micafungin', 'anidulafungin']
MODES = ['fixed', 'rebuilt']
MODE_LABELS = {'fixed': 'Original variants per gene retained', 'rebuilt': '5% rare-variant cutoff reapplied'}
EXPECTED = {'fluconazole': (871, 839), 'voriconazole': (659, 554),
            'caspofungin': (766, 758), 'micafungin': (667, 665), 'anidulafungin': (589, 585),
            'posaconazole': (611, 611)}
RANGES = {'fluconazole': '2 < MIC ≤ 4', 'voriconazole': '0.125 < MIC ≤ 0.5',
          'caspofungin': '0.25 < MIC ≤ 0.5', 'micafungin': '0.25 < MIC ≤ 0.5',
          'anidulafungin': '0.25 < MIC ≤ 0.5'}
KEY_SYMBOLS = {'ERG11', 'TAC1', 'UPC2', 'ERG3', 'FKS1', 'ALI1', 'INO80', 'ATP8', 'C5_00080C', 'C1_14090W'}
MEASURES = ['rank', 'p', 'neglog10P', 'carrier_vs_noncarrier_log2mic', 'carrier_vs_noncarrier_se',
            'n_carriers', 'n_noncarriers', 'n_missing', 'n_collapsed_variants', 'cohort_n',
            'tested_genes', 'bonferroni_p', 'suggestive_p', 'significance']


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def table(frame, path):
    frame.to_csv(path, sep='\t', index=False, na_rep='', float_format='%.17g')


def normalized_results(path, expected_drugs, mode):
    frame = pd.read_csv(path, sep='\t')
    required = {'drug', 'gene_id', 'gene_symbol', 'chromosome', 'bp', 'p', 'b', 'se', 'neglog10P',
                'cohort_n', 'n_carriers', 'n_noncarriers', 'n_missing', 'n_called', 'n_collapsed_variants',
                'carrier_vs_noncarrier_log2mic', 'carrier_vs_noncarrier_se'}
    if not required.issubset(frame):
        raise ValueError(f'Missing result columns in {path}: {sorted(required-set(frame))}')
    if set(frame.drug) != set(expected_drugs) or frame.duplicated(['drug', 'gene_id']).any():
        raise ValueError(f'Drug set or distinct gene identifiers differ in {path}')
    if not np.isfinite(frame[['p', 'b', 'se', 'neglog10P']]).all().all() or not frame.p.between(0, 1).all() or not frame.se.gt(0).all():
        raise ValueError(f'Invalid association statistic in {path}')
    if not (frame.n_carriers.gt(0) & frame.n_noncarriers.gt(0)).all():
        raise ValueError(f'Nonpolymorphic tested marker in {path}')
    if not (frame.n_carriers + frame.n_noncarriers).eq(frame.n_called).all() or not (frame.n_called + frame.n_missing).eq(frame.cohort_n).all():
        raise ValueError(f'Carrier/missing counts differ from cohort size in {path}')
    if not np.allclose(frame.carrier_vs_noncarrier_log2mic, 2 * frame.b, rtol=1e-10, atol=1e-12):
        raise ValueError(f'Carrier contrast differs from 2*b in {path}')
    if not np.allclose(frame.carrier_vs_noncarrier_se, 2 * frame.se, rtol=1e-10, atol=1e-12):
        raise ValueError(f'Carrier contrast standard error differs from 2*se in {path}')
    nonzero = frame.p.gt(0)
    if not np.allclose(frame.loc[nonzero, 'neglog10P'], -np.log10(frame.loc[nonzero, 'p']), rtol=1e-9, atol=1e-9):
        raise ValueError(f'Negative-log P values differ from P in {path}')
    for drug, part in frame.groupby('drug'):
        expected_n = EXPECTED[drug][0 if mode == 'original' else 1]
        if set(part.cohort_n) != {expected_n}:
            raise ValueError(f'{mode}/{drug}: expected cohort {expected_n}, received {part.cohort_n.unique()}')
    frame = frame.sort_values(['drug', 'p', 'chromosome', 'bp', 'gene_id']).reset_index(drop=True)
    frame['rank'] = frame.groupby('drug').cumcount() + 1
    frame['tested_genes'] = frame.groupby('drug').gene_id.transform('size')
    frame['bonferroni_p'] = .05 / frame.tested_genes
    frame['suggestive_p'] = 1. / frame.tested_genes
    frame['significance'] = np.select([frame.p.lt(frame.bonferroni_p), frame.p.lt(frame.suggestive_p)],
                                      ['Bonferroni', 'suggestive'], default='below_suggestive')
    for old, new in [('per_drug_bonferroni_p', 'bonferroni_p'), ('per_drug_suggestive_p', 'suggestive_p')]:
        if old in frame and not np.allclose(frame[old], frame[new], rtol=1e-10, atol=1e-15):
            raise ValueError(f'Test-count threshold differs from actual tested genes: {path}, {old}')
    frame['display_gene'] = frame.gene_symbol.fillna(frame.gene_id).astype(str)
    frame.loc[frame.display_gene.str.strip().eq(''), 'display_gene'] = frame.loc[frame.display_gene.str.strip().eq(''), 'gene_id']
    frame['mode'] = mode
    return frame


def comparisons(original, filtered):
    joined = []
    summaries = []
    for mode in MODES:
        for drug in DRUGS:
            a = original[original.drug.eq(drug)]
            b = filtered[mode][filtered[mode].drug.eq(drug)]
            meta = pd.concat([a[['gene_id', 'display_gene', 'chromosome', 'bp']],
                              b[['gene_id', 'display_gene', 'chromosome', 'bp']]]).drop_duplicates('gene_id')
            x = a[['gene_id'] + MEASURES].merge(b[['gene_id'] + MEASURES], on='gene_id', how='outer',
                                               suffixes=('_original', '_new'), validate='one_to_one').merge(meta, on='gene_id', validate='one_to_one')
            x.insert(0, 'drug', drug)
            x.insert(0, 'mode', mode)
            x['tested_original'] = x.rank_original.notna()
            x['tested_new'] = x.rank_new.notna()
            x['significance_original'] = x.significance_original.fillna('untested')
            x['significance_new'] = x.significance_new.fillna('untested')
            x['rank_improvement'] = x.rank_original - x.rank_new
            x['neglog10P_change'] = x.neglog10P_new - x.neglog10P_original
            x['carrier_effect_change'] = x.carrier_vs_noncarrier_log2mic_new - x.carrier_vs_noncarrier_log2mic_original
            x['effect_direction_agrees'] = np.where(x.tested_original & x.tested_new,
                np.sign(x.carrier_vs_noncarrier_log2mic_original).eq(np.sign(x.carrier_vs_noncarrier_log2mic_new)).astype(object), None)
            x['top20_original'] = x.rank_original.le(20)
            x['top20_new'] = x.rank_new.le(20)
            x['top20_transition'] = np.select([x.top20_original & x.top20_new, x.top20_new, x.top20_original],
                                               ['retained', 'entered', 'left'], default='outside_both')
            for criterion, sa, sb in [('bonferroni', x.significance_original.eq('Bonferroni'), x.significance_new.eq('Bonferroni')),
                                      ('suggestive_or_stronger', x.significance_original.isin(['Bonferroni', 'suggestive']), x.significance_new.isin(['Bonferroni', 'suggestive']))]:
                x[criterion + '_transition'] = np.select([sa & sb, ~sa & sb, sa & ~sb], ['retained', 'gained', 'lost'], default='below_both')
            common = x[x.tested_original & x.tested_new]
            row = {'mode': mode, 'drug': drug, 'original_n': int(a.cohort_n.iloc[0]), 'new_n': int(b.cohort_n.iloc[0]),
                   'removed_n': int(a.cohort_n.iloc[0] - b.cohort_n.iloc[0]),
                   'original_tested_genes': len(a), 'new_tested_genes': len(b), 'shared_tested_genes': len(common),
                   'original_bonferroni_p': .05 / len(a), 'new_bonferroni_p': .05 / len(b),
                   'original_suggestive_p': 1 / len(a), 'new_suggestive_p': 1 / len(b),
                   'original_top20_retained': int((x.top20_original & x.top20_new).sum()),
                   'original_top20_tested_new': int((x.top20_original & x.tested_new).sum()),
                   'top20_jaccard': float((x.top20_original & x.top20_new).sum() / (x.top20_original | x.top20_new).sum()),
                   'common_gene_rank_spearman': float(common.rank_original.corr(common.rank_new, method='spearman')),
                   'original_bonferroni_count': int(a.significance.eq('Bonferroni').sum()), 'new_bonferroni_count': int(b.significance.eq('Bonferroni').sum()),
                   'original_suggestive_only_count': int(a.significance.eq('suggestive').sum()), 'new_suggestive_only_count': int(b.significance.eq('suggestive').sum())}
            for criterion in ['bonferroni', 'suggestive_or_stronger']:
                for change in ['retained', 'gained', 'lost']:
                    subset = x[x[criterion + '_transition'].eq(change)]
                    row[criterion + '_' + change + '_count'] = len(subset)
                    row[criterion + '_' + change + '_genes'] = '; '.join(subset.sort_values('rank_new' if change != 'lost' else 'rank_original').display_gene)
            summaries.append(row)
            joined.append(x)
    combined = pd.concat(joined, ignore_index=True)
    summary = pd.DataFrame(summaries)
    return combined, summary


def score_figures(compare, out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9, 'axes.titlesize': 11,
                         'pdf.fonttype': 42, 'ps.fonttype': 42, 'axes.spines.right': False, 'axes.spines.top': False})
    palette = {'below_both': '#B9C0C5', 'retained': '#B22875', 'gained': '#008A8A', 'lost': '#DE9408'}
    for mode in MODES:
        fig, axes = plt.subplots(2, 3, figsize=(13, 8.3))
        fig.subplots_adjust(left=.075, right=.985, top=.91, bottom=.095, hspace=.30, wspace=.26)
        for ax, drug in zip(axes.flat, DRUGS):
            frame = compare[compare['mode'].eq(mode) & compare.drug.eq(drug) & compare.tested_original & compare.tested_new].copy()
            lim = max(7., float(frame[['neglog10P_original', 'neglog10P_new']].max().max()) * 1.10)
            for category in ['below_both', 'retained', 'lost', 'gained']:
                part = frame[frame.bonferroni_transition.eq(category)]
                ax.scatter(part.neglog10P_original, part.neglog10P_new,
                           s=5 if category == 'below_both' else 20, color=palette[category],
                           alpha=.45 if category == 'below_both' else .9, linewidths=0, rasterized=True)
            ax.plot([0, lim], [0, lim], color='#515B64', lw=.7, ls='--', zorder=0)
            thresholds = compare[compare['mode'].eq(mode) & compare.drug.eq(drug)]
            ax.axvline(-np.log10(thresholds.bonferroni_p_original.dropna().iloc[0]), color='#9A9A9A', ls=':', lw=.7)
            ax.axhline(-np.log10(thresholds.bonferroni_p_new.dropna().iloc[0]), color='#9A9A9A', ls=':', lw=.7)
            leading = frame.rank_original.le(5) | frame.rank_new.le(5)
            selected_key = frame.display_gene.isin(KEY_SYMBOLS) & (frame.rank_original.le(20) | frame.rank_new.le(20))
            labels = frame[frame.bonferroni_transition.ne('below_both') | leading | selected_key].copy()
            labels['priority'] = labels[['neglog10P_original', 'neglog10P_new']].max(axis=1)
            labels = labels.sort_values('priority', ascending=False).head(7)
            occupied = []
            for row in labels.itertuples():
                xx, yy = float(row.neglog10P_original), float(row.neglog10P_new)
                tx, ty = min(xx + .04 * lim, .80 * lim), yy + .025 * lim
                for _ in range(20):
                    if all(abs(ty - prior) > .058 * lim for prior in occupied):
                        break
                    ty += .060 * lim
                    if ty > .94 * lim:
                        ty = .12 * lim
                occupied.append(ty)
                label = row.display_gene + (' (n=1)' if row.display_gene == 'ATP8' and row.n_carriers_new == 1 else '')
                ax.annotate(label, (xx, yy), (tx, ty), fontsize=7.2,
                            arrowprops={'arrowstyle': '-', 'color': '#777777', 'lw': .45},
                            bbox={'boxstyle': 'round,pad=0.1', 'fc': 'white', 'ec': 'none', 'alpha': .8})
            ax.set(xlim=(-.01 * lim, lim), ylim=(-.01 * lim, lim), title=drug.capitalize(),
                   xlabel='Original −log₁₀P', ylabel='After exclusion −log₁₀P')
        legend_ax = axes.flat[-1]
        legend_ax.axis('off')
        handles = [Line2D([], [], marker='o', ls='none', color=palette[c], markersize=5,
                          label=t) for c, t in [('retained', 'Bonferroni significant in both'), ('gained', 'Gained Bonferroni significance'),
                                               ('lost', 'Lost Bonferroni significance'), ('below_both', 'Other tested genes')]]
        legend_ax.legend(handles=handles, frameon=False, loc='upper left', fontsize=9)
        legend_ax.text(.025, .54, 'Each point is a gene tested in both cohorts.\nDashed diagonal: equal evidence.\nDotted lines: each scan’s Bonferroni threshold.\n\nThresholds use the actual number of tested genes.\nCarrier–noncarrier effects and ranks are supplied\nin the comparison tables.', transform=legend_ax.transAxes, va='top', fontsize=9, linespacing=1.6)
        fig.suptitle('Gene GWAS after I/SDD-range exclusion\n' + MODE_LABELS[mode], fontsize=13, y=.99, fontweight='bold')
        fig.text(.5, .018, 'Quantitative endpoint: log₂MIC. Gene indicators use qualifying HIGH/MODERATE-impact variants.', ha='center', fontsize=9)
        for ext in ['png', 'pdf']:
            fig.savefig(out / f'gene_signal_comparison_{mode}.{ext}', dpi=300, bbox_inches='tight', facecolor='white')
        plt.close(fig)


def control_check(control_path, original):
    control = pd.read_csv(control_path, sep='\t')
    cols = ['drug', 'gene_id', 'p', 'b', 'se', 'neglog10P', 'n_carriers', 'n_noncarriers', 'n_missing']
    joined = control[cols].merge(original[cols], on=['drug', 'gene_id'], suffixes=('_control', '_original'), validate='one_to_one')
    if len(joined) != len(control):
        raise ValueError('Control contains genes outside the original result set')
    for field in ['n_carriers', 'n_noncarriers', 'n_missing']:
        if not joined[field + '_control'].eq(joined[field + '_original']).all():
            raise ValueError('Original-cohort control carrier counts differ: ' + field)
    return {'path': str(control_path), 'sha256': digest(control_path), 'genes_compared': len(joined),
            'drugs': sorted(control.drug.unique().tolist()),
            'max_absolute_neglog10P_difference': float((joined.neglog10P_control - joined.neglog10P_original).abs().max()),
            'max_absolute_raw_effect_difference': float((joined.b_control - joined.b_original).abs().max()),
            'max_absolute_raw_se_difference': float((joined.se_control - joined.se_original).abs().max())}


def write_report(summary, compare, top5, original, out, audit):
    lines = ['GENE-LEVEL GWAS: MIC-DEFINED I/SDD-RANGE SENSITIVITY ANALYSIS', '',
             'This analysis tests the robustness of associations with quantitative log₂MIC after removing isolates in drug-specific intermediate or susceptible-dose-dependent MIC ranges. These retrospective exclusions were inferred from recorded MICs, independently of the stored binary R/S phenotype. The exclusions alter the sampled MIC distribution; differences are interpreted as cohort sensitivity.', '',
             'Procedure',
             'The primary analysis retains the original gene indicators and removes the selected isolates (fixed). The additional analysis rebuilds indicators from HIGH/MODERATE-impact variants carried by at least one isolate and at most 5% of the remaining drug cohort (rebuilt). The latter includes changes in rare-variant eligibility as well as cohort membership. The 5% denominator is the total drug cohort. Gene presence takes precedence over missingness; absence requires all qualifying alleles to be called and absent.',
             'Both analyses use the same retained isolates and the original common-SNP relationship matrices, whose allele frequencies were defined in the full 1,153-isolate cohort: chromosome-excluded matrices for nuclear genes and the full nuclear matrix for mitochondrial genes. The quantitative endpoint and mixed-model testing follow the original scans. Effects are carrier–noncarrier contrasts in log₂MIC (2 × the raw GCTA coefficient); corresponding standard errors are also doubled.',
             'Gene ranks follow ascending P, then chromosome, midpoint and gene ID to break ties. Bonferroni and suggestive thresholds are 0.05/M and 1/M, with M recalculated from the tested genes in each scan. Top-20 overlap is the number of original top-20 genes retained in the new top 20.', '',
             'Threshold source: CLSI M27M44S, 3rd edition (2022), Table 1 and the convention assigning between-dilution MIC measurements to the next higher twofold-dilution category: https://clsi.org/media/osthhxax/m27m44sed3e_sample.pdf', '',
             'Cohorts']
    for drug in DRUGS:
        a, b = EXPECTED[drug]
        lines.append(f'{drug.capitalize()}: {a} → {b} isolates; {a-b} removed using {RANGES[drug]} µg/ml.')
    lines += ['Posaconazole retains its original 611-isolate analysis and is presented as the unchanged reference.', '', 'Findings']
    for mode in MODES:
        lines += ['', MODE_LABELS[mode]]
        for drug in DRUGS:
            s = summary[summary['mode'].eq(mode) & summary.drug.eq(drug)].iloc[0]
            lines.append(f'{drug.capitalize()}: {s.original_top20_retained}/20 original leading genes retained; tested genes {s.original_tested_genes} → {s.new_tested_genes}; Bonferroni significant {s.original_bonferroni_count} → {s.new_bonferroni_count}; suggestive-only {s.original_suggestive_only_count} → {s.new_suggestive_only_count}.')
            for change in ['retained', 'gained', 'lost']:
                names = s['bonferroni_' + change + '_genes']
                lines.append(f'  Bonferroni {change}: {names or "none"}.')
            leading = top5[top5['mode'].eq(mode) & top5.drug.eq(drug)]
            bits = [f'{r.display_gene} (rank {int(r.rank)}, P={r.p:.3g}, effect={r.carrier_vs_noncarrier_log2mic:+.3f}, carriers={int(r.n_carriers)})' for r in leading.itertuples()]
            lines.append('  New top five: ' + '; '.join(bits) + '.')
            changes = compare[compare['mode'].eq(mode) & compare.drug.eq(drug)]
            for state, label in [('entered', 'Entered top 20'), ('left', 'Left top 20')]:
                moved = changes[changes.top20_transition.eq(state)].sort_values('rank_new' if state == 'entered' else 'rank_original')
                parts = []
                for r in moved.itertuples():
                    old_rank = str(int(r.rank_original)) if pd.notna(r.rank_original) else 'untested'
                    new_rank = str(int(r.rank_new)) if pd.notna(r.rank_new) else 'untested'
                    parts.append(f'{r.display_gene} ({old_rank} → {new_rank})')
                lines.append('  ' + label + ': ' + ('; '.join(parts) or 'none') + '.')
    atp = compare[compare.drug.eq('anidulafungin') & compare.display_gene.eq('ATP8')]
    if len(atp) and atp.n_carriers_new.eq(1).all():
        lines += ['', 'ATP8 carrier count', 'The anidulafungin ATP8 association is supported by one carrier in each post-exclusion analysis. Its higher ranking and Bonferroni significance therefore describe a signal based on that single isolate, whose carrier count remains essential to interpretation.']
    lines += ['', 'Original posaconazole reference']
    pos = original[original.drug.eq('posaconazole')].sort_values('rank').head(5)
    lines.append('; '.join(f'{r.display_gene} (rank {int(r.rank)}, P={r.p:.3g}, effect={r.carrier_vs_noncarrier_log2mic:+.3f})' for r in pos.itertuples()) + '.')
    if audit.get('control'):
        c = audit['control']
        lines += ['', 'Original-cohort numerical control', f'Compared {c["genes_compared"]} gene tests across {", ".join(c["drugs"])}. Maximum absolute difference in −log₁₀P: {c["max_absolute_neglog10P_difference"]:.3g}; raw GCTA coefficient: {c["max_absolute_raw_effect_difference"]:.3g}; raw standard error: {c["max_absolute_raw_se_difference"]:.3g}. Carrier counts matched.']
    if audit.get('independent_validation'):
        v = audit['independent_validation']
        lines += ['', 'Independent validation', f'{v["status"]}: {v["scans_validated"]} scans and {v["total_gene_tests"]:,} gene tests validated, with {v["original_control_genes_exactly_reproduced"]:,} original-cohort control tests reproduced exactly. Evidence: ../validation/final_sensitivity_audit.json.']
    lines += ['', 'Reading the tables',
              'A gene marked untested lacks an eligible polymorphic indicator in that scan. Rank improvement is original rank minus new rank, so positive values indicate upward movement. The significance tables retain all original and new Bonferroni/suggestive signals, including transitions and genes that become untested. Effects retain their signs and carrier counts. Fixed and rebuilt comparisons each refer to the same original result set.', '',
              'Outputs', 'comparison_summary.tsv', 'all_gene_comparisons.tsv.gz', 'top20_comparison.tsv', 'top5_results.tsv',
              'bonferroni_changes.tsv', 'suggestive_changes.tsv', 'key_gene_comparison.tsv', 'original_posaconazole_reference.tsv',
              'gene_signal_comparison_fixed.png / .pdf', 'gene_signal_comparison_rebuilt.png / .pdf', 'source_audit.json']
    (out / 'comparison_report.txt').write_text('\n'.join(lines) + '\n')
    sections = {'Procedure', 'Cohorts', 'Findings', 'ATP8 carrier count', 'Original posaconazole reference', 'Original-cohort numerical control', 'Independent validation', 'Reading the tables', 'Outputs'}
    body = []
    for line in lines[1:]:
        if not line:
            continue
        tag = 'h2' if line in sections else ('h3' if line in MODE_LABELS.values() else 'p')
        body.append(f'<{tag}>{html.escape(line)}</{tag}>')
    images = ''.join(f'<h2>{html.escape(MODE_LABELS[m])}</h2><a href="gene_signal_comparison_{m}.pdf"><img src="gene_signal_comparison_{m}.png" alt="Original versus post-exclusion gene association evidence for all five affected drugs; points are coloured by Bonferroni-significance transition."></a>' for m in MODES)
    content = '<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Gene GWAS category sensitivity</title><style>body{max-width:1150px;margin:34px auto;padding:0 24px;font:16px/1.6 system-ui,sans-serif;color:#193344}h1{line-height:1.2}h2{margin-top:2rem;font-size:1.35rem}h3{color:#007d88}img{width:100%;height:auto}p{overflow-wrap:anywhere}a{color:#007d88}</style><h1>Gene GWAS after MIC-defined I/SDD-range exclusions</h1>' + ''.join(body) + images + '</html>'
    content = content.replace('../validation/final_sensitivity_audit.json', '<a href="../validation/final_sensitivity_audit.json">independent validation evidence</a>')
    (out / 'comparison_report.html').write_text(content)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root', type=Path, required=True, help='Original Gene_GWAS directory containing results/.')
    parser.add_argument('--analysis-root', type=Path, required=True, help='Sensitivity analysis directory containing fixed/ and rebuilt/.')
    parser.add_argument('--outdir', type=Path, help='Defaults to <analysis-root>/comparison.')
    parser.add_argument('--validate-source-only', action='store_true')
    args = parser.parse_args()
    original_path = args.source_root / 'results/all_gene_level_variant_results.tsv.gz'
    original = normalized_results(original_path, DRUGS + ['posaconazole'], 'original')
    if args.validate_source_only:
        print(original.groupby('drug').agg(cohort_n=('cohort_n', 'first'), tested_genes=('gene_id', 'size')).to_string())
        return
    sources = [original_path]
    filtered = {}
    for mode in MODES:
        path = args.analysis_root / mode / 'results/all_gene_level_variant_results.tsv.gz'
        filtered[mode] = normalized_results(path, DRUGS, mode)
        sources.append(path)
    compare, summary = comparisons(original, filtered)
    top5 = pd.concat([original[original['rank'].le(5)]] + [filtered[m][filtered[m]['rank'].le(5)] for m in MODES], ignore_index=True)
    top5 = top5[['mode', 'drug', 'gene_id', 'display_gene', 'chromosome', 'bp'] + MEASURES]
    audit = {'endpoint': 'quantitative log2 MIC', 'primary_comparison': 'fixed indicators versus original cohort',
             'additional_comparison': 'rebuilt rare-variant indicators versus original cohort',
             'rank_tie_breaking': ['p', 'chromosome', 'bp', 'gene_id'],
             'effect_definition': 'carrier-versus-noncarrier log2MIC = 2*b; SE = 2*se',
             'thresholds': 'P < 0.05/M (Bonferroni); P < 1/M (suggestive), actual tested genes per scan',
             'expected_cohort_counts': EXPECTED, 'sources': {str(p): digest(p) for p in sources}}
    control_path = args.analysis_root / 'control/results/all_gene_level_variant_results.tsv.gz'
    if control_path.exists():
        audit['control'] = control_check(control_path, original)
    validation_path = args.analysis_root / 'validation/final_sensitivity_audit.json'
    if validation_path.exists():
        evidence = json.loads(validation_path.read_text())
        if evidence['status'] != 'PASS' or evidence['total_gene_tests'] != sum(len(f) for f in filtered.values()):
            raise ValueError('Independent validation status or tested-gene count differs')
        audit['independent_validation'] = {key: evidence[key] for key in ['status', 'scans_validated', 'total_gene_tests', 'original_control_genes_exactly_reproduced']}
        audit['independent_validation'].update(path=str(validation_path), sha256=digest(validation_path))
    for path in sorted(args.analysis_root.glob('*.json')):
        audit.setdefault('analysis_provenance', {})[str(path)] = {'sha256': digest(path), 'content': json.loads(path.read_text())}
    out = args.outdir or args.analysis_root / 'comparison'
    out.mkdir(parents=True, exist_ok=True)
    table(summary, out / 'comparison_summary.tsv')
    table(compare, out / 'all_gene_comparisons.tsv.gz')
    table(compare[compare.top20_original | compare.top20_new].sort_values(['mode', 'drug', 'rank_original', 'rank_new']), out / 'top20_comparison.tsv')
    table(top5, out / 'top5_results.tsv')
    table(compare[compare.bonferroni_transition.ne('below_both')], out / 'bonferroni_changes.tsv')
    table(compare[compare.suggestive_or_stronger_transition.ne('below_both')], out / 'suggestive_changes.tsv')
    keys = compare.display_gene.isin(KEY_SYMBOLS) | compare.gene_id.str.replace('_A$', '', regex=True).isin(KEY_SYMBOLS)
    table(compare[keys], out / 'key_gene_comparison.tsv')
    table(original[original.drug.eq('posaconazole')][['drug', 'gene_id', 'display_gene', 'chromosome', 'bp'] + MEASURES], out / 'original_posaconazole_reference.tsv')
    score_figures(compare, out)
    write_report(summary, compare, top5, original, out, audit)
    (out / 'source_audit.json').write_text(json.dumps(audit, indent=2) + '\n')
    print(json.dumps({'outdir': str(out), 'comparison_rows': len(compare), 'summaries': len(summary), 'top5_rows': len(top5)}, indent=2))


if __name__ == '__main__':
    main()
