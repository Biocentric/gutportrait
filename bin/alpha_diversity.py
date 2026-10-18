#!/usr/bin/env python3
"""
Alpha diversity at a fixed sequencing depth, and where each personal sample sits.

WHY RAREFY HERE
---------------
Observed richness grows with sequencing depth, and the samples of a co-processed run
differ several-fold in depth. Comparing them without a common depth would rank people
by how deeply they were sequenced. Every sample is therefore subsampled to the same
depth, repeatedly with fixed seeds, and the metrics are averaged over the iterations.
Samples below the depth are reported, not silently dropped.

Shannon is computed on ASVs with the natural log; Pielou's evenness is Shannon divided by
the log of observed richness.

WHERE A SAMPLE SITS
-------------------
For every personal sample and metric: its percentile in the pooled panel (share of panel
samples at or below it), the panel median, the range of the study medians, and whether it
lies below the 5th (or above the 95th) percentile of EVERY panel study with enough
samples. That last test is the between-study-spread rule of DECISIONS.md D12: a value is
only unusual if no single reference study would call it ordinary.
"""

import argparse
import sys

import numpy as np
import pandas as pd

METRICS = ('shannon', 'observed', 'pielou')


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--asv-table', required=True, help='ampliseq dada2/ASV_table.tsv')
    p.add_argument('--samplesheet', required=True, help='samplesheet CSV (sample, role, study, ...)')
    p.add_argument('--prefix', default='out')
    p.add_argument('--depth', type=int, default=5000, help='rarefaction depth (reads)')
    p.add_argument('--iterations', type=int, default=20, help='rarefaction repeats to average')
    p.add_argument('--seed', type=int, default=20261018)
    p.add_argument('--min-study-n', type=int, default=10,
                   help='panel studies smaller than this are left out of the per-study bounds')
    p.add_argument('--low-pct', type=float, default=5.0)
    p.add_argument('--high-pct', type=float, default=95.0)
    return p.parse_args(argv)


def read_inputs(a):
    asv = pd.read_csv(a.asv_table, sep='\t', index_col=0)
    sheet = pd.read_csv(a.samplesheet, dtype=str).fillna('')
    for col in ('sample', 'role', 'study'):
        if col not in sheet.columns:
            sys.exit("FATAL: samplesheet has no '%s' column" % col)
    missing = sorted(set(sheet['sample']) - set(asv.columns))
    if missing:
        sys.exit('FATAL: %d samplesheet sample(s) have no column in the ASV table: %s'
                 % (len(missing), ', '.join(missing[:10])))
    if 'context' not in sheet.columns:
        sheet['context'] = ''
    return asv[list(sheet['sample'])], sheet.set_index('sample')


def rarefied_metrics(counts, depth, iterations, seed):
    """Mean Shannon / observed / Pielou over `iterations` rarefactions of one sample."""
    acc = np.zeros(3)
    for i in range(iterations):
        rng = np.random.default_rng(seed + i)
        sub = rng.multivariate_hypergeometric(counts, depth)
        nz = sub[sub > 0]
        p = nz / depth
        shannon = float(-(p * np.log(p)).sum())
        observed = float(len(nz))
        pielou = shannon / np.log(observed) if observed > 1 else 0.0
        acc += (shannon, observed, pielou)
    return acc / iterations


def percentile(ref, value):
    ref = np.sort(np.asarray(ref, dtype=float))
    return 100.0 * np.searchsorted(ref, value, side='right') / len(ref)


def main(argv=None):
    a = parse_args(argv)
    asv, sheet = read_inputs(a)

    rows = []
    for sample in asv.columns:
        counts = asv[sample].to_numpy(dtype=np.int64)
        total = int(counts.sum())
        row = {'sample': sample, 'role': sheet.at[sample, 'role'], 'study': sheet.at[sample, 'study'],
               'context': sheet.at[sample, 'context'], 'reads': total}
        if total < a.depth:
            row.update({'rarefied': False, 'shannon': np.nan, 'observed': np.nan, 'pielou': np.nan})
        else:
            sh, ob, pi = rarefied_metrics(counts, a.depth, a.iterations, a.seed)
            row.update({'rarefied': True, 'shannon': round(sh, 4), 'observed': round(ob, 2),
                        'pielou': round(pi, 4)})
        rows.append(row)
    alpha = pd.DataFrame(rows)
    alpha.to_csv('%s.alpha_diversity.tsv' % a.prefix, sep='\t', index=False)

    low = alpha[~alpha['rarefied']]
    for _, r in low.iterrows():
        print('note: %s (%s) has %d reads, below the rarefaction depth %d; no alpha diversity'
              % (r['sample'], r['role'], r['reads'], a.depth), file=sys.stderr)

    panel = alpha[(alpha['role'] == 'panel') & alpha['rarefied']]
    personal = alpha[(alpha['role'] == 'personal') & alpha['rarefied']]
    if panel.empty:
        sys.exit('FATAL: no panel sample reaches the rarefaction depth %d' % a.depth)
    study_sizes = panel.groupby('study').size()
    studies = [s for s, n in study_sizes.items() if n >= a.min_study_n]

    pos = []
    for _, r in personal.iterrows():
        for m in METRICS:
            v = float(r[m])
            ref = panel[m].to_numpy()
            per_study = {s: panel.loc[panel['study'] == s, m].to_numpy() for s in studies}
            medians = [float(np.median(x)) for x in per_study.values()]
            below_all = bool(studies) and all(v < np.percentile(x, a.low_pct) for x in per_study.values())
            above_all = bool(studies) and all(v > np.percentile(x, a.high_pct) for x in per_study.values())
            pos.append({
                'sample': r['sample'], 'context': r['context'], 'metric': m, 'value': round(v, 4),
                'panel_percentile': round(percentile(ref, v), 1), 'panel_n': len(ref),
                'panel_median': round(float(np.median(ref)), 4),
                'study_median_min': round(min(medians), 4) if medians else np.nan,
                'study_median_max': round(max(medians), 4) if medians else np.nan,
                'n_studies_in_bounds': len(studies),
                'below_every_study_p%g' % a.low_pct: below_all,
                'above_every_study_p%g' % a.high_pct: above_all,
            })
    pd.DataFrame(pos).to_csv('%s.alpha_position.tsv' % a.prefix, sep='\t', index=False)

    pd.DataFrame([{
        'depth': a.depth, 'iterations': a.iterations, 'seed': a.seed, 'min_study_n': a.min_study_n,
        'samples_below_depth': len(low), 'panel_samples_used': len(panel),
        'personal_samples_used': len(personal), 'studies_in_bounds': len(studies),
    }]).to_csv('%s.alpha_params.tsv' % a.prefix, sep='\t', index=False)
    print('alpha diversity at %d reads x %d iterations: %d panel, %d personal, %d below depth'
          % (a.depth, a.iterations, len(panel), len(personal), len(low)))


if __name__ == '__main__':
    main()
