#!/usr/bin/env python3
"""
Position personal samples against a reference panel, per genus, in CLR space.

WHY CLR AND NOT RELATIVE ABUNDANCE
----------------------------------
A relative abundance is a ratio to an unmeasured total. Two samples can differ in one
genus's percentage purely because a *different* genus changed, with no biological change
in the genus being reported. Percentiles computed on raw proportions therefore manufacture
findings out of compositional arithmetic -- this is the mechanism behind "your Akkermansia
is low" style claims. The centred-log-ratio removes the shared scale by dividing each part
by the sample's own geometric mean, so a value only moves when that part moves relative to
the rest of the composition.

WHAT THIS DOES NOT FIX
----------------------
CLR removes the *compositional* artefact. It does nothing about the wet-lab confound: if
the panel was amplified with different primers than the personal samples, the two are still
measuring different things and the percentile is orientation only. That is what Stage 3's
batch-dominance gate is for. This module reports positions; it does not certify them.

ZERO HANDLING
-------------
CLR is undefined at zero. We use a simple additive pseudocount (default 0.5) applied to the
aligned count matrix before the transform, which is the transparent option. Bayesian
multiplicative replacement is better-behaved for very sparse data; it is deliberately not
used here because its behaviour is harder to explain in a report that has to be honest
about its own uncertainty. The pseudocount used is written into the output.
"""

import argparse
import sys

import numpy as np
import pandas as pd


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--personal', required=True, help='personal genus counts TSV (genus x sample)')
    p.add_argument('--panel', required=True, help='panel genus counts TSV (genus x sample)')
    p.add_argument('--prefix', default='out')
    p.add_argument('--pseudocount', type=float, default=0.5)
    p.add_argument('--min-prevalence', type=float, default=0.10,
                   help='drop genera present in fewer than this fraction of PANEL samples '
                        '(default 0.10). Rare taxa have replicate ICC <= 0.48 and are noise '
                        'at n-of-1.')
    p.add_argument('--drop-unassigned', action='store_true',
                   help='exclude the Unassigned row from the composition. Off by default: '
                        'it is part of the composition and dropping it silently inflates '
                        'everything else.')
    return p.parse_args(argv)


def clr(counts, pseudocount):
    """CLR by sample. counts: genera (rows) x samples (cols)."""
    x = counts.astype(float) + pseudocount
    logx = np.log(x)
    return logx - logx.mean(axis=0)          # subtract each SAMPLE's geometric mean


def main(argv=None):
    a = parse_args(argv)

    personal = pd.read_csv(a.personal, sep='\t', index_col=0)
    panel = pd.read_csv(a.panel, sep='\t', index_col=0)

    if a.drop_unassigned:
        personal = personal.drop(index='Unassigned', errors='ignore')
        panel = panel.drop(index='Unassigned', errors='ignore')

    # Prevalence filter on the PANEL only -- filtering on the personal samples would
    # select features for being present in you, which is circular.
    prevalence = (panel > 0).sum(axis=1) / panel.shape[1]
    keep = prevalence[prevalence >= a.min_prevalence].index
    if len(keep) == 0:
        sys.exit('FATAL: no genera pass the %.0f%% panel prevalence filter'
                 % (100 * a.min_prevalence))

    # Align to a shared feature set. Genera absent from one side become zeros, so both
    # sides' geometric means are taken over the SAME parts -- otherwise the CLR values
    # are not on a common scale and the percentiles are meaningless.
    shared = keep.union(personal.index).intersection(keep)
    panel_a = panel.reindex(shared).fillna(0)
    personal_a = personal.reindex(shared).fillna(0)

    panel_clr = clr(panel_a, a.pseudocount)
    personal_clr = clr(personal_a, a.pseudocount)

    rows = []
    for genus in shared:
        ref = panel_clr.loc[genus].to_numpy()
        ref_sorted = np.sort(ref)
        for sample in personal_clr.columns:
            val = personal_clr.at[genus, sample]
            pct = 100.0 * np.searchsorted(ref_sorted, val, side='right') / len(ref_sorted)
            rows.append({
                'genus': genus,
                'sample': sample,
                'clr': round(float(val), 4),
                'percentile': round(float(pct), 1),
                'panel_median_clr': round(float(np.median(ref)), 4),
                'panel_n': len(ref),
                'panel_prevalence': round(float(prevalence[genus]), 3),
                'personal_counts': int(personal_a.at[genus, sample]),
            })

    out = pd.DataFrame(rows).sort_values(['sample', 'percentile'],
                                         ascending=[True, False])
    out.to_csv('%s.taxon_percentiles.tsv' % a.prefix, sep='\t', index=False)

    meta = pd.DataFrame([{
        'pseudocount': a.pseudocount,
        'min_prevalence': a.min_prevalence,
        'drop_unassigned': a.drop_unassigned,
        'n_genera_tested': len(shared),
        'n_genera_panel_total': panel.shape[0],
        'panel_samples': panel.shape[1],
        'personal_samples': personal.shape[1],
    }])
    meta.to_csv('%s.percentile_params.tsv' % a.prefix, sep='\t', index=False)

    print('tested %d genera (of %d in panel) against %d panel samples'
          % (len(shared), panel.shape[0], panel.shape[1]))
    print('wrote %s.taxon_percentiles.tsv' % a.prefix)


if __name__ == '__main__':
    main()
