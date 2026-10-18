#!/usr/bin/env python3
"""
Extremes-only deviation flags for each personal sample (DESIGN.md Stage 5, DECISIONS.md D12).

The question is "does anything look grossly abnormal against people processed identically?",
never "is this healthy?". Every check is reported, whether raised or not, with its value and
threshold, so that "nothing stands out" is visible as a result rather than an absence.

CHECKS
------
low_alpha                  Shannon below the low percentile of the pooled panel AND of every
                           panel study (the between-study-spread rule).
single_taxon_dominance     one assigned genus above the dominance fraction (default 50%).
enterobacteriaceae_bloom   family Enterobacteriaceae above the high percentile of the pooled
                           panel AND of every study.
proteobacteria_bloom       phylum Proteobacteria (SILVA) / Pseudomonadota (GTDB), same rule.
faecalibacterium_depletion genus Faecalibacterium below the low percentile, same rule.
core_genus_missing         a genus present in >= core-prevalence of the panel has no reads,
                           although the panel median abundance at this sample's depth predicts
                           at least core-min-expected reads (so depth alone cannot explain it).
beta_outlier               Aitchison distance to the panel centroid above the high percentile of
                           the panel's own (leave-one-out) distances to that centroid.
excess_extremes            more genera in the outer percentile tails (at or below the low, at or
                           above the high percentile) than the high percentile of what panel
                           people have. Everyone has some extreme genera: with ~10% of the
                           distribution in the two tails, a typical person has about one genus in
                           ten there by chance alone. Only an excess over that is a flag.

CALIBRATION
-----------
How often would each check fire for an ordinary person? Every panel sample is treated as if
it were a personal sample and checked against the panel without it (leave-one-out), with
exactly the same code. The share of panel people for whom a check fires is reported as its
`reference_fire_rate`, next to every flag, and per check in flag_calibration.tsv. A check
that fires for 5% of ordinary people is weak evidence on its own; one that fires for none
of them is strong. The leave-one-out counts of extreme genera are also the reference for
excess_extremes.

Abundance checks use relative abundance, so that a flag reads as "31% against a panel 95th
percentile of 2.6%". Distances use CLR, with the same panel-only prevalence filter and
pseudocount as the percentiles.

A flag from one arm is provisional: DESIGN.md requires agreement between the native and the
V4-projected arm before a flag is raised in the report. The `agreement` column says so.
"""

import argparse
import sys

import numpy as np
import pandas as pd

PROTEO = ('Proteobacteria', 'Pseudomonadota')


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--asv-table', required=True)
    p.add_argument('--taxonomy', required=True, help='named-rank ASV_tax_species.<db>.tsv')
    p.add_argument('--alpha', required=True, help='alpha_diversity.tsv from ALPHA_DIVERSITY')
    p.add_argument('--samplesheet', required=True)
    p.add_argument('--prefix', default='out')
    p.add_argument('--arm', default='native')
    p.add_argument('--low-pct', type=float, default=5.0)
    p.add_argument('--high-pct', type=float, default=95.0)
    p.add_argument('--dominance', type=float, default=0.5)
    p.add_argument('--core-prevalence', type=float, default=0.95)
    p.add_argument('--core-min-expected', type=float, default=10.0)
    p.add_argument('--min-prevalence', type=float, default=0.1)
    p.add_argument('--pseudocount', type=float, default=0.5)
    p.add_argument('--min-study-n', type=int, default=10)
    return p.parse_args(argv)


def collapse(asv, tax, rank):
    """Counts per taxon at `rank`; ASVs without a call there become 'Unassigned'."""
    labels = tax.reindex(asv.index)[rank].fillna('').astype(str).str.strip()
    labels = labels.where(labels != '', 'Unassigned')
    return asv.groupby(labels).sum()


def pct(ref, v):
    ref = np.sort(np.asarray(ref, dtype=float))
    return 100.0 * np.searchsorted(ref, v, side='right') / len(ref)


class Panel:
    """Pooled and per-study reference distributions for one quantity."""

    def __init__(self, values, study, min_n):
        self.values = values
        self.groups = {s: values[study == s] for s in np.unique(study)
                       if (study == s).sum() >= min_n}

    def below_all(self, v, q):
        pooled = np.percentile(self.values, q)
        return v < pooled and all(v < np.percentile(x, q) for x in self.groups.values()), pooled

    def above_all(self, v, q):
        pooled = np.percentile(self.values, q)
        return v > pooled and all(v > np.percentile(x, q) for x in self.groups.values()), pooled


def main(argv=None):
    a = parse_args(argv)
    asv = pd.read_csv(a.asv_table, sep='\t', index_col=0)
    tax = pd.read_csv(a.taxonomy, sep='\t', index_col=0)
    sheet = pd.read_csv(a.samplesheet, dtype=str).fillna('').set_index('sample')
    if 'context' not in sheet.columns:
        sheet['context'] = ''
    alpha = pd.read_csv(a.alpha, sep='\t').set_index('sample')
    for rank in ('Phylum', 'Family', 'Genus'):
        if rank not in tax.columns:
            sys.exit("FATAL: taxonomy has no '%s' column; a named-rank table is required" % rank)
    missing = sorted(set(sheet.index) - set(asv.columns))
    if missing:
        sys.exit('FATAL: samplesheet samples missing from the ASV table: %s' % ', '.join(missing[:10]))

    asv = asv[list(sheet.index)]
    panel_ids = list(sheet.index[sheet['role'] == 'panel'])
    personal_ids = list(sheet.index[sheet['role'] == 'personal'])

    ranks = {r: collapse(asv, tax, r) for r in ('Phylum', 'Family', 'Genus')}
    rel = {r: c / c.sum(axis=0).replace(0, np.nan) for r, c in ranks.items()}
    genus = ranks['Genus']
    depth = genus.sum(axis=0)
    present = genus > 0

    def rank_share(rank, names):
        tbl = rel[rank]
        hit = [n for n in names if n in tbl.index]
        return tbl.loc[hit].sum(axis=0) if hit else pd.Series(0.0, index=tbl.columns)

    entero = rank_share('Family', ['Enterobacteriaceae'])
    proteo = rank_share('Phylum', PROTEO)
    faecali = rank_share('Genus', ['Faecalibacterium'])
    shannon = alpha['shannon'].reindex(sheet.index)

    def evaluate(s, ref):
        """Every check for sample `s` against the reference samples `ref`.

        Used for the personal samples against the whole panel and, for calibration, for
        each panel sample against the panel without it. Returns {flag: (raised, value,
        threshold, detail)}, the number of extreme genera and the number of genera tested.
        """
        study = sheet.loc[ref, 'study'].to_numpy()
        res = {}
        # low alpha
        if pd.isna(shannon[s]):
            res['low_alpha'] = ('not_evaluable', '', '', 'below the rarefaction depth')
        else:
            ok = shannon[ref].notna().to_numpy()
            p_sh = Panel(shannon[ref].to_numpy(dtype=float)[ok], study[ok], a.min_study_n)
            hit, thr = p_sh.below_all(shannon[s], a.low_pct)
            res['low_alpha'] = ('yes' if hit else 'no', round(float(shannon[s]), 3),
                                'p%g = %.3f' % (a.low_pct, thr),
                                'pooled-panel percentile %.1f' % pct(p_sh.values, shannon[s]))
        # dominance
        g = rel['Genus'][s].drop('Unassigned', errors='ignore')
        top, frac = g.idxmax(), float(g.max())
        res['single_taxon_dominance'] = ('yes' if frac > a.dominance else 'no', round(frac, 4),
                                         '> %g' % a.dominance, 'top genus %s' % top)
        # blooms and depletion
        for flag, series, high in (('enterobacteriaceae_bloom', entero, True),
                                   ('proteobacteria_bloom', proteo, True),
                                   ('faecalibacterium_depletion', faecali, False)):
            ref_panel = Panel(series[ref].to_numpy(dtype=float), study, a.min_study_n)
            v = float(series[s])
            hit, thr = (ref_panel.above_all(v, a.high_pct) if high else ref_panel.below_all(v, a.low_pct))
            res[flag] = ('yes' if hit else 'no', round(v, 4),
                         'p%g = %.4f' % (a.high_pct if high else a.low_pct, thr),
                         'pooled-panel percentile %.1f' % pct(ref_panel.values, v))
        # core genera
        prevalence = present[ref].mean(axis=1)
        core = prevalence[prevalence >= a.core_prevalence].index.difference(['Unassigned'])
        median_rel = rel['Genus'][ref].median(axis=1)
        gone = [c for c in core if genus.at[c, s] == 0 and median_rel[c] * depth[s] >= a.core_min_expected]
        res['core_genus_missing'] = ('yes' if gone else 'no', len(gone),
                                     'prevalence >= %g, expected >= %g reads'
                                     % (a.core_prevalence, a.core_min_expected),
                                     ', '.join(gone) if gone else '%d core genera all present' % len(core))
        # CLR on genus, reference-only prevalence filter, as in TAXON_PERCENTILES
        kept = prevalence[prevalence >= a.min_prevalence].index
        logx = np.log(genus.loc[kept, ref + [s]].astype(float) + a.pseudocount)
        clr = logx - logx.mean(axis=0)
        # A sample's distance to a centroid it is part of is shorter than to one it is not:
        # leaving it out scales the distance by exactly n / (n - 1). The reference distances
        # are scaled so, to be comparable with a sample that is never part of the centroid.
        centroid = clr[ref].mean(axis=1)
        ref_dist = np.sqrt(((clr[ref].sub(centroid, axis=0)) ** 2).sum(axis=0)).to_numpy()
        ref_dist = ref_dist * len(ref) / (len(ref) - 1)
        v = float(np.sqrt(((clr[s] - centroid) ** 2).sum()))
        thr = float(np.percentile(ref_dist, a.high_pct))
        res['beta_outlier'] = ('yes' if v > thr else 'no', round(v, 3), 'p%g = %.3f' % (a.high_pct, thr),
                               'distance-to-centroid percentile %.1f' % pct(ref_dist, v))
        # extreme genera, with the percentile and category conventions of TAXON_PERCENTILES
        # and PERCENTILE_STABILITY (share of the reference at or below; <= low, >= high)
        rows_g = [k for k in kept if k != 'Unassigned']
        R = clr.loc[rows_g, ref].to_numpy()
        p = 100.0 * (R <= clr.loc[rows_g, s].to_numpy()[:, None]).sum(axis=1) / R.shape[1]
        n_extreme = int(((p <= a.low_pct) | (p >= a.high_pct)).sum())
        return res, n_extreme, len(rows_g)

    # calibration: each panel person against the panel without them
    loo = []
    for j in panel_ids:
        res, n_ext, _ = evaluate(j, [i for i in panel_ids if i != j])
        loo.append((res, n_ext))
    loo_ext = np.array([n for _, n in loo], dtype=float)
    q_ext = float(np.percentile(loo_ext, a.high_pct))
    for k, (res, n_ext) in enumerate(loo):
        others = np.delete(loo_ext, k)
        res['excess_extremes'] = ('yes' if n_ext > np.percentile(others, a.high_pct) else 'no',)

    flags = ['low_alpha', 'single_taxon_dominance', 'enterobacteriaceae_bloom', 'proteobacteria_bloom',
             'faecalibacterium_depletion', 'core_genus_missing', 'beta_outlier', 'excess_extremes']
    nominal = {'low_alpha': '<= %g%%' % a.low_pct,
               'enterobacteriaceae_bloom': '<= %g%%' % (100 - a.high_pct),
               'proteobacteria_bloom': '<= %g%%' % (100 - a.high_pct),
               'faecalibacterium_depletion': '<= %g%%' % a.low_pct,
               'beta_outlier': '%g%%' % (100 - a.high_pct),
               'excess_extremes': '%g%%' % (100 - a.high_pct)}
    calib, rate = [], {}
    for f in flags:
        states = [res[f][0] for res, _ in loo]
        n = sum(st in ('yes', 'no') for st in states)
        fired = sum(st == 'yes' for st in states)
        rate[f] = round(fired / n, 4) if n else ''
        calib.append({'flag': f, 'reference_n': n, 'reference_fired': fired, 'reference_fire_rate': rate[f],
                      'nominal_rate': nominal.get(f, ''), 'note': ''})

    rows, n_genera = [], 0

    def add(sample, flag, raised, value, threshold, detail):
        rows.append({'sample': sample, 'context': sheet.at[sample, 'context'], 'flag': flag,
                     'raised': raised, 'value': value, 'threshold': threshold, 'detail': detail,
                     'reference_fire_rate': rate[flag], 'arm': a.arm,
                     'agreement': 'single arm: needs the other arm to confirm'})

    for s in personal_ids:
        res, n_ext, n_genera = evaluate(s, panel_ids)
        for f in flags[:-1]:
            add(s, f, *res[f])
        add(s, 'excess_extremes', 'yes' if n_ext > q_ext else 'no', n_ext, '> p%g = %g' % (a.high_pct, q_ext),
            'of %d genera; a reference person has a median of %g in the tails' % (n_genera, np.median(loo_ext)))

    chance = n_genera * (a.low_pct + 100 - a.high_pct) / 100.0
    notes = {
        'single_taxon_dominance': 'absolute threshold, no nominal rate',
        'core_genus_missing': 'depends on how many core genera the panel has',
        'beta_outlier': "percentile of the panel's own distances to its centroid",
        'excess_extremes': ('reference people have a median of %g extreme genera of %d (p%g = %g); '
                            'chance alone, if genera were independent, would give about %.0f'
                            % (np.median(loo_ext), n_genera, a.high_pct, q_ext, chance)),
    }
    for c in calib:
        c['note'] = notes.get(c['flag'], 'pooled AND every study, so at most the nominal rate')
    pd.DataFrame(calib).to_csv('%s.flag_calibration.tsv' % a.prefix, sep='\t', index=False)

    out = pd.DataFrame(rows)
    out.to_csv('%s.deviation_flags.tsv' % a.prefix, sep='\t', index=False)
    raised = out[out['raised'] == 'yes']
    print('%d checks on %d personal samples, %d raised (single-arm, provisional)'
          % (len(out), len(personal_ids), len(raised)))
    for _, r in raised.iterrows():
        print('  %s %s: %s (%s; fires for %s of reference people)'
              % (r['sample'], r['flag'], r['value'], r['detail'], r['reference_fire_rate']))


if __name__ == '__main__':
    main()
