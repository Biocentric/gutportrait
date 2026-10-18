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
                           the panel's own distances to that centroid.

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
    study = sheet.loc[panel_ids, 'study'].to_numpy()

    ranks = {r: collapse(asv, tax, r) for r in ('Phylum', 'Family', 'Genus')}
    rel = {r: c / c.sum(axis=0).replace(0, np.nan) for r, c in ranks.items()}
    genus = ranks['Genus']
    depth = genus.sum(axis=0)

    # CLR on genus, panel-only prevalence filter, as in TAXON_PERCENTILES
    prevalence = (genus[panel_ids] > 0).mean(axis=1)
    kept = prevalence[prevalence >= a.min_prevalence].index
    logx = np.log(genus.loc[kept].astype(float) + a.pseudocount)
    clr = logx - logx.mean(axis=0)
    centroid = clr[panel_ids].mean(axis=1)
    dist = np.sqrt(((clr.sub(centroid, axis=0)) ** 2).sum(axis=0))
    panel_dist = dist[panel_ids].to_numpy()

    def panel_of(series):
        return Panel(series[panel_ids].to_numpy(dtype=float), study, a.min_study_n)

    def rank_share(rank, names):
        tbl = rel[rank]
        hit = [n for n in names if n in tbl.index]
        return tbl.loc[hit].sum(axis=0) if hit else pd.Series(0.0, index=tbl.columns)

    entero = rank_share('Family', ['Enterobacteriaceae'])
    proteo = rank_share('Phylum', PROTEO)
    faecali = rank_share('Genus', ['Faecalibacterium'])
    shannon = alpha['shannon'].reindex(sheet.index)
    core = prevalence[prevalence >= a.core_prevalence].index.difference(['Unassigned'])
    median_rel = rel['Genus'][panel_ids].median(axis=1)

    rows = []

    def add(sample, flag, raised, value, threshold, detail):
        rows.append({'sample': sample, 'context': sheet.at[sample, 'context'], 'flag': flag,
                     'raised': raised, 'value': value, 'threshold': threshold, 'detail': detail,
                     'arm': a.arm, 'agreement': 'single arm: needs the other arm to confirm'})

    p_entero, p_proteo, p_faecali = panel_of(entero), panel_of(proteo), panel_of(faecali)
    shannon_ok = shannon[panel_ids].notna()
    p_shannon = Panel(shannon[panel_ids][shannon_ok].to_numpy(),
                      study[shannon_ok.to_numpy()], a.min_study_n)

    for s in personal_ids:
        # low alpha
        if pd.isna(shannon[s]):
            add(s, 'low_alpha', 'not_evaluable', '', '', 'below the rarefaction depth')
        else:
            hit, thr = p_shannon.below_all(shannon[s], a.low_pct)
            add(s, 'low_alpha', 'yes' if hit else 'no', round(float(shannon[s]), 3),
                'p%g = %.3f' % (a.low_pct, thr), 'pooled-panel percentile %.1f'
                % pct(p_shannon.values, shannon[s]))
        # dominance
        g = rel['Genus'][s].drop('Unassigned', errors='ignore')
        top, frac = g.idxmax(), float(g.max())
        add(s, 'single_taxon_dominance', 'yes' if frac > a.dominance else 'no', round(frac, 4),
            '> %g' % a.dominance, 'top genus %s' % top)
        # blooms and depletion
        for flag, series, ref, high in (('enterobacteriaceae_bloom', entero, p_entero, True),
                                        ('proteobacteria_bloom', proteo, p_proteo, True),
                                        ('faecalibacterium_depletion', faecali, p_faecali, False)):
            v = float(series[s])
            hit, thr = (ref.above_all(v, a.high_pct) if high else ref.below_all(v, a.low_pct))
            add(s, flag, 'yes' if hit else 'no', round(v, 4),
                'p%g = %.4f' % (a.high_pct if high else a.low_pct, thr),
                'pooled-panel percentile %.1f' % pct(ref.values, v))
        # core genera
        gone = [c for c in core if genus.at[c, s] == 0 and median_rel[c] * depth[s] >= a.core_min_expected]
        add(s, 'core_genus_missing', 'yes' if gone else 'no', len(gone),
            'prevalence >= %g, expected >= %g reads' % (a.core_prevalence, a.core_min_expected),
            ', '.join(gone) if gone else '%d core genera all present' % len(core))
        # beta outlier
        thr = float(np.percentile(panel_dist, a.high_pct))
        v = float(dist[s])
        add(s, 'beta_outlier', 'yes' if v > thr else 'no', round(v, 3), 'p%g = %.3f' % (a.high_pct, thr),
            'distance-to-centroid percentile %.1f' % pct(panel_dist, v))

    out = pd.DataFrame(rows)
    out.to_csv('%s.deviation_flags.tsv' % a.prefix, sep='\t', index=False)
    raised = out[out['raised'] == 'yes']
    print('%d checks on %d personal samples, %d raised (single-arm, provisional)'
          % (len(out), len(personal_ids), len(raised)))
    for _, r in raised.iterrows():
        print('  %s %s: %s (%s)' % (r['sample'], r['flag'], r['value'], r['detail']))


if __name__ == '__main__':
    main()
