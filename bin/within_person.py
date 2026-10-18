#!/usr/bin/env python3
"""
How much do YOUR samples differ from each other, against how much one person's samples
normally differ? And are they from one person at all?

THE BORROWED NOISE FLOOR
------------------------
A gut microbiome is never sampled twice the same: day-to-day variation, stool
heterogeneity, extraction and sequencing all move a profile. Whether a change between
two of your samples means anything depends on how large that ordinary variation is,
and three samples of one person cannot measure it. So it is borrowed from the panel:
some panel studies sampled the same people twice (samplesheet role `repeat`, with
`subject` naming the person). Their within-person differences are the reference for
"ordinary change"; differences between DIFFERENT people of the same study are the
reference for "this looks like someone else".

Both references come from within one study, so they compare samples processed in the
same lab. Your samples are compared with each other only, within your own batch.

DISTANCES
---------
Aitchison distance: Euclidean distance between CLR genus profiles, the whole-profile
difference used throughout the pipeline. It adds up every change.
Trimmed distance: the same, after removing any shift shared by all genera (the CLR
differences are centred on their median) and leaving out the 5% of genera that changed
most. One genus that blooms or crashes moves the Aitchison distance a lot and the trimmed
distance hardly at all; a different person differs in most genera and moves both.

OUTPUTS
-------
within_person_reference.tsv   one row per reference pair (two samples of one panel person)
within_person_summary.tsv     the within- and between-person distributions of both distances
personal_changes.tsv          one row per pair of your samples: the distance, its percentile
                              among within- and between-person pairs, the change in Shannon
                              diversity, and how many genera changed more than they usually do
personal_genus_changes.tsv    per pair and genus: the change in CLR abundance in units of the
                              ordinary within-person change of that genus (z)
identity.tsv                  per pair: same person, same person with a large change in a few
                              genera, or possibly different people; and one overall verdict

GENUS CHANGES
-------------
For each genus, the ordinary change is the root-mean-square change in CLR abundance over
the reference pairs in which the genus was seen in at least one of the two samples (a
genus absent from both samples did not change, and counting those pairs would make every
appearance look extraordinary). A change is reported when |z| > --change-z. Some genera
pass that by chance in anyone, so the same count is made for every reference pair, each
measured against the other pairs (leave-one-out): `genera_changed_expected` is what a
typical person shows between two samples.

IDENTITY RULE
-------------
You sent in your own samples, so "one person" is the starting assumption, and the check
looks for evidence against it, as a guard against a mix-up at the lab. A distance is
ORDINARY for one person when it is at or below the --identity-within-pct percentile of the
within-person pairs, or below the --identity-between-pct percentile of the between-person
pairs (closer than nearly all pairs of different people).
same_person                 the Aitchison distance is ordinary
same_person_large_change    it is not, but the trimmed distance is: the two samples differ
                            in a few genera only, as after a bloom or a crash of one group,
                            which happens within one person
possibly_different_people   neither is ordinary: most genera differ, as between two people
Without enough reference pairs only the between-person half of "ordinary" can be applied,
nothing is called a different person (`unclear` instead), and the detail says so.

On real panels the Aitchison distance separates the two references best; the trimmed
distance is the second step only, because trimming also removes real differences.

CLR uses the panel-only prevalence filter and pseudocount used everywhere else in the
pipeline.
"""

import argparse
import itertools
import sys

import numpy as np
import pandas as pd


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--counts', required=True, help='genus counts (uncorrected) for all samples')
    p.add_argument('--alpha', required=True, help='alpha_diversity.tsv from ALPHA_DIVERSITY')
    p.add_argument('--samplesheet', required=True)
    p.add_argument('--prefix', default='out')
    p.add_argument('--min-prevalence', type=float, default=0.1)
    p.add_argument('--pseudocount', type=float, default=0.5)
    p.add_argument('--change-z', type=float, default=2.0)
    p.add_argument('--min-pairs', type=int, default=10,
                   help='fewer reference pairs than this: no within-person reference')
    p.add_argument('--identity-within-pct', type=float, default=95.0)
    p.add_argument('--identity-between-pct', type=float, default=5.0)
    p.add_argument('--trim', type=float, default=0.05,
                   help='share of the most-changed genera left out of the trimmed distance')
    return p.parse_args(argv)


def pct(ref, v):
    """Share of the reference at or below v, in percent (the pipeline-wide convention)."""
    ref = np.sort(np.asarray(ref, dtype=float))
    return 100.0 * np.searchsorted(ref, v, side='right') / len(ref)


def rnd(v, k=3):
    return '' if v is None or not np.isfinite(v) else round(float(v), k)


def read_sheet(path):
    sheet = pd.read_csv(path, dtype=str).fillna('')
    for col in ('sample', 'role', 'study', 'run'):
        if col not in sheet.columns:
            sys.exit("FATAL: samplesheet has no '%s' column" % col)
    for col in ('context', 'subject'):
        if col not in sheet.columns:
            sheet[col] = ''
    return sheet.set_index('sample')


def reference_people(sheet):
    """{(study, subject): [samples]} for every panel person with at least one repeat."""
    people = {}
    for sid, r in sheet[sheet['role'] == 'panel'].iterrows():
        people[(r['study'], r['subject'] or sid)] = [sid]
    for sid, r in sheet[sheet['role'] == 'repeat'].iterrows():
        key = (r['study'], r['subject'])
        if key not in people:
            sys.exit("FATAL: repeat %s names subject '%s', which is no panel person of %s"
                     % (sid, r['subject'], r['study']))
        people[key].append(sid)
    return {k: v for k, v in people.items() if len(v) > 1}


class Reference:
    """Within- and between-person distributions of one distance, and the 'ordinary' test."""

    def __init__(self, within, between, a):
        self.within = np.asarray(within, dtype=float)
        self.between = np.asarray(between, dtype=float)
        self.have_within = len(self.within) >= a.min_pairs
        self.within_thr = float(np.percentile(self.within, a.identity_within_pct)) if self.have_within else np.nan
        self.between_thr = float(np.percentile(self.between, a.identity_between_pct))

    def within_pct(self, v):
        return pct(self.within, v) if self.have_within else np.nan

    def between_pct(self, v):
        return pct(self.between, v)

    def ordinary(self, v):
        return v < self.between_thr or (self.have_within and v <= self.within_thr)


def main(argv=None):
    a = parse_args(argv)
    sheet = read_sheet(a.samplesheet)
    counts = pd.read_csv(a.counts, sep='\t', index_col=0)
    missing = sorted(set(sheet.index) - set(counts.columns))
    if missing:
        sys.exit('FATAL: samplesheet samples missing from the counts table: %s' % ', '.join(missing[:10]))
    counts = counts[list(sheet.index)]
    shannon = pd.read_csv(a.alpha, sep='\t').set_index('sample')['shannon'].reindex(sheet.index)

    panel_ids = list(sheet.index[sheet['role'] == 'panel'])
    personal_ids = list(sheet.index[sheet['role'] == 'personal'])
    prev = (counts[panel_ids] > 0).mean(axis=1)
    kept = prev[prev >= a.min_prevalence].index
    logx = np.log(counts.loc[kept].astype(float) + a.pseudocount)
    clr = logx - logx.mean(axis=0)
    genera = [g for g in kept if g != 'Unassigned']
    C = clr.loc[genera]
    P = counts.loc[genera] > 0
    n_trim = max(1, int(round(a.trim * len(kept))))

    def dist(x, y):
        return float(np.linalg.norm(clr[x] - clr[y]))

    def trimmed(x, y):
        d = (clr[y] - clr[x]).to_numpy()
        d2 = np.sort((d - np.median(d)) ** 2)
        return float(np.sqrt(d2[:-n_trim].sum()))

    # --- the two references ---------------------------------------------------------
    people = reference_people(sheet)
    ref = []
    for (study, subject), members in sorted(people.items()):
        for x, y in itertools.combinations(members, 2):
            ref.append({'study': study, 'subject': subject, 'sample_a': x, 'sample_b': y,
                        'same_run': sheet.at[x, 'run'] == sheet.at[y, 'run'],
                        'aitchison': dist(x, y), 'trimmed_distance': trimmed(x, y),
                        'delta_shannon': float(shannon[y] - shannon[x])})
    ref = pd.DataFrame(ref, columns=['study', 'subject', 'sample_a', 'sample_b', 'same_run', 'aitchison',
                                     'trimmed_distance', 'delta_shannon'])
    have_ref = len(ref) >= a.min_pairs

    btw = []
    for study, grp in sheet.loc[panel_ids].groupby('study'):
        for x, y in itertools.combinations(list(grp.index), 2):
            btw.append((study, dist(x, y), trimmed(x, y)))
    if not btw:
        sys.exit('FATAL: no two panel people share a study; there is no between-person reference')
    btw = pd.DataFrame(btw, columns=['study', 'aitchison', 'trimmed_distance'])
    full = Reference(ref['aitchison'], btw['aitchison'], a)
    trim = Reference(ref['trimmed_distance'], btw['trimmed_distance'], a)

    # --- per-genus ordinary change, and the leave-one-out count for each reference pair --
    n_tot = np.zeros(len(genera), dtype=int)
    sd = np.full(len(genera), np.nan)
    ref['genera_changed_loo'] = 0
    if len(ref):
        A, B = list(ref['sample_a']), list(ref['sample_b'])
        dclr = C[B].to_numpy() - C[A].to_numpy()                       # genus x pair
        seen = P[A].to_numpy() | P[B].to_numpy()
        sq = np.where(seen, dclr ** 2, 0.0)
        s_tot, n_tot = sq.sum(axis=1), seen.sum(axis=1)
        ok = n_tot >= a.min_pairs
        sd[ok] = np.sqrt(s_tot[ok] / n_tot[ok])
        n_loo = n_tot[:, None] - seen
        sd_loo = np.full(sq.shape, np.nan)
        m = n_loo >= a.min_pairs
        sd_loo[m] = np.sqrt((s_tot[:, None] - sq)[m] / n_loo[m])
        testable = seen & np.isfinite(sd_loo) & (np.nan_to_num(sd_loo) > 0)
        z_loo = np.zeros(sq.shape)
        z_loo[testable] = np.abs(dclr[testable]) / sd_loo[testable]
        ref['genera_changed_loo'] = (z_loo > a.change_z).sum(axis=0)
    dsh = ref['delta_shannon'].dropna().to_numpy()
    sh_sd = float(np.sqrt(np.mean(dsh ** 2))) if len(dsh) else np.nan
    loo_counts = ref['genera_changed_loo'].to_numpy(dtype=float)
    same_run_share = float(ref['same_run'].mean()) if len(ref) else np.nan

    out = ref.copy()
    for col in ('aitchison', 'trimmed_distance', 'delta_shannon'):
        out[col] = out[col].round(4)
    out.to_csv('%s.within_person_reference.tsv' % a.prefix, sep='\t', index=False)

    summ = []
    for metric in ('aitchison', 'trimmed_distance'):
        groups = [('within_person', 'all', ref[metric])] + \
                 [('within_person', s, g[metric]) for s, g in ref.groupby('study')] + \
                 [('between_person', 'all', btw[metric])] + \
                 [('between_person', s, g[metric]) for s, g in btw.groupby('study')]
        for scope, study, values in groups:
            v = np.asarray(values, dtype=float)
            summ.append({'metric': metric, 'reference': scope, 'study': study, 'n_pairs': len(v),
                         'median': rnd(np.median(v)) if len(v) else '',
                         'q05': rnd(np.percentile(v, 5)) if len(v) else '',
                         'q95': rnd(np.percentile(v, 95)) if len(v) else '',
                         'max': rnd(v.max()) if len(v) else ''})
    pd.DataFrame(summ).to_csv('%s.within_person_summary.tsv' % a.prefix, sep='\t', index=False)

    # --- your pairs -------------------------------------------------------------------
    changes, gchanges, ident = [], [], []
    pairs = list(itertools.combinations(personal_ids, 2))
    for x, y in pairs:
        d, t = dist(x, y), trimmed(x, y)
        same_run = sheet.at[x, 'run'] == sheet.at[y, 'run']
        caveat = ''
        if not have_ref:
            caveat = ('fewer than %d reference pairs: no within-person reference (add repeat samples)'
                      % a.min_pairs)
        elif not same_run and same_run_share > 0.5:
            caveat = ('sequenced in different runs, while %.0f%% of the reference pairs share a run: '
                      'part of this difference may be run-to-run technical variation'
                      % (100 * same_run_share))
        dsx = float(shannon[y] - shannon[x])
        n_changed, n_tested = 0, 0
        if have_ref:
            for gi, g in enumerate(genera):
                if not (P.at[g, x] or P.at[g, y]) or not np.isfinite(sd[gi]) or sd[gi] == 0:
                    continue
                n_tested += 1
                delta = float(C.at[g, y] - C.at[g, x])
                z = delta / sd[gi]
                changed = abs(z) > a.change_z
                n_changed += int(changed)
                gchanges.append({'sample_a': x, 'sample_b': y, 'genus': g,
                                 'changed': 'yes' if changed else 'no',
                                 'direction': 'up' if delta > 0 else 'down' if delta < 0 else '',
                                 'z': rnd(z, 2), 'delta_clr': rnd(delta), 'ordinary_change_sd': rnd(sd[gi]),
                                 'rel_abundance_a': round(float(counts.at[g, x] / max(counts[x].sum(), 1)), 5),
                                 'rel_abundance_b': round(float(counts.at[g, y] / max(counts[y].sum(), 1)), 5),
                                 'reference_pairs': int(n_tot[gi])})
        sh_ok = np.isfinite(dsx) and len(dsh) >= a.min_pairs and sh_sd > 0
        changes.append({
            'sample_a': x, 'sample_b': y,
            'context_a': sheet.at[x, 'context'], 'context_b': sheet.at[y, 'context'],
            'same_run': same_run, 'aitchison': rnd(d),
            'within_person_percentile': rnd(full.within_pct(d), 1),
            'between_person_percentile': rnd(full.between_pct(d), 1),
            'reference_pairs': len(ref),
            'reference_within_median': rnd(np.median(full.within)) if have_ref else '',
            'reference_between_median': rnd(np.median(full.between)),
            'delta_shannon': rnd(dsx),
            'delta_shannon_z': rnd(dsx / sh_sd, 2) if sh_ok else '',
            'delta_shannon_percentile': rnd(pct(np.abs(dsh), abs(dsx)), 1) if sh_ok else '',
            'genera_tested': n_tested if have_ref else '',
            'genera_changed': n_changed if have_ref else '',
            'genera_changed_expected': rnd(np.median(loo_counts), 1) if have_ref else '',
            'genera_changed_percentile': rnd(pct(loo_counts, n_changed), 1) if have_ref else '',
            'caveat': caveat,
        })

        if full.ordinary(d):
            verdict = 'same_person'
        elif trim.ordinary(t):
            verdict = 'same_person_large_change'
        else:
            verdict = 'possibly_different_people' if have_ref else 'unclear'
        detail = ('distance %.3f (one person: <= %s or < %.3f); trimmed %.3f (<= %s or < %.3f)'
                  % (d, rnd(full.within_thr), full.between_thr, t, rnd(trim.within_thr), trim.between_thr))
        if not have_ref:
            detail += ('; no within-person reference, so only the between-person half of the rule '
                       'was applied')
        ident.append({'sample_a': x, 'sample_b': y, 'aitchison': rnd(d),
                      'within_person_percentile': rnd(full.within_pct(d), 1),
                      'between_person_percentile': rnd(full.between_pct(d), 1),
                      'trimmed_distance': rnd(t),
                      'trimmed_within_person_percentile': rnd(trim.within_pct(t), 1),
                      'trimmed_between_person_percentile': rnd(trim.between_pct(t), 1),
                      'verdict': verdict, 'detail': detail})

    # overall verdict: one person, one odd sample out, or unclear
    one = ('same_person', 'same_person_large_change')
    if not pairs:
        overall, why = 'not_assessable', 'needs at least two personal samples'
    else:
        verdicts = {(r['sample_a'], r['sample_b']): r['verdict'] for r in ident}
        odd = [s for s in personal_ids
               if all(v == 'possibly_different_people' for k, v in verdicts.items() if s in k)
               and all(v in one for k, v in verdicts.items() if s not in k)]
        if all(v in one for v in verdicts.values()):
            large = ['%s-%s' % k for k, v in verdicts.items() if v == 'same_person_large_change']
            overall = 'one_person'
            why = 'every pair of samples is consistent with one person' + (
                '; %s differ in a few genera by more than one person ordinarily does' % ', '.join(large)
                if large else '')
        elif len(odd) == 1 and len(personal_ids) > 2:
            overall, why = 'one_sample_differs', ('%s looks like a different person than the other '
                                                 'samples: check for a mix-up' % odd[0])
        elif all(v == 'possibly_different_people' for v in verdicts.values()):
            overall, why = 'possibly_different_people', 'no pair of samples looks like one person'
        else:
            overall, why = 'unclear', ', '.join('%s-%s %s' % (k[0], k[1], v) for k, v in verdicts.items())
    ident.append({'sample_a': 'all', 'sample_b': '', 'aitchison': '', 'within_person_percentile': '',
                  'between_person_percentile': '', 'trimmed_distance': '',
                  'trimmed_within_person_percentile': '', 'trimmed_between_person_percentile': '',
                  'verdict': overall, 'detail': why})

    pd.DataFrame(changes, columns=[
        'sample_a', 'sample_b', 'context_a', 'context_b', 'same_run', 'aitchison',
        'within_person_percentile', 'between_person_percentile', 'reference_pairs',
        'reference_within_median', 'reference_between_median', 'delta_shannon', 'delta_shannon_z',
        'delta_shannon_percentile', 'genera_tested', 'genera_changed', 'genera_changed_expected',
        'genera_changed_percentile', 'caveat']).to_csv('%s.personal_changes.tsv' % a.prefix, sep='\t', index=False)
    gc = pd.DataFrame(gchanges, columns=[
        'sample_a', 'sample_b', 'genus', 'changed', 'direction', 'z', 'delta_clr', 'ordinary_change_sd',
        'rel_abundance_a', 'rel_abundance_b', 'reference_pairs'])
    if len(gc):
        gc['abs_z'] = gc['z'].astype(float).abs()
        gc = gc.sort_values(['sample_a', 'sample_b', 'abs_z', 'genus'],
                            ascending=[True, True, False, True]).drop(columns='abs_z')
    gc.to_csv('%s.personal_genus_changes.tsv' % a.prefix, sep='\t', index=False)
    pd.DataFrame(ident).to_csv('%s.identity.tsv' % a.prefix, sep='\t', index=False)

    print('within-person reference: %d pairs from %d people%s; between-person: %d pairs'
          % (len(ref), len(people), '' if have_ref else ' (too few, reference not used)', len(btw)))
    for r, i in zip(changes, ident):
        print('  %s -> %s: distance %s (within-person pct %s, between-person pct %s), %s genera changed '
              '(expected %s); identity %s' % (r['sample_a'], r['sample_b'], r['aitchison'],
                                              r['within_person_percentile'], r['between_person_percentile'],
                                              r['genera_changed'], r['genera_changed_expected'], i['verdict']))
    print('identity: %s (%s)' % (overall, why))


if __name__ == '__main__':
    main()
