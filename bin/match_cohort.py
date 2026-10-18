#!/usr/bin/env python3
"""
Compared with people like you: a matched cohort from the panel (DECISIONS.md D12).

THE LADDER
----------
The personal samples' bands (age_cat, sex, bmi_cat, from the samplesheet) select panel
people who share them. Exact matching comes first; when it yields fewer than --min-matched
people, the match is widened one step at a time, in the order D12 fixed:
  exact               same age band, sex and BMI band
  adjacent_age        age within one band either side, same sex and BMI band
  drop_bmi            age within one band either side, same sex
  drop_sex            age within one band either side
The first level with at least --min-matched people is the matched cohort. If none reaches
it, there is no matched comparison: the ladder is still reported, with its counts, so it is
visible how far the panel's metadata goes. A band missing for the personal samples drops
that criterion from every level.

WHO COUNTS AS A MATCH
---------------------
Only a panel person whose band is KNOWN to satisfy the criterion. An empty cell never
matches. A cell may hold a set of possible bands separated by "|" (a study's inclusion
criteria, e.g. "30s|40s|50s" for "adults 30 to 50"); it matches only if every band in the
set satisfies the criterion, so a match is never a guess.

WHAT IS COMPARED
----------------
Per genus, the CLR percentile of each personal sample among the matched people, next to its
percentile among the whole panel, and whether its band (low / middle / high, at --low-pct
and --high-pct) changes. Alpha diversity likewise. The genus set, pseudocount and CLR are
those of TAXON_PERCENTILES (panel-wide prevalence filter), so the two percentiles differ
only in who is in the reference.

OUTPUTS
-------
cohort_ladder.tsv         one row per level: criteria, matched people, per study, chosen
matched_samples.tsv       the matched panel people and their bands
matched_percentiles.tsv   per genus and personal sample: matched vs whole-panel percentile
matched_alpha.tsv         per metric and personal sample: matched vs whole-panel percentile
"""

import argparse
import sys
from collections import Counter

import numpy as np
import pandas as pd

AGE = ['20s', '30s', '40s', '50s', '60s', '70+']
LEVELS = [
    ('exact', 'same age band, sex and BMI band', {'age': 0, 'sex': True, 'bmi': True}),
    ('adjacent_age', 'age within one band either side, same sex and BMI band', {'age': 1, 'sex': True, 'bmi': True}),
    ('drop_bmi', 'age within one band either side, same sex', {'age': 1, 'sex': True, 'bmi': False}),
    ('drop_sex', 'age within one band either side', {'age': 1, 'sex': False, 'bmi': False}),
]


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--counts', required=True, help='genus counts (uncorrected) for all samples')
    p.add_argument('--alpha', required=True, help='alpha_diversity.tsv from ALPHA_DIVERSITY')
    p.add_argument('--samplesheet', required=True)
    p.add_argument('--prefix', default='out')
    p.add_argument('--min-matched', type=int, default=50)
    p.add_argument('--min-prevalence', type=float, default=0.1)
    p.add_argument('--pseudocount', type=float, default=0.5)
    p.add_argument('--low-pct', type=float, default=10)
    p.add_argument('--high-pct', type=float, default=90)
    return p.parse_args(argv)


def bands(cell):
    """The set of possible bands in a cell; None when unknown."""
    vals = {v.strip() for v in str(cell).split('|') if v.strip()}
    return vals or None


def pct(ref, v):
    ref = np.sort(np.asarray(ref, dtype=float))
    return 100.0 * np.searchsorted(ref, v, side='right') / len(ref)


def category(p, lo, hi):
    return 'low' if p <= lo else 'high' if p >= hi else 'middle'


def personal_bands(sheet, personal_ids):
    """The person's band per variable: the most common non-empty value over their samples."""
    out, notes = {}, []
    for col, key in (('age_cat', 'age'), ('sex', 'sex'), ('bmi_cat', 'bmi')):
        vals = [v for v in sheet.loc[personal_ids, col] if v] if col in sheet.columns else []
        if not vals:
            out[key] = None
            notes.append('%s not given for the personal samples, so not matched on' % col)
            continue
        common = Counter(vals).most_common()
        out[key] = common[0][0]
        if len(common) > 1:
            notes.append('%s differs between the personal samples (%s); matched on %s'
                         % (col, ', '.join(sorted(set(vals))), out[key]))
    return out, notes


def allowed(level, me):
    """{variable: allowed set} for one ladder level, given the person's bands."""
    rule = level[2]
    out = {}
    if me['age'] in AGE:
        i = AGE.index(me['age'])
        out['age_cat'] = set(AGE[max(0, i - rule['age']):i + rule['age'] + 1])
    if rule['sex'] and me['sex']:
        out['sex'] = {me['sex']}
    if rule['bmi'] and me['bmi']:
        out['bmi_cat'] = {me['bmi']}
    return out


def main(argv=None):
    a = parse_args(argv)
    sheet = pd.read_csv(a.samplesheet, dtype=str).fillna('').set_index('sample')
    for c in ('age_cat', 'sex', 'bmi_cat', 'context'):
        if c not in sheet.columns:
            sheet[c] = ''
    counts = pd.read_csv(a.counts, sep='\t', index_col=0)
    alpha = pd.read_csv(a.alpha, sep='\t').set_index('sample')
    panel_ids = list(sheet.index[sheet['role'] == 'panel'])
    personal_ids = list(sheet.index[sheet['role'] == 'personal'])
    missing = sorted(set(panel_ids + personal_ids) - set(counts.columns))
    if missing:
        sys.exit('FATAL: samplesheet samples missing from the counts table: %s' % ', '.join(missing[:10]))

    me, notes = personal_bands(sheet, personal_ids)
    ladder, chosen = [], None
    for level in LEVELS:
        rule = allowed(level, me)
        if not rule:
            ladder.append({'level': level[0], 'rule': level[1], 'criteria': 'none available', 'n_matched': 0,
                           'per_study': '', 'chosen': 'no'})
            continue
        members = [s for s in panel_ids
                   if all((b := bands(sheet.at[s, col])) is not None and b <= ok for col, ok in rule.items())]
        per_study = Counter(sheet.loc[members, 'study'])
        row = {'level': level[0], 'rule': level[1],
               'criteria': '; '.join('%s in {%s}' % (c, ', '.join(sorted(v, key=lambda x: AGE.index(x)
                                                                          if x in AGE else 0)))
                                     for c, v in rule.items()),
               'n_matched': len(members),
               'per_study': ', '.join('%s %d' % kv for kv in sorted(per_study.items())),
               'chosen': 'no'}
        if chosen is None and len(members) >= a.min_matched:
            row['chosen'] = 'yes'
            chosen = (level, members)
        ladder.append(row)
    status = ('matched at %s: %d people' % (chosen[0][0], len(chosen[1])) if chosen else
              'insufficient: no level reaches %d matched people' % a.min_matched)
    ladder = pd.DataFrame(ladder)
    ladder['status'] = status
    ladder['min_matched'] = a.min_matched
    ladder['personal_bands'] = ', '.join('%s %s' % (k, v) for k, v in
                                         (('age_cat', me['age']), ('sex', me['sex']), ('bmi_cat', me['bmi'])) if v)
    ladder['note'] = '; '.join(notes)
    ladder.to_csv('%s.cohort_ladder.tsv' % a.prefix, sep='\t', index=False)

    members = chosen[1] if chosen else []
    sheet.loc[members, ['study', 'age_cat', 'sex', 'bmi_cat']].rename_axis('sample').to_csv(
        '%s.matched_samples.tsv' % a.prefix, sep='\t')

    # percentiles: same genus set and CLR as TAXON_PERCENTILES, two references
    prev = (counts[panel_ids] > 0).mean(axis=1)
    kept = prev[prev >= a.min_prevalence].index
    logx = np.log(counts.loc[kept].astype(float) + a.pseudocount)
    clr = logx - logx.mean(axis=0)
    prow = []
    if members:
        for g in kept:
            if g == 'Unassigned':
                continue
            ref_all, ref_m = clr.loc[g, panel_ids].to_numpy(), clr.loc[g, members].to_numpy()
            for s in personal_ids:
                v = clr.at[g, s]
                p_all, p_m = pct(ref_all, v), pct(ref_m, v)
                c_all, c_m = category(p_all, a.low_pct, a.high_pct), category(p_m, a.low_pct, a.high_pct)
                prow.append({'genus': g, 'sample': s, 'context': sheet.at[s, 'context'],
                             'percentile_matched': round(p_m, 1), 'percentile_panel': round(p_all, 1),
                             'category_matched': c_m, 'category_panel': c_all,
                             'changed': 'yes' if c_m != c_all else 'no', 'shift': round(p_m - p_all, 1)})
    pd.DataFrame(prow, columns=['genus', 'sample', 'context', 'percentile_matched', 'percentile_panel',
                                'category_matched', 'category_panel', 'changed', 'shift']).to_csv(
        '%s.matched_percentiles.tsv' % a.prefix, sep='\t', index=False)

    arow = []
    if members:
        for m in ('shannon', 'observed', 'pielou'):
            vals = pd.to_numeric(alpha[m], errors='coerce')
            ref_all = vals.reindex(panel_ids).dropna().to_numpy()
            ref_m = vals.reindex(members).dropna().to_numpy()
            for s in personal_ids:
                v = vals.get(s, np.nan)
                if not np.isfinite(v) or not len(ref_m):
                    continue
                arow.append({'sample': s, 'context': sheet.at[s, 'context'], 'metric': m, 'value': round(float(v), 4),
                             'percentile_matched': round(pct(ref_m, v), 1),
                             'percentile_panel': round(pct(ref_all, v), 1),
                             'matched_median': round(float(np.median(ref_m)), 4), 'matched_n': len(ref_m)})
    pd.DataFrame(arow, columns=['sample', 'context', 'metric', 'value', 'percentile_matched', 'percentile_panel',
                                'matched_median', 'matched_n']).to_csv('%s.matched_alpha.tsv' % a.prefix,
                                                                       sep='\t', index=False)

    for _, r in ladder.iterrows():
        print('%-13s n=%-4d %s%s' % (r['level'], r['n_matched'], r['per_study'], '  <- chosen' if r['chosen'] == 'yes'
                                     else ''))
    print(status)
    if prow:
        changed = sum(r['changed'] == 'yes' for r in prow)
        print('%d of %d genus positions change band against the matched cohort' % (changed, len(prow)))


if __name__ == '__main__':
    main()
