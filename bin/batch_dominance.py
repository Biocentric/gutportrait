#!/usr/bin/env python3
"""
Batch-dominance gate (DECISIONS.md D12): which claims may be read as biology at all.

Three rules, each reported with its numbers:

study_banner           R2(study) in the panel above the banner threshold (from BATCH_DIAGNOSE):
                       every population-facing claim is orientation only, said at the top.
contrast:<variable>    a lifestyle or metadata contrast counts as a finding only if its
                       marginal R2 exceeds R2(study) in the same model and p < --alpha.
personal_positioning   is a personal sample further from the pooled panel centroid than an
                       ORDINARY person in a lab- or cohort-SHIFTED study could be? By the
                       triangle inequality, a person within the 95th percentile of their own
                       study's spread (q95 of within-study distances), in a study whose centroid
                       is shifted by at most the largest shift seen in the panel (D_max), lies at
                       most q95 + D_max from the pooled centroid. Only a sample beyond that bound
                       is "interpretable"; otherwise it is within individual plus between-study
                       spread, and its deviation cannot be read as biology. Needs at least two
                       panel studies; with one, it is not assessable.

                       (An earlier version compared the personal centroid with the study
                       centroids directly. That compared one person with group averages, which
                       always sit close together, so nearly everyone came out "interpretable".)

Distances are Aitchison (Euclidean on CLR genus profiles, panel-only prevalence filter and
the same pseudocount as everywhere else).
"""

import argparse

import numpy as np
import pandas as pd


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--permanova', required=True)
    p.add_argument('--batch-summary', required=True)
    p.add_argument('--counts', required=True, help='genus counts (uncorrected) for all samples')
    p.add_argument('--samplesheet', required=True)
    p.add_argument('--prefix', default='out')
    p.add_argument('--alpha', type=float, default=0.05)
    p.add_argument('--min-prevalence', type=float, default=0.1)
    p.add_argument('--pseudocount', type=float, default=0.5)
    p.add_argument('--min-study-n', type=int, default=10)
    return p.parse_args(argv)


def main(argv=None):
    a = parse_args(argv)
    perm = pd.read_csv(a.permanova, sep='\t').fillna('')
    summ = pd.read_csv(a.batch_summary, sep='\t').iloc[0]
    sheet = pd.read_csv(a.samplesheet, dtype=str).fillna('').set_index('sample')
    counts = pd.read_csv(a.counts, sep='\t', index_col=0)[list(sheet.index)]
    rows = []

    def add(gate, status, value, threshold, detail):
        rows.append({'gate': gate, 'status': status, 'value': value, 'threshold': threshold, 'detail': detail})

    r2s = summ['r2_study']
    if pd.isna(r2s) or str(r2s) in ('', 'NA'):
        add('study_banner', 'not_assessable', '', summ['banner_threshold'], 'fewer than two panel studies')
    else:
        add('study_banner', 'BANNER' if str(summ['study_banner']).upper() == 'TRUE' else 'clear',
            r2s, summ['banner_threshold'],
            'R2(study) %s the banner threshold: population claims are %s'
            % ('exceeds' if str(summ['study_banner']).upper() == 'TRUE' else 'is below',
               'orientation only' if str(summ['study_banner']).upper() == 'TRUE' else 'not demoted by it'))

    def num(x):
        """A statistic as float, or None when PERMANOVA could not produce one (NA)."""
        try:
            v = float(x)
        except (TypeError, ValueError):
            return None
        return v if np.isfinite(v) else None

    for model, sub in perm[perm['model'].str.startswith('study + ') | ~perm['model'].isin(['study'])].groupby('model'):
        v = model.replace('study + ', '')
        mine = sub[sub['term'] == v]
        r2v = num(mine.iloc[0]['R2']) if not mine.empty else None
        pv = num(mine.iloc[0]['p']) if not mine.empty else None
        if r2v is None or pv is None:
            add('contrast:' + v, 'not_testable', '', '',
                (mine.iloc[0]['note'] if not mine.empty else '') or 'no test statistic (NA)')
            continue
        st = sub[sub['term'] == 'study']
        if 'confounded' in str(mine.iloc[0]['note']) or st.empty:
            add('contrast:' + v, 'orientation_only', r2v, 'R2(%s) > R2(study)' % v,
                'not separable from study: %s' % mine.iloc[0]['note'])
            continue
        r2st = num(st.iloc[0]['R2'])
        if r2st is None:
            add('contrast:' + v, 'not_testable', r2v, '', 'R2(study) not estimable in this model')
            continue
        finding = r2v > r2st and pv < a.alpha
        add('contrast:' + v, 'finding' if finding else 'orientation_only', r2v, 'R2(%s) > R2(study)' % v,
            'R2(%s) = %.4f (p = %.3g) vs R2(study) = %.4f' % (v, r2v, pv, r2st))

    panel_ids = list(sheet.index[sheet['role'] == 'panel'])
    personal_ids = list(sheet.index[sheet['role'] == 'personal'])
    prev = (counts[panel_ids] > 0).mean(axis=1)
    kept = prev[prev >= a.min_prevalence].index
    logx = np.log(counts.loc[kept].astype(float) + a.pseudocount)
    clr = logx - logx.mean(axis=0)
    pooled = clr[panel_ids].mean(axis=1)
    studies = [s for s, n in sheet.loc[panel_ids, 'study'].value_counts().items() if n >= a.min_study_n]
    if len(studies) < 2:
        add('personal_positioning', 'not_assessable', '', '',
            'needs at least two panel studies of >= %d samples to measure between-study spread' % a.min_study_n)
    else:
        shift, within = {}, []
        for s in studies:
            members = [i for i in panel_ids if sheet.at[i, 'study'] == s]
            centre = clr[members].mean(axis=1)
            shift[s] = float(np.linalg.norm(centre - pooled))
            within += list(np.sqrt(((clr[members].sub(centre, axis=0)) ** 2).sum(axis=0)))
        d_max, q95 = max(shift.values()), float(np.percentile(within, 95))
        bound = q95 + d_max
        detail = ('bound = q95 within-study distance %.3f + largest study shift %.3f (%s)'
                  % (q95, d_max, ', '.join('%s %.3f' % kv for kv in sorted(shift.items()))))
        verdicts = []
        for sid in personal_ids:
            d = float(np.linalg.norm(clr[sid] - pooled))
            ok = d > bound
            verdicts.append(ok)
            ctx = sheet.at[sid, 'context'] if 'context' in sheet.columns else ''
            add('personal_positioning:%s' % sid, 'interpretable' if ok else 'within_spread', round(d, 3),
                '%.3f' % bound, ('%s; ' % ctx if ctx else '') + detail)
        overall = ('interpretable' if all(verdicts) else 'within_spread' if not any(verdicts)
                   else 'mixed: %d of %d samples interpretable' % (sum(verdicts), len(verdicts)))
        add('personal_positioning', overall, '', '%.3f' % bound, detail)

    out = pd.DataFrame(rows)
    out.to_csv('%s.gates.tsv' % a.prefix, sep='\t', index=False)
    for _, r in out.iterrows():
        print('%-36s %-28s %s' % (r['gate'], r['status'], r['value']))


if __name__ == '__main__':
    main()
