#!/usr/bin/env python3
"""
Batch-dominance gate (DECISIONS.md D12): which claims may be read as biology at all.

Three rules, each reported with its numbers:

study_banner           R2(study) in the panel above the banner threshold (from BATCH_DIAGNOSE):
                       every population-facing claim is orientation only, said at the top.
contrast:<variable>    a lifestyle or metadata contrast counts as a finding only if its
                       marginal R2 exceeds R2(study) in the same model and p < --alpha.
personal_positioning   the personal samples' centroid must lie further from the pooled panel
                       centroid than the furthest panel STUDY centroid does. Otherwise the
                       personal deviation is "within between-study spread": a lab or cohort
                       difference of that size is ordinary, so it cannot be read as biology.
                       Needs at least two panel studies; with one, it is not assessable.

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

    for model, sub in perm[perm['model'].str.startswith('study + ') | ~perm['model'].isin(['study'])].groupby('model'):
        v = model.replace('study + ', '')
        mine = sub[sub['term'] == v]
        if mine.empty or mine.iloc[0]['R2'] == '':
            add('contrast:' + v, 'not_testable', '', '', mine.iloc[0]['note'] if not mine.empty else '')
            continue
        r2v, pv = float(mine.iloc[0]['R2']), float(mine.iloc[0]['p'])
        st = sub[sub['term'] == 'study']
        if 'confounded' in str(mine.iloc[0]['note']) or st.empty:
            add('contrast:' + v, 'orientation_only', r2v, 'R2(%s) > R2(study)' % v,
                'not separable from study: %s' % mine.iloc[0]['note'])
            continue
        r2st = float(st.iloc[0]['R2'])
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
        study_dist = {s: float(np.linalg.norm(clr[[i for i in panel_ids if sheet.at[i, 'study'] == s]].mean(axis=1) - pooled))
                      for s in studies}
        spread = max(study_dist.values())
        personal_dist = float(np.linalg.norm(clr[personal_ids].mean(axis=1) - pooled))
        add('personal_positioning', 'interpretable' if personal_dist > spread else 'within_between_study_spread',
            round(personal_dist, 3), 'largest study-centroid distance %.3f' % spread,
            'study-centroid distances: ' + ', '.join('%s %.3f' % kv for kv in sorted(study_dist.items())))
        for sid in personal_ids:
            d = float(np.linalg.norm(clr[sid] - pooled))
            add('personal_distance:%s' % sid, 'info', round(d, 3), '', sheet.at[sid, 'context'] if 'context' in sheet.columns else '')

    out = pd.DataFrame(rows)
    out.to_csv('%s.gates.tsv' % a.prefix, sep='\t', index=False)
    for _, r in out.iterrows():
        if not r['gate'].startswith('personal_distance'):
            print('%-24s %-28s %s' % (r['gate'], r['status'], r['detail']))


if __name__ == '__main__':
    main()
