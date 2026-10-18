#!/usr/bin/env python3
"""
Were your samples processed like the reference samples? Technical comparability from
ampliseq's read tracking.

A profile can look unusual for technical reasons: a sample that lost most of its reads at
quality filtering, merged poorly or came out mostly chimeric has a different composition
for reasons that have nothing to do with the gut. Every step of ampliseq's read tracking
(overall_summary.tsv) is therefore turned into a retention rate, and each of your samples
is placed in the panel's distribution of that rate.

METRICS  (all from overall_summary.tsv except the last two)
  input_reads           read pairs entering the run (cutadapt input, else DADA2 input)
  primer_retention      share of reads in which cutadapt found the primers
  filter_retention      share kept by DADA2 quality filtering
  denoise_retention     share of filtered reads denoised in both directions
  merge_rate            share of denoised pairs that merged
  chimera_retention     share of merged reads not removed as chimeras
  end_to_end_retention  final reads / input reads
  analysed_reads        reads in the ASV table, after ampliseq's own filters (COLLAPSE_GENUS)
  unassigned_fraction   share of reads without a genus call (COLLAPSE_GENUS)
For every metric lower is worse, except unassigned_fraction, where higher is worse.

RULE
A metric is `atypical` when, on its worse side, it is further than --max-z robust standard
deviations from the median: the modified z-score of Iglewicz and Hoaglin (1993),
0.6745 * (value - median) / MAD, beyond 3.5 by default. That must hold against the pooled
panel AND against every panel study with at least --min-study-n samples (the
between-study-spread rule of the deviation flags: a value is only atypical if no single
reference study would call it ordinary). Read counts are compared on a log scale. A rate
that differs from the median by less than 5 percentage points is never atypical, however
tightly the panel clusters: such a difference cannot move a profile. Metrics ampliseq did
not report (e.g. no merging for single-end data) are `not_evaluable`.

A robust z-score rather than a tail percentile: with nine metrics per sample, a 2.5% tail
would mark one sample in five as atypical by chance alone.

Atypical is not an error. It says that the sample's composition may partly reflect how it
was processed, which the report states next to every result for that sample.
"""

import argparse
import sys

import numpy as np
import pandas as pd

LOWER_IS_WORSE = ['input_reads', 'primer_retention', 'filter_retention', 'denoise_retention',
                  'merge_rate', 'chimera_retention', 'end_to_end_retention', 'analysed_reads']
HIGHER_IS_WORSE = ['unassigned_fraction']
COUNTS = ('input_reads', 'analysed_reads')
MIN_RATE_GAP = 0.05


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--summary', required=True, help="ampliseq's overall_summary.tsv")
    p.add_argument('--collapse-stats', required=True, help='collapse_stats.tsv from COLLAPSE_GENUS')
    p.add_argument('--samplesheet', required=True)
    p.add_argument('--prefix', default='out')
    p.add_argument('--max-z', type=float, default=3.5, help='modified z-score beyond which a metric is atypical')
    p.add_argument('--min-study-n', type=int, default=10)
    return p.parse_args(argv)


def number(series):
    """ampliseq writes '12,345', '93.8%' and 'NA'; return floats with NaN for missing."""
    s = series.astype(str).str.replace(',', '', regex=False).str.replace('%', '', regex=False).str.strip()
    return pd.to_numeric(s.where(~s.isin(['', 'NA', 'nan', 'None'])), errors='coerce')


def ratio(num, den):
    return (num / den.where(den > 0)).astype(float)


def metrics_table(summary, stats):
    col = {c: number(summary[c]) for c in summary.columns if c != 'sample'}

    def get(name):
        return col.get(name, pd.Series(np.nan, index=summary.index))

    inp = get('cutadapt_total_processed').fillna(get('DADA2_input'))
    pairs = get('denoisedPairs')
    if pairs.isna().all():
        pairs = pd.concat([get('denoisedF'), get('denoisedR')], axis=1).min(axis=1, skipna=False)
    m = pd.DataFrame({
        'input_reads': inp,
        'primer_retention': ratio(get('cutadapt_passing_filters'), get('cutadapt_total_processed')),
        'filter_retention': ratio(get('filtered'), get('DADA2_input')),
        'denoise_retention': ratio(pairs, get('filtered')),
        'merge_rate': ratio(get('merged'), pairs),
        'chimera_retention': ratio(get('nonchim'), get('merged')),
        'end_to_end_retention': ratio(get('nonchim'), inp),
    })
    m.index = summary['sample'].astype(str)
    m['analysed_reads'] = stats['total_reads'].reindex(m.index).astype(float)
    m['unassigned_fraction'] = stats['unassigned_fraction'].reindex(m.index).astype(float)
    return m


def pct(ref, v):
    ref = np.sort(np.asarray(ref, dtype=float))
    return 100.0 * np.searchsorted(ref, v, side='right') / len(ref)


def modified_z(ref, v):
    """Iglewicz-Hoaglin modified z-score of v against ref; +-inf when ref has no spread."""
    med = float(np.median(ref))
    mad = float(np.median(np.abs(ref - med)))
    if mad == 0:
        return 0.0 if v == med else float(np.sign(v - med)) * np.inf
    return 0.6745 * (v - med) / mad


def main(argv=None):
    a = parse_args(argv)
    sheet = pd.read_csv(a.samplesheet, dtype=str).fillna('').set_index('sample')
    if 'context' not in sheet.columns:
        sheet['context'] = ''
    summary = pd.read_csv(a.summary, sep='\t', dtype=str)
    if 'sample' not in summary.columns:
        sys.exit("FATAL: %s has no 'sample' column; is it ampliseq's overall_summary.tsv?" % a.summary)
    stats = pd.read_csv(a.collapse_stats, sep='\t', index_col=0)
    stats.index = stats.index.astype(str)
    m = metrics_table(summary, stats)
    missing = sorted(set(sheet.index) - set(m.index))
    if missing:
        sys.exit('FATAL: samplesheet samples missing from overall_summary.tsv: %s' % ', '.join(missing[:10]))
    m = m.loc[list(sheet.index)]

    table = m.copy()
    for c in table.columns:
        table[c] = table[c].round(0 if c in COUNTS else 4)
    table.insert(0, 'run', sheet['run'])
    table.insert(0, 'study', sheet['study'])
    table.insert(0, 'role', sheet['role'])
    table.index.name = 'sample'
    table.to_csv('%s.technical_metrics.tsv' % a.prefix, sep='\t')

    panel_ids = list(sheet.index[sheet['role'] == 'panel'])
    personal_ids = list(sheet.index[sheet['role'] == 'personal'])
    study = sheet.loc[panel_ids, 'study']
    studies = [s for s, n in study.value_counts().items() if n >= a.min_study_n]

    def worse(ref, v, metric):
        """(modified z, atypical against this reference) on the metric's worse side."""
        if metric in COUNTS:
            ref, v = np.log10(ref + 1), np.log10(v + 1)
        z = modified_z(ref, v)
        sign = 1 if metric in HIGHER_IS_WORSE else -1
        gap_ok = metric in COUNTS or abs(v - np.median(ref)) >= MIN_RATE_GAP
        return z, sign * z > a.max_z and gap_ok

    rows = []
    for sid in personal_ids:
        bad = []
        for metric in LOWER_IS_WORSE + HIGHER_IS_WORSE:
            v = m.at[sid, metric]
            ref = m.loc[panel_ids, metric].dropna().to_numpy(dtype=float)
            row = {'sample': sid, 'context': sheet.at[sid, 'context'], 'metric': metric,
                   'worse_when': 'higher' if metric in HIGHER_IS_WORSE else 'lower'}
            if not np.isfinite(v) or len(ref) == 0:
                row.update({'value': '', 'panel_percentile': '', 'panel_median': '', 'modified_z': '',
                            'atypical': 'not_evaluable', 'detail': 'not reported for this sample or the panel'})
                rows.append(row)
                continue
            z, pooled_bad = worse(ref, v, metric)
            per_study = {}
            for s in studies:
                sref = m.loc[study.index[study == s], metric].dropna().to_numpy(dtype=float)
                if len(sref):
                    per_study[s] = worse(sref, v, metric)
            atypical = pooled_bad and all(b for _, b in per_study.values())
            if atypical:
                bad.append(metric)
            k = 0 if metric in COUNTS else 4
            row.update({'value': round(float(v), k), 'panel_percentile': round(pct(ref, v), 1),
                        'panel_median': round(float(np.median(ref)), k),
                        'modified_z': round(z, 2) if np.isfinite(z) else str(z),
                        'atypical': 'yes' if atypical else 'no',
                        'detail': ('far outside the pooled panel and every study' if atypical else
                                   'within the usual range' if not pooled_bad else
                                   'far outside the pooled panel, but ordinary for at least one study')
                        + ('; modified z per study: %s' % ', '.join(
                            '%s %s' % (s, round(zz, 2) if np.isfinite(zz) else zz)
                            for s, (zz, _) in sorted(per_study.items())) if per_study else '')})
            rows.append(row)
        rows.append({'sample': sid, 'context': sheet.at[sid, 'context'], 'metric': 'overall',
                     'worse_when': '', 'value': len(bad), 'panel_percentile': '', 'panel_median': '',
                     'modified_z': '', 'atypical': 'yes' if bad else 'no',
                     'detail': ('technically atypical in: %s' % ', '.join(bad)) if bad
                     else 'processed like the reference samples'})

    out = pd.DataFrame(rows, columns=['sample', 'context', 'metric', 'worse_when', 'value', 'panel_percentile',
                                      'panel_median', 'modified_z', 'atypical', 'detail'])
    out.to_csv('%s.technical_qc.tsv' % a.prefix, sep='\t', index=False)
    for _, r in out[out['metric'] == 'overall'].iterrows():
        print('%s: %s' % (r['sample'], r['detail']))


if __name__ == '__main__':
    main()
