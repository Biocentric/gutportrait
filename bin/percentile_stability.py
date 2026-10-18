#!/usr/bin/env python3
"""
Which genus positions survive batch correction? (DESIGN.md: every number is computed on the
corrected AND the uncorrected panel; anything that flips between them is unstable.)

Each personal-sample genus percentile is put in a category, low (<= --low-pct), high
(>= --high-pct) or middle, once against the uncorrected panel and once against the
corrected one. A position is stable when both agree. An unstable extreme was produced, or
hidden, by between-study differences, and must not be reported as a finding.
"""

import argparse

import numpy as np
import pandas as pd


def category(p, lo, hi):
    return np.where(p <= lo, 'low', np.where(p >= hi, 'high', 'middle'))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--uncorrected', required=True)
    p.add_argument('--corrected', required=True)
    p.add_argument('--prefix', default='out')
    p.add_argument('--low-pct', type=float, default=10.0)
    p.add_argument('--high-pct', type=float, default=90.0)
    a = p.parse_args(argv)

    u = pd.read_csv(a.uncorrected, sep='\t').fillna({'context': ''})
    c = pd.read_csv(a.corrected, sep='\t')
    keys = ['genus', 'sample']
    m = u[keys + ['context', 'percentile']].merge(c[keys + ['percentile']], on=keys, how='outer',
                                                  suffixes=('_uncorrected', '_corrected'))
    m['category_uncorrected'] = category(m['percentile_uncorrected'], a.low_pct, a.high_pct)
    m['category_corrected'] = category(m['percentile_corrected'], a.low_pct, a.high_pct)
    m['stable'] = np.where(m['category_uncorrected'] == m['category_corrected'], 'yes', 'no')
    m['shift'] = (m['percentile_corrected'] - m['percentile_uncorrected']).round(1)
    m = m.sort_values(['sample', 'stable', 'genus'])
    m.to_csv('%s.percentile_stability.tsv' % a.prefix, sep='\t', index=False)

    ext = m[(m['category_uncorrected'] != 'middle') | (m['category_corrected'] != 'middle')]
    print('%d genus positions; %d extreme in at least one version, of which %d stable'
          % (len(m), len(ext), int((ext['stable'] == 'yes').sum())))


if __name__ == '__main__':
    main()
