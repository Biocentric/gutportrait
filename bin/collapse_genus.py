#!/usr/bin/env python3
"""
Collapse an ampliseq ASV table to genus-level counts.

WHY THIS EXISTS RATHER THAN READING rel-table-6.tsv
---------------------------------------------------
ampliseq's `qiime2/rel_abundance_tables/rel-table-N.tsv` files are numbered by *column
position* in the taxonomy string, not by rank name. SBDI-GTDB emits both a `Domain` and a
`Kingdom` column (both "Bacteria"), which shifts every rank down one relative to SILVA:
under SBDI-GTDB `rel-table-6` is FAMILY, and ampliseq writes only tables 2-6, so no genus
table is produced at all.

So we collapse from the *named-rank* taxonomy file (`dada2/ASV_tax_species.*.tsv`), whose
columns are labelled Domain/Kingdom/Phylum/.../Genus/Species. Named ranks do not shift when
the database changes. Verified empirically on ampliseq 2.18.0 + sbdi-gtdb.

UNASSIGNED GENERA ARE NOT DROPPED
---------------------------------
ASVs without a genus call are aggregated into a single `Unassigned` row rather than
discarded, and the fraction they represent is written to the stats file. Silently dropping
them would inflate every other genus's relative abundance by an unreported amount -- which
is exactly the class of error this pipeline exists to avoid.
"""

import argparse
import sys

import pandas as pd


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--asv-table', required=True,
                   help='ampliseq dada2/ASV_table.tsv (ASV_ID + one column per sample)')
    p.add_argument('--taxonomy', required=True,
                   help='ampliseq dada2/ASV_tax_species.<db>.tsv (named rank columns)')
    p.add_argument('--rank', default='Genus',
                   help='rank column to collapse on (default: Genus)')
    p.add_argument('--id-col', default='ASV_ID', help='ASV id column name')
    p.add_argument('--prefix', default='out', help='output file prefix')
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    asv = pd.read_csv(args.asv_table, sep='\t')
    tax = pd.read_csv(args.taxonomy, sep='\t')

    for name, df, path in (('ASV table', asv, args.asv_table),
                           ('taxonomy', tax, args.taxonomy)):
        if args.id_col not in df.columns:
            sys.exit("FATAL: '%s' has no column '%s' (found: %s)"
                     % (path, args.id_col, ', '.join(map(str, df.columns[:8]))))

    if args.rank not in tax.columns:
        sys.exit("FATAL: taxonomy file has no rank column '%s'. Available: %s\n"
                 "       Named-rank columns are required -- do not pass a positional\n"
                 "       rel-table-N.tsv here, see the module docstring."
                 % (args.rank, ', '.join(map(str, tax.columns))))

    sample_cols = [c for c in asv.columns if c != args.id_col]
    if not sample_cols:
        sys.exit("FATAL: no sample columns found in %s" % args.asv_table)

    merged = asv.merge(tax[[args.id_col, args.rank]], on=args.id_col, how='left')

    n_asv_total = len(merged)
    unassigned_mask = merged[args.rank].isna() | (merged[args.rank].astype(str).str.strip() == '')
    n_asv_unassigned = int(unassigned_mask.sum())

    merged[args.rank] = merged[args.rank].where(~unassigned_mask, 'Unassigned')

    collapsed = merged.groupby(args.rank, dropna=False)[sample_cols].sum()
    collapsed.index.name = args.rank
    collapsed = collapsed.sort_index()
    collapsed.to_csv('%s.genus_counts.tsv' % args.prefix, sep='\t')

    # Per-sample accounting: how much of each sample ended up unassigned.
    totals = merged[sample_cols].sum()
    unassigned_reads = merged.loc[unassigned_mask, sample_cols].sum()
    stats = pd.DataFrame({
        'total_reads': totals,
        'unassigned_reads': unassigned_reads,
        'unassigned_fraction': (unassigned_reads / totals.replace(0, pd.NA)).round(5),
        'n_taxa': (collapsed > 0).sum(),
    })
    stats.index.name = 'sample'
    stats.to_csv('%s.collapse_stats.tsv' % args.prefix, sep='\t')

    worst = stats['unassigned_fraction'].max()
    print('collapsed %d ASVs -> %d %s-level taxa across %d samples'
          % (n_asv_total, len(collapsed), args.rank, len(sample_cols)))
    print('%d ASVs (%.1f%%) had no %s call; highest per-sample unassigned fraction: %.1f%%'
          % (n_asv_unassigned, 100.0 * n_asv_unassigned / max(n_asv_total, 1),
             args.rank, 100.0 * (worst if pd.notna(worst) else 0)))


if __name__ == '__main__':
    main()
