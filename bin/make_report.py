#!/usr/bin/env python3
"""
The gutportrait report: every table the pipeline produced, in one document, in plain
language, as HTML and as PDF.

WHO IT IS FOR
-------------
Someone who knows what a gut microbiome is, but not what a centred log-ratio is. So every
section says, in this order: what it shows, your result in words, how to read it, and what
it does NOT mean. The caveats come first, before any result. A glossary closes the report.

CONFIDENCE LABELS
-----------------
Every section carries one, decided by the pipeline's own gates, not by taste:
  Solid             a technical result that is unlikely to change (identity, processing)
  Indicative        a real comparison with a known limit, stated next to it
  Orientation only  a position that differences between laboratories or studies could
                    explain (the batch-dominance gate decides this)
  Screening         a check for gross problems, provisional until a second, independent
                    comparison agrees

HOW IT IS BUILT
---------------
The sections are built once, as a list of simple blocks (headings, paragraphs, tables,
charts, call-outs), and rendered twice: to a self-contained HTML page (charts as inline
SVG, no external files or scripts) and to an A4 PDF typeset with matplotlib, with the same
charts. Nothing in either depends on the date, so both are reproducible.

Inline markup in the text: **bold**, [[n]] for the marker of personal sample n, and
<<level>> for a confidence label.
"""

import argparse
import glob
import html
import io
import itertools
import os
import re
import sys
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt                         # noqa: E402
from matplotlib import font_manager                     # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages    # noqa: E402
from matplotlib.ft2font import FT2Font                  # noqa: E402
from matplotlib.patches import Rectangle                # noqa: E402

matplotlib.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 8.5, 'svg.fonttype': 'none',
    'svg.hashsalt': 'gutportrait', 'pdf.fonttype': 42,
    'axes.spines.top': False, 'axes.spines.right': False, 'axes.spines.left': False,
    'axes.edgecolor': '#8c959f', 'xtick.color': '#57606a', 'ytick.color': '#57606a',
    'axes.labelcolor': '#57606a',
})

SAMPLE_COLOURS = ['#0072B2', '#E69F00', '#009E73', '#CC79A7', '#D55E00', '#56B4E9', '#000000']
REF_BAND, REF_DOT, REF_MEDIAN, INK, MUTED = '#d0d7de', '#afb8c1', '#57606a', '#1f2328', '#57606a'
LEVELS = {
    'solid': ('Solid', '#1a7f37'),
    'indicative': ('Indicative', '#0969da'),
    'orientation': ('Orientation only', '#9a6700'),
    'screening': ('Screening', '#6e7781'),
    'na': ('Not assessable', '#8c959f'),
}
ARMS = {'native': 'native V3–V4 region', 'v4projected': 'V4 region (projected)'}

FLAG_TEXT = {
    'low_alpha': ('Very low diversity', 'Shannon diversity below the {low}th percentile of the whole reference '
                  'panel and of every reference study'),
    'single_taxon_dominance': ('One genus dominating', 'a single genus making up more than {dom:.0%} of the sample'),
    'enterobacteriaceae_bloom': ('Enterobacteriaceae bloom', 'the family that includes E. coli above the {high}th '
                                 'percentile of the panel and of every study'),
    'proteobacteria_bloom': ('Proteobacteria bloom', 'the phylum Proteobacteria above the {high}th percentile of '
                             'the panel and of every study'),
    'faecalibacterium_depletion': ('Faecalibacterium very low', 'Faecalibacterium, one of the most common gut '
                                   'genera, below the {low}th percentile of the panel and of every study'),
    'core_genus_missing': ('Common genus missing', 'a genus found in at least {core:.0%} of reference people is '
                           'absent, although the sample was sequenced deeply enough to find it'),
    'beta_outlier': ('Profile far from everyone', 'the whole profile further from the average reference person '
                     'than {high_rest}% of reference people are'),
    'excess_extremes': ('Unusually many extreme genera', 'more genera in the outer tails than nearly every '
                        'reference person has'),
}
QC_TEXT = {
    'input_reads': ('Reads delivered', 'read pairs that entered the analysis'),
    'primer_retention': ('Primers found', 'reads that start with the expected primer, so come from the intended '
                         'gene region'),
    'filter_retention': ('Passed quality filter', 'reads good enough to keep'),
    'denoise_retention': ('Kept by error correction', 'filtered reads that error correction resolved in both '
                          'read directions'),
    'merge_rate': ('Read pairs joined', 'read pairs whose two halves overlapped and joined into one sequence'),
    'chimera_retention': ('Not chimeric', 'joined reads that are not PCR artefacts stitched from two sequences'),
    'end_to_end_retention': ('Kept overall', 'delivered reads that survived every step'),
    'analysed_reads': ('Reads analysed', 'reads the comparison is based on'),
    'unassigned_fraction': ('No genus name', 'reads that could not be named at genus level'),
}
GLOSSARY = [
    ('16S rRNA sequencing', 'Reading one gene that all bacteria carry (the 16S ribosomal RNA gene) from every '
     'bacterium in a sample at once. Its sequence tells which bacteria are present, usually down to the genus. '
     'It cannot see viruses or fungi and says nothing about what the bacteria are doing.'),
    ('Genus (plural genera)', 'A group of closely related bacterial species, such as Bacteroides or '
     'Faecalibacterium. This report works at the genus level, because 16S cannot reliably tell most species '
     'apart.'),
    ('ASV', 'Amplicon sequence variant: one exact 16S sequence found in the data, roughly one kind of '
     'bacterium, finer than a genus.'),
    ('Read', 'One DNA sequence produced by the sequencing machine. More reads give a more detailed look at a '
     'sample.'),
    ('Relative abundance', 'The share of a sample\'s reads that belong to a genus. Shares add up to 100%, so '
     'when one genus rises, all the others fall, even if nothing else changed.'),
    ('CLR (centred log-ratio)', 'A way to compare shares fairly: each genus is expressed relative to the '
     'average genus of the same sample, on a log scale. All percentiles and distances in this report use it.'),
    ('Percentile', 'The share of reference people with a value at or below yours. The 90th percentile means '
     'higher than 9 in 10 reference people; the 50th is the middle.'),
    ('Reference panel', 'The people you are compared with: participants of public studies whose raw '
     'sequencing data were processed in the same run as your samples.'),
    ('Study', 'A set of samples collected and sequenced together by one laboratory. Samples from one study '
     'resemble each other partly for technical reasons.'),
    ('Batch effect', 'A difference between studies caused by how samples were handled (collection kit, DNA '
     'extraction, sequencing), not by biology.'),
    ('Batch correction', 'A statistical adjustment (here MMUPHin) that removes study differences from the '
     'reference panel. Your own samples are never adjusted.'),
    ('Shannon diversity', 'A diversity measure that combines how many kinds of bacteria there are and how '
     'evenly they share the community.'),
    ('Observed ASVs', 'The number of different sequence variants found after subsampling every sample to the '
     'same number of reads: a count of kinds of bacteria.'),
    ('Evenness (Pielou)', 'How evenly the community is shared, from 0 (one kind takes nearly everything) to 1 '
     '(all kinds equally common).'),
    ('Rarefaction', 'Subsampling every sample to the same number of reads before measuring diversity, so a '
     'sample that was sequenced more deeply does not look more diverse.'),
    ('Aitchison distance', 'A distance between two whole profiles, computed on CLR values: 0 means identical, '
     'larger means more different. Its units have no everyday meaning; compare it only with the reference '
     'distances shown next to it.'),
    ('Within-person variation', 'How much two samples of the same person differ. It is never zero: diet, '
     'stool consistency, sampling and the laboratory all move a profile from day to day.'),
    ('Prevalence', 'The share of reference people in whom a genus is found at all.'),
    ('Core genus', 'A genus found in nearly every reference person (at least 95% of them).'),
    ('Positive control', 'A test with a known answer: an effect is planted on purpose, to check that the '
     'analysis can still find it.'),
    ('PERMANOVA and R²', 'A test of how much of the variation between people is explained by one factor, '
     'such as the study they come from. R² is that share: 0.06 means 6%.'),
    ('Leave-one-out', 'Repeating a calculation with each reference person set aside in turn and treated as '
     'if they were you, to see how often a check fires for an ordinary person.'),
    ('Modified z-score', 'How many typical deviations a value lies from the middle of the reference values. '
     'Beyond 3.5 counts as far outside the usual range.'),
    ('Enterobacteriaceae, Proteobacteria', 'A bacterial family that includes E. coli, and the larger group '
     '(phylum) it belongs to. Normally a small part of the gut community; they can bloom after an infection, '
     'antibiotics or during gut inflammation.'),
    ('Faecalibacterium', 'One of the most common genera in the gut of healthy adults. It produces butyrate, a '
     'fatty acid that the cells lining the gut use as fuel.'),
    ('Chimera', 'An artefact of the copying (PCR) step in which two different 16S sequences are joined into '
     'one. Chimeras are removed before analysis.'),
]


# ------------------------------------------------------------------------------ inputs

def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--tables', required=True, help='directory with every table the pipeline produced')
    p.add_argument('--samplesheet', required=True)
    p.add_argument('--prefix', default='out')
    p.add_argument('--arm', default='native')
    p.add_argument('--pipeline-version', default='')
    p.add_argument('--taxonomy-db', default='silva')
    p.add_argument('--flag-low-pct', type=float, default=5)
    p.add_argument('--flag-high-pct', type=float, default=95)
    p.add_argument('--dominance', type=float, default=0.5)
    p.add_argument('--core-prevalence', type=float, default=0.95)
    p.add_argument('--stability-low-pct', type=float, default=10)
    p.add_argument('--stability-high-pct', type=float, default=90)
    p.add_argument('--change-z', type=float, default=2)
    p.add_argument('--qc-max-z', type=float, default=3.5)
    p.add_argument('--no-pdf', action='store_true', help='write the HTML report only')
    return p.parse_args(argv)


class Tables:
    """Every pipeline table, found by its file-name suffix; None when absent."""

    def __init__(self, directory):
        self.files = sorted(glob.glob(os.path.join(directory, '*.tsv')))

    def get(self, suffix, exclude=None):
        hits = [f for f in self.files if f.endswith(suffix) and not (exclude and exclude in os.path.basename(f))]
        if not hits:
            return None
        return pd.read_csv(hits[0], sep='\t', dtype=str, keep_default_na=False)


def num(x, default=np.nan):
    try:
        v = float(str(x).replace(',', ''))
    except (TypeError, ValueError):
        return default
    return v if np.isfinite(v) else default


def col(df, name):
    return pd.to_numeric(df[name], errors='coerce') if df is not None and name in df.columns else pd.Series(dtype=float)


# ------------------------------------------------------------------------------ words

def ordinal(v):
    n = int(round(v))
    suffix = 'th' if 10 <= n % 100 <= 20 else {1: 'st', 2: 'nd', 3: 'rd'}.get(n % 10, 'th')
    return '%d%s' % (n, suffix)


def share(v):
    if not np.isfinite(v):
        return '–'
    if v == 0:
        return '0%'
    if v < 0.0001:
        return '<0.01%'
    if v < 0.01:
        return '%.2f%%' % (100 * v)
    if v < 0.1:
        return '%.1f%%' % (100 * v)
    return '%.0f%%' % (100 * v)


def pvalue(p):
    if not np.isfinite(p):
        return 'p not available'
    return 'p < 0.001' if p < 0.001 else 'p = %.3f' % p


def thousands(v):
    return '{:,}'.format(int(round(v))) if np.isfinite(v) else '–'


def plural(n, word, many=None):
    return '%s %s' % (thousands(n) if isinstance(n, (int, np.integer)) or float(n).is_integer() else n,
                      word if n == 1 else (many or word + 's'))


def listing(items, last='and'):
    items = list(items)
    if len(items) <= 1:
        return ''.join(items)
    return '%s %s %s' % (', '.join(items[:-1]), last, items[-1])


def pct_range(values):
    vals = [v for v in values if np.isfinite(v)]
    if not vals:
        return '–'
    lo, hi = ordinal(min(vals)), ordinal(max(vals))
    return lo if lo == hi else '%s to %s' % (lo, hi)


# ------------------------------------------------------------------------------ blocks

@dataclass
class Heading:
    text: str
    level: int = 2
    anchor: str = ''
    badge: tuple = ()          # (level, reason)


@dataclass
class Para:
    text: str
    style: str = 'normal'      # normal, lead, small, meta


@dataclass
class Bullets:
    items: list


@dataclass
class Table:
    header: list
    rows: list
    widths: list
    align: list = field(default_factory=list)
    note: str = ''


@dataclass
class Figure:
    draw: object               # draw(fig, rect)
    height: float              # inches, at the full text width
    caption: str = ''


@dataclass
class Callout:
    kind: str                  # warn, info, note
    title: str
    items: list
    bullets: bool = False


@dataclass
class Details:
    summary: str
    blocks: list
    pdf_note: str = ''


def runs(text):
    """Split inline markup into (text, style) runs."""
    out = []
    for part in re.split(r'(\*\*.+?\*\*|\[\[\d+\]\]|<<\w+>>)', text):
        if not part:
            continue
        if part.startswith('**') and part.endswith('**'):
            out.append((part[2:-2], {'bold': True}))
        elif part.startswith('[['):
            n = int(part[2:-2])
            out.append((str(n), {'marker': n, 'atomic': True}))
        elif part.startswith('<<'):
            out.append((LEVELS[part[2:-2]][0], {'badge': part[2:-2], 'atomic': True}))
        else:
            out.append((part, {}))
    return out


def colour(n):
    return SAMPLE_COLOURS[(n - 1) % len(SAMPLE_COLOURS)]


# ------------------------------------------------------------------------------ charts

def axes_in(fig, rect, l, b, w, h):
    """An axes at relative position (l, b, w, h) inside rect (figure fractions)."""
    x0, y0, rw, rh = rect
    return fig.add_axes([x0 + l * rw, y0 + b * rh, w * rw, h * rh])


def range_rows(rows, xlabel='', vlines=(), shared=True):
    """Horizontal range chart.

    rows: dicts with 'label', 'values' (reference values) or 'stats' (q05, q50, q95), and
    'marks': list of (text, value, colour). The grey bar spans the middle 90% of the
    reference, the dark tick is its median, small dots are individual reference values.
    Layout is in inches, so labels keep their size whatever the chart height.
    """
    def draw(fig, rect):
        fw, fh = fig.get_size_inches()
        h_in, w_in = rect[3] * fh, rect[2] * fw
        top, bottom = 0.30 if vlines else 0.05, 0.30 + (0.18 if xlabel else 0)
        n = len(rows)
        row_in = (h_in - top - bottom) / n
        allv = []
        for r in rows:
            if r.get('values') is not None:
                allv += list(np.asarray(r['values'], dtype=float))
            if r.get('stats'):
                allv += [r['stats'][k] for k in ('q05', 'q50', 'q95', 'max') if np.isfinite(r['stats'].get(k, np.nan))]
            allv += [m[1] for m in r['marks'] if np.isfinite(m[1])]
        allv += [v[0] for v in vlines]
        for i, r in enumerate(rows):
            band_top = h_in - top - i * row_in - 0.34          # room for the row title above
            ax_h = min(0.42, row_in - 0.40 - (0 if shared else 0.24))   # and for its own ticks below
            ax = axes_in(fig, rect, 0.03, (band_top - ax_h) / h_in, 0.94, ax_h / h_in)
            vals = r.get('values')
            if vals is not None and len(vals):
                vals = np.asarray(vals, dtype=float)
                q05, q50, q95 = np.percentile(vals, [5, 50, 95])
            else:
                vals = None
                q05, q50, q95 = (r['stats'][k] for k in ('q05', 'q50', 'q95'))
            if shared and allv:
                lo, hi = min(allv), max(allv)
            else:
                here = [q05, q95] + [m[1] for m in r['marks'] if np.isfinite(m[1])]
                lo, hi = min(here), max(here)
            pad = (hi - lo) * 0.05 or 1
            lo, hi = lo - pad, hi + pad
            ax.set_xlim(lo, hi)
            ax.add_patch(Rectangle((q05, 0.33), q95 - q05, 0.34, color=REF_BAND, lw=0, zorder=1))
            if vals is not None:
                jitter = np.random.default_rng(7).uniform(0.38, 0.62, len(vals))
                ax.scatter(vals, jitter, s=4, color=REF_DOT, lw=0, zorder=2)
            ax.plot([q50, q50], [0.25, 0.75], color=REF_MEDIAN, lw=1.6, zorder=3)
            for v, text, style in vlines:
                c = '#cf222e' if style == 'limit' else MUTED
                ax.axvline(v, color=c, lw=1, ls='--', zorder=2)
                if i == 0:
                    frac = (v - lo) / (hi - lo)
                    ax.text(v, 1.0 + 0.30 / ax_h, text, fontsize=7.5, color=c, va='bottom',
                            ha='right' if frac > 0.7 else 'left' if frac < 0.3 else 'center',
                            transform=ax.get_xaxis_transform())
            # stagger markers that would touch: greedy over three levels
            levels, placed = (0.5, 0.86, 0.14), []
            for text, v, c in sorted(r['marks'], key=lambda m: m[1]):
                diam = (0.26 if len(text) > 2 else 0.19) / (0.94 * w_in) * (hi - lo)
                y = next((lv for lv in levels if all(abs(v - pv) >= diam or py != lv for pv, py in placed)), 0.5)
                placed.append((v, y))
                ax.scatter([v], [y], s=150 if len(text) <= 2 else 290, color=c, zorder=4, lw=0.8, edgecolor='white',
                           clip_on=False)
                ax.text(v, y, text, ha='center', va='center', fontsize=7 if len(text) <= 2 else 6.2,
                        color='white', fontweight='bold', zorder=5, clip_on=False)
            ax.set_ylim(0, 1)
            ax.set_yticks([])
            ax.tick_params(axis='x', labelsize=7.5, length=3)
            ax.set_title(r['label'], loc='left', fontsize=8.5, color=INK, pad=11)
            if shared and i < n - 1:
                ax.set_xticklabels([])
            if xlabel and i == n - 1:
                ax.set_xlabel(xlabel, fontsize=7.5)
    return draw


def percentile_dots(genera, values, low, high):
    """values: {genus: [(sample_no, percentile)]}; tails shaded."""
    def draw(fig, rect):
        ax = axes_in(fig, rect, 0.30, 0.10, 0.66, 0.86)
        n = len(genera)
        ax.axvspan(0, low, color='#fff1e5', lw=0, zorder=0)
        ax.axvspan(high, 100, color='#ddf4ff', lw=0, zorder=0)
        ax.axvline(50, color=REF_BAND, lw=1, zorder=1)
        samples = sorted({s for g in genera for s, _ in values[g]})
        offsets = {s: (k - (len(samples) - 1) / 2) * 0.22 for k, s in enumerate(samples)}
        for i, g in enumerate(genera):
            y = n - 1 - i
            ax.axhline(y, color='#eaeef2', lw=0.6, zorder=1)
            for s, v in values[g]:
                ax.scatter([v], [y + offsets[s]], s=34, color=colour(s), zorder=3, lw=0.5, edgecolor='white')
        ax.set_yticks(range(n))
        ax.set_yticklabels(list(reversed(genera)), fontsize=8, fontstyle='italic', color=INK)
        ax.tick_params(axis='y', length=0)
        ax.set_xlim(0, 100)
        ax.set_ylim(-0.7, n - 0.3)
        ax.set_xticks([0, low, 25, 50, 75, high, 100])
        ax.tick_params(axis='x', labelsize=7.5)
        ax.set_xlabel('percentile among the reference people', fontsize=7.5)
        ax.text(low / 2, n - 0.35, 'low', ha='center', va='bottom', fontsize=7, color='#bc4c00')
        ax.text((high + 100) / 2, n - 0.35, 'high', ha='center', va='bottom', fontsize=7, color='#0969da')
    return draw


# ------------------------------------------------------------------------------ the data

class Data:
    def __init__(self, a):
        self.a = a
        t = Tables(a.tables)
        self.sheet = pd.read_csv(a.samplesheet, dtype=str, keep_default_na=False).set_index('sample')
        for c in ('context', 'subject'):
            if c not in self.sheet.columns:
                self.sheet[c] = ''
        self.personal = list(self.sheet.index[self.sheet['role'] == 'personal'])
        self.panel = list(self.sheet.index[self.sheet['role'] == 'panel'])
        self.repeats = list(self.sheet.index[self.sheet['role'] == 'repeat'])
        self.no = {s: i + 1 for i, s in enumerate(self.personal)}
        self.studies = self.sheet.loc[self.panel, 'study'].value_counts().sort_index()
        self.collapse = t.get('.collapse_stats.tsv')
        self.counts = t.get('.genus_counts.tsv')
        self.pct = t.get('.taxon_percentiles.tsv', exclude='.corrected.')
        self.pct_params = t.get('.percentile_params.tsv', exclude='.corrected.')
        self.stability = t.get('.percentile_stability.tsv')
        self.batch = t.get('.batch_summary.tsv')
        self.permanova = t.get('.permanova.tsv')
        self.correction = t.get('.correction_summary.tsv')
        self.pc = t.get('.positive_control.tsv')
        self.gates = t.get('.gates.tsv')
        self.alpha = t.get('.alpha_diversity.tsv')
        self.alpha_pos = t.get('.alpha_position.tsv')
        self.alpha_params = t.get('.alpha_params.tsv')
        self.flags = t.get('.deviation_flags.tsv')
        self.calib = t.get('.flag_calibration.tsv')
        self.wp_ref = t.get('.within_person_reference.tsv')
        self.wp_sum = t.get('.within_person_summary.tsv')
        self.changes = t.get('.personal_changes.tsv')
        self.gchanges = t.get('.personal_genus_changes.tsv')
        self.identity = t.get('.identity.tsv')
        self.qc_metrics = t.get('.technical_metrics.tsv')
        self.qc = t.get('.technical_qc.tsv')

    def label(self, s, context=True):
        c = self.sheet.at[s, 'context'] if s in self.sheet.index else ''
        return 'Sample %d%s' % (self.no[s], ' (%s)' % c if context and c else '')

    def mark(self, s):
        return '[[%d]]' % self.no[s]

    def gate(self, name):
        if self.gates is None:
            return None
        hit = self.gates[self.gates['gate'] == name]
        return hit.iloc[0] if len(hit) else None

    def positioning(self):
        g = self.gate('personal_positioning')
        return (g['status'] if g is not None else 'not_assessable'), g

    def banner(self):
        g = self.gate('study_banner')
        return g is not None and g['status'] == 'BANNER', g

    def qc_atypical(self):
        """{sample: [metric, ...]} of atypical metrics."""
        out = {}
        if self.qc is None:
            return out
        for s in self.personal:
            sub = self.qc[(self.qc['sample'] == s) & (self.qc['metric'] != 'overall') & (self.qc['atypical'] == 'yes')]
            out[s] = list(sub['metric'])
        return out

    def stable_extremes(self):
        """{genus: {sample: 'high'|'low'>> for positions extreme before AND after correction."""
        out = {}
        if self.stability is None:
            return out
        st = self.stability[self.stability['genus'] != 'Unassigned']
        for _, r in st.iterrows():
            if r['category_uncorrected'] == r['category_corrected'] and r['category_uncorrected'] in ('high', 'low'):
                out.setdefault(r['genus'], {})[r['sample']] = r['category_uncorrected']
        return out


# ------------------------------------------------------------------------------ sections

def section_title(d):
    a = d.a
    n_gen = int(num(d.pct_params.iloc[0]['n_genera_tested'])) if d.pct_params is not None else None
    blocks = [Heading('Your gut microbiome, compared', 1, 'top'),
              Para('%s of yours, set against %s from %s, all processed together in one run of the same '
                   'software.' % (plural(len(d.personal), 'sample').capitalize(),
                                  plural(len(d.panel), 'reference person', 'reference people'),
                                  plural(len(d.studies), 'reference study', 'reference studies')), 'lead'),
              Para(' · '.join(x for x in [
                  'gutportrait %s' % a.pipeline_version if a.pipeline_version else 'gutportrait',
                  'compared on the %s' % ARMS.get(a.arm, a.arm),
                  '%d genera compared' % n_gen if n_gen else '',
                  'not medical advice'] if x), 'meta')]
    return blocks


def caveats(d):
    out = []
    status, g = d.positioning()
    if status == 'within_spread':
        out.append('**Your overall profile is within the range that differences between laboratories '
                   'alone can produce.** Where a genus of yours sits high or low against the reference '
                   'people, part of that may be your laboratory rather than your gut. Read every position '
                   'as orientation, not as a finding (see [Overall position](#position)).')
    elif status == 'interpretable':
        out.append('Your overall profile is further from the reference people than laboratory and study '
                   'differences alone can explain, so where you differ, it is more likely to be biology.')
    elif status.startswith('mixed'):
        out.append('**Only some of your samples lie further from the reference people than laboratory '
                   'differences can explain** (%s). Read positions as orientation.' % status.split(': ')[-1])
    elif d.gates is None:
        out.append('**The study-difference checks are missing from this run**, so how much laboratories differ '
                   'is unknown. Read every position as orientation.')
    else:
        out.append('**There is only one usable reference study**, so the pipeline cannot measure how much '
                   'laboratories differ. Read every position as orientation.')
    up, bg = d.banner()
    if up:
        out.append('**The reference studies differ strongly from each other**: %.0f%% of the variation '
                   'among the reference people is explained by which study they come from, above the '
                   'alarm level of %.0f%%. Every comparison with the panel is orientation only.'
                   % (100 * num(bg['value']), 100 * num(bg['threshold'])))
    out.append('**The reference people were not matched to you** on age, sex or body-mass index. The public '
               'studies publish little about each person, and matching is not yet part of the pipeline. The '
               'panel describes these %s, not people like you in particular.'
               % plural(len(d.panel), 'person', 'people'))
    atyp = d.qc_atypical()
    bad = sorted({m for v in atyp.values() for m in v}, key=list(QC_TEXT).index)
    if bad:
        who = [s for s, v in atyp.items() if v]
        out.append('**%s processed differently from the reference samples** at some steps (%s). Their '
                   'composition may partly reflect that (see [Technical comparability](#technical)).'
                   % ('Your samples were' if len(who) == len(d.personal) and len(who) > 1 else
                      listing([d.label(s, False) for s in who]) + (' was' if len(who) == 1 else ' were'),
                      listing([QC_TEXT[m][0].lower() for m in bad])))
    if d.flags is not None and len(d.flags):
        out.append('**The abnormality checks are provisional.** They come from one comparison, and the design '
                   'asks for a second, independent comparison (the V4-region arm) to agree before a flag '
                   'counts.')
    if d.changes is not None and len(d.changes) and 'different runs' in ' '.join(d.changes['caveat']):
        out.append('**Your samples were sequenced in different runs**, while the reference pairs were not. Part '
                   'of the change between your samples may be technical (see [Changes](#changes)).')
    return out


def section_first(d):
    blocks = [Heading('Read this first', 2, 'first'),
              Para('This report compares the bacteria in your stool samples with those of %s. Their raw '
                   'sequencing data were processed together with yours, in one run of the same software, so '
                   'the comparison is not distorted by different analysis methods. Their samples were '
                   'collected and sequenced by other laboratories, though, and that difference cannot be '
                   'fully removed afterwards. It is the main limit on everything below.'
                   % plural(len(d.panel), 'other person', 'other people')),
              Para('**What 16S sequencing can and cannot do.** It reads one gene that all bacteria carry. That '
                   'tells which bacteria are present and roughly in what proportions, usually down to the '
                   'genus. It cannot see viruses or fungi, cannot tell most species and strains apart, and '
                   'says nothing about what the bacteria are doing.'),
              Callout('warn', 'What this report does not do', [
                  'It does not diagnose anything, and nothing in it is medical advice.',
                  'It gives no health or wellness score: no validated score exists for this kind of data.',
                  'It gives no dietary advice.',
                  '"Nothing stands out" is the most reassuring result it can give, and it is not the same '
                  'as "healthy". 16S finds gross abnormalities far better than it certifies health.'], True),
              Callout('info', 'Caveats for this comparison', caveats(d), True),
              Para('**Confidence labels.** Every section carries one: <<solid>> a technical result that is '
                   'unlikely to change; <<indicative>> a real comparison with a known limit, stated next to '
                   'it; <<orientation>> a position that differences between laboratories or studies could '
                   'explain; <<screening>> a check for gross problems, provisional until confirmed.')]
    return blocks


def identity_answer(d):
    if d.identity is None or not len(d.identity):
        return 'Not assessed (needs at least two samples).', 'na'
    last = d.identity.iloc[-1]
    v = last['verdict']
    if v == 'one_person':
        large = d.identity[d.identity['verdict'] == 'same_person_large_change']
        txt = 'Yes. Every pair of samples is as alike as two samples of one person.'
        if len(large):
            txt += ' %s differ in a few genera by more than usual.' % listing(
                ['%s and %s' % (d.mark(r['sample_a']), d.mark(r['sample_b'])) for _, r in large.iterrows()])
        level = 'solid' if d.wp_ref is not None and len(d.wp_ref) else 'indicative'
        return txt, level
    if v == 'one_sample_differs':
        odd = [s for s in d.personal if s in last['detail']]
        return ('Probably not: %s looks like a different person\'s sample. Check for a mix-up before '
                'reading anything else.' % (d.mark(odd[0]) if odd else 'one sample')), 'solid'
    if v == 'possibly_different_people':
        return 'Possibly not: no pair of samples looks like one person. Check for a mix-up.', 'solid'
    if v == 'not_assessable':
        return 'Not assessed (needs at least two samples).', 'na'
    return 'Unclear: see [Changes](#changes).', 'indicative'


def flag_name(f):
    return FLAG_TEXT.get(f, (f, ''))[0]


def flags_answer(d):
    if d.flags is None or not len(d.flags):
        return 'Not assessed.', 'na'
    raised = d.flags[d.flags['raised'] == 'yes']
    if not len(raised):
        return ('Nothing stands out on any of the %d checks.' % d.flags['flag'].nunique()), 'screening'
    parts = []
    for f, grp in raised.groupby('flag', sort=False):
        rate = num(grp.iloc[0]['reference_fire_rate'])
        parts.append('%s in %s (fires for %s of ordinary people)'
                     % (flag_name(f).lower(), listing([d.mark(s) for s in grp['sample']]), share(rate)))
    weak = all(num(r) >= 0.1 for r in raised['reference_fire_rate'])
    head = 'Nothing that stands out: ' if weak else ''
    return head + '%s raised: %s.' % ('one check' if raised['flag'].nunique() == 1 else
                                      '%d checks' % raised['flag'].nunique(), '; '.join(parts)), 'screening'


def section_glance(d):
    rows = []
    t, lvl = identity_answer(d)
    rows.append(['Are your samples from one person?', t, '<<%s>>' % lvl])
    t, lvl = flags_answer(d)
    rows.append(['Does anything look grossly unusual?', t, '<<%s>>' % lvl])
    if d.alpha_pos is not None and len(d.alpha_pos):
        sh = d.alpha_pos[d.alpha_pos['metric'] == 'shannon']
        low = sh[sh['below_every_study_p%g' % d.a.flag_low_pct] == 'True'] if (
            'below_every_study_p%g' % d.a.flag_low_pct) in sh.columns else sh.iloc[0:0]
        rows.append(['How diverse is your gut community?',
                     'Shannon diversity at the %s percentile of the reference people: %s.'
                     % (pct_range(col(sh, 'panel_percentile')),
                        'below the usual range in %s' % listing([d.mark(s) for s in low['sample']]) if len(low)
                        else 'within the usual range'), '<<indicative>>'])
    status, g = d.positioning()
    rows.append(['How different are you from the reference people overall?',
                 {'within_spread': 'Within the range that laboratory and study differences alone produce.',
                  'interpretable': 'Further from them than laboratory and study differences can explain.'}.get(
                     status, 'Not assessable with fewer than two reference studies.' if status == 'not_assessable'
                     else 'Mixed: %s.' % status.split(': ')[-1]),
                 '<<%s>>' % ('indicative' if status == 'interpretable' else 'orientation')])
    ext = d.stable_extremes()
    n = len(d.personal)
    hi = sorted(g for g, v in ext.items() if len(v) == n and set(v.values()) == {'high'})
    lo = sorted(g for g, v in ext.items() if len(v) == n and set(v.values()) == {'low'})
    if d.stability is not None:
        txt = ('%s in the outer %g%% in every sample, before and after batch correction%s'
               % (plural(len(hi) + len(lo), 'genus', 'genera').capitalize(), d.a.stability_low_pct,
                  ': ' + '; '.join(x for x in ['high **%s**' % ', '.join(hi[:6]) if hi else '',
                                               'low **%s**' % ', '.join(lo[:6]) if lo else ''] if x) + '.'
                  if hi or lo else '.'))
        rows.append(['Which genera stand out?', txt, '<<%s>>' % ('indicative' if status == 'interpretable'
                                                                and not d.banner()[0] else 'orientation')])
    if d.changes is not None and len(d.changes) and str(d.changes.iloc[0]['within_person_percentile']):
        w = col(d.changes, 'within_person_percentile')
        over = d.changes[col(d.changes, 'genera_changed_percentile') >= 95]
        txt = ('About as much as two samples of one person usually differ (%s percentile).' % pct_range(w)
               if (w <= 95).all() else 'More than two samples of one person usually differ, for %s.'
               % listing(['%s→%s' % (d.mark(r['sample_a']), d.mark(r['sample_b']))
                          for _, r in d.changes[w > 95].iterrows()]))
        if len(over):
            txt += ' %s changed in more genera than %s of reference people.' % (
                listing(['%s→%s' % (d.mark(r['sample_a']), d.mark(r['sample_b'])) for _, r in over.iterrows()]),
                listing(['%s%%' % ordinal(num(r['genera_changed_percentile']))[:-2] for _, r in over.iterrows()]))
        rows.append(['How much did your microbiome change between samples?', txt, '<<indicative>>'])
    elif len(d.personal) > 1:
        rows.append(['How much did your microbiome change between samples?',
                     'Not assessable: no reference people sampled twice.', '<<na>>'])
    atyp = d.qc_atypical()
    if atyp:
        bad = sorted({m for v in atyp.values() for m in v}, key=list(QC_TEXT).index)
        rows.append(['Were your samples processed like the reference samples?',
                     'Yes, at every step.' if not bad else 'Not entirely: atypical in %s.'
                     % listing([QC_TEXT[m][0].lower() for m in bad]), '<<solid>>'])
    return [Heading('At a glance', 2, 'glance'),
            Table(['Question', 'Answer', 'Confidence'], rows, [0.30, 0.52, 0.18])]


def section_samples(d):
    rows = []
    reads = dict(zip(d.collapse['sample'], col(d.collapse, 'total_reads'))) if d.collapse is not None else {}
    qc_over = {}
    if d.qc is not None:
        for _, r in d.qc[d.qc['metric'] == 'overall'].iterrows():
            ms = [m for m in d.qc_atypical().get(r['sample'], [])]
            qc_over[r['sample']] = 'atypical: %s' % ', '.join(QC_TEXT[m][0].lower() for m in ms) if ms \
                else 'like the reference'
    for s in d.personal:
        rows.append(['[[%d]] Sample %d' % (d.no[s], d.no[s]), d.sheet.at[s, 'context'] or '–', s,
                     thousands(reads.get(s, np.nan)), qc_over.get(s, '–')])
    return [Heading('Your samples', 2, 'samples'),
            Para('Your samples are numbered in the order of your samplesheet. The numbered markers are used in '
                 'every chart.'),
            Table(['Sample', 'When', 'Sample name', 'Reads analysed', 'Processing'], rows,
                  [0.15, 0.30, 0.17, 0.14, 0.24], ['l', 'l', 'l', 'r', 'l'])]


def section_diversity(d):
    blocks = [Heading('Diversity', 2, 'diversity', ('indicative', 'compared with the whole panel and with '
                                                    'every study'))]
    if d.alpha is None or d.alpha_pos is None:
        return blocks + [Para('Not available.')]
    depth = int(num(d.alpha_params.iloc[0]['depth'])) if d.alpha_params is not None else None
    blocks += [Para('Diversity describes how many different kinds of bacteria a sample holds and how evenly '
                    'they share the space. Diverse communities are common in healthy adults, but there is no '
                    'healthy cut-off: the range among ordinary people is wide.'),
               Bullets(['**Shannon diversity** combines the number of kinds and how evenly they are spread. It is '
                        'the most widely used single measure.',
                        '**Observed ASVs** counts the different sequence variants: roughly, the kinds of '
                        'bacteria.',
                        '**Evenness** runs from 0 (one kind takes nearly everything) to 1 (all equally common).']),
               Para('Every sample was first subsampled to %s reads, so that a sample sequenced more deeply does '
                    'not look more diverse.' % (thousands(depth) if depth else 'the same number of')
                    + ' A percentile says where you sit: the 30th percentile means 30% of the reference people '
                    'have a value at or below yours.')]
    panel = d.alpha[d.alpha['role'] == 'panel']
    rows = []
    for metric, label in (('shannon', 'Shannon diversity'), ('observed', 'Observed ASVs'), ('pielou', 'Evenness')):
        sub = d.alpha_pos[d.alpha_pos['metric'] == metric]
        rows.append({'label': label, 'values': col(panel, metric).dropna().to_numpy(),
                     'marks': [(str(d.no[r['sample']]), num(r['value']), colour(d.no[r['sample']]))
                               for _, r in sub.iterrows()]})
    blocks.append(Figure(range_rows(rows, shared=False), 3.3,
                         'Grey bar: the middle 90%% of the %d reference people; dark tick: their median; small dots: '
                         'individual reference people; numbered markers: your samples.' % len(panel)))
    trows = []
    for s in d.personal:
        cells = ['[[%d]] Sample %d' % (d.no[s], d.no[s])]
        for metric, k in (('shannon', 2), ('observed', 0), ('pielou', 2)):
            r = d.alpha_pos[(d.alpha_pos['sample'] == s) & (d.alpha_pos['metric'] == metric)]
            if len(r):
                r = r.iloc[0]
                cells.append('%s (%s)' % (('%.' + str(k) + 'f') % num(r['value']), ordinal(num(r['panel_percentile']))))
            else:
                cells.append('not measured (too few reads)')
        trows.append(cells)
    med = {m: num(d.alpha_pos[d.alpha_pos['metric'] == m].iloc[0]['panel_median']) for m in ('shannon', 'observed', 'pielou')
           if len(d.alpha_pos[d.alpha_pos['metric'] == m])}
    trows.append(['Reference median', '%.2f' % med.get('shannon', np.nan), '%.0f' % med.get('observed', np.nan),
                  '%.2f' % med.get('pielou', np.nan)])
    blocks.append(Table(['', 'Shannon (percentile)', 'Observed ASVs (percentile)', 'Evenness (percentile)'],
                        trows, [0.22, 0.26, 0.26, 0.26], ['l', 'r', 'r', 'r']))
    sh = d.alpha_pos[d.alpha_pos['metric'] == 'shannon']
    lowcol = 'below_every_study_p%g' % d.a.flag_low_pct
    low = sh[sh[lowcol] == 'True'] if lowcol in sh.columns else sh.iloc[0:0]
    blocks.append(Para('Your Shannon diversity sits at the %s percentile of the reference people. %s'
                       % (pct_range(col(sh, 'panel_percentile')),
                          ('%s below the %sth percentile of every reference study, the pipeline\'s threshold for '
                           '"unusually low".' % (listing([d.mark(s) for s in low['sample']]) +
                                                 (' is' if len(low) == 1 else ' are'), int(d.a.flag_low_pct)))
                          if len(low) else 'None of your samples is below the %sth percentile of every reference '
                          'study, the pipeline\'s threshold for "unusually low".' % int(d.a.flag_low_pct))))
    blocks.append(Callout('note', 'What it does not mean', [
        'Lower diversity than most reference people is not a disease marker on its own. Diversity moves with '
        'diet, age, recent illness, medication and stool consistency, and differs between laboratories. Only a '
        'value below every study\'s usual range is flagged, and even that is a reason to look, not a '
        'diagnosis.']))
    return blocks


def panel_distances(d):
    """Distances of panel people and your samples to the average reference person (CLR)."""
    if d.counts is None or d.pct_params is None:
        return None, {}
    counts = d.counts.set_index(d.counts.columns[0]).apply(pd.to_numeric, errors='coerce').fillna(0)
    ids = [s for s in d.panel + d.personal if s in counts.columns]
    counts = counts[ids]
    pc = num(d.pct_params.iloc[0]['pseudocount'], 0.5)
    mp = num(d.pct_params.iloc[0]['min_prevalence'], 0.1)
    panel = [s for s in d.panel if s in counts.columns]
    prev = (counts[panel] > 0).mean(axis=1)
    kept = prev[prev >= mp].index
    logx = np.log(counts.loc[kept] + pc)
    clr = logx - logx.mean(axis=0)
    centre = clr[panel].mean(axis=1)
    dist = np.sqrt(((clr.sub(centre, axis=0)) ** 2).sum(axis=0))
    return dist[panel].to_numpy(), {s: float(dist[s]) for s in d.personal if s in dist.index}


def section_position(d):
    status, g = d.positioning()
    level = 'indicative' if status == 'interpretable' else 'orientation'
    blocks = [Heading('Overall position', 2, 'position', (level, 'decided by the batch-dominance gate'))]
    blocks.append(Para('This compares your whole profile at once, all genera together, with the average '
                       'reference person. The distance used (Aitchison distance) is 0 for identical profiles and '
                       'grows as profiles differ; its units only mean something next to the reference distances.'))
    ref, mine = panel_distances(d)
    rows = [d.gate('personal_positioning:%s' % s) for s in d.personal]
    bound = num(g['threshold']) if g is not None else np.nan
    if ref is not None:
        marks = [(str(d.no[s]), mine[s], colour(d.no[s])) for s in d.personal if s in mine]
        vlines = [(bound, 'limit for an ordinary person in another lab', 'limit')] if np.isfinite(bound) else []
        blocks.append(Figure(range_rows([{'label': 'Distance to the average reference person', 'values': ref,
                                          'marks': marks}], xlabel='Aitchison distance', vlines=vlines), 1.6,
                             'Grey bar and dots: the %d reference people, each measured against the average of them '
                             'all; numbered markers: your samples; dashed line: how far an ordinary person measured '
                             'in a different laboratory or cohort could be.' % len(ref)))
    if status in ('within_spread', 'interpretable') or status.startswith('mixed'):
        vals = [num(r['value']) for r in rows if r is not None]
        detail = g['detail'] if g is not None else ''
        m = re.search(r'q95 within-study distance ([\d.]+) \+ largest study shift ([\d.]+)', detail)
        blocks.append(Para('Your samples are %s from the average reference person. The limit, %.1f, adds up two '
                           'things an ordinary person could show: the spread within one study (%s; 95%% of reference '
                           'people are closer to their own study\'s average) and the largest difference between the '
                           'averages of two studies (%s), which is what a different laboratory or cohort can add.'
                           % (listing(['%.1f' % v for v in vals]) if vals else '–', bound,
                              '%.1f' % float(m.group(1)) if m else '–', '%.1f' % float(m.group(2)) if m else '–')))
        blocks.append(Para({'within_spread': '**All your samples are within that limit.** So your overall '
                                             'difference from the reference people is no larger than what laboratory '
                                             'and study differences alone produce, and the positions in the next '
                                             'section are orientation, not findings.',
                            'interpretable': '**Your samples are beyond that limit.** Your overall difference from '
                                             'the reference people is larger than laboratory and study differences '
                                             'alone produce.'}.get(status, 'Some of your samples are beyond that '
                                                                           'limit and some within it.')))
    elif d.gates is None:
        blocks.append(Para('The study-difference checks did not run, so this comparison is not assessable.'))
    else:
        blocks.append(Para('With fewer than two reference studies of enough people, the pipeline cannot measure how '
                           'much a laboratory can shift a profile, so this comparison is not assessable.'))
    blocks.append(Callout('note', 'What it does not mean', [
        'Being close to the average reference person is not "healthier", and being far from it is not '
        '"unhealthier". Healthy people differ a lot from each other; the average is a reference point, not a '
        'target.']))
    return blocks


def section_genera(d):
    status, _ = d.positioning()
    level = 'indicative' if status == 'interpretable' and not d.banner()[0] else 'orientation'
    blocks = [Heading('Genus by genus', 2, 'genera', (level, 'positions, not findings, while your profile is within '
                                                       'the laboratory spread' if level == 'orientation' else
                                                       'stable under batch correction'))]
    if d.pct is None:
        return blocks + [Para('Not available.')]
    lo_t, hi_t = d.a.stability_low_pct, d.a.stability_high_pct
    prev = num(d.pct_params.iloc[0]['min_prevalence']) if d.pct_params is not None else np.nan
    blocks.append(Para('For every genus found in at least %s of the reference people, its share in your sample '
                       'is placed among its shares in the reference people, after the log-ratio transformation '
                       '(CLR) that makes shares comparable. A high percentile means the genus takes up more of '
                       'your community than in most reference people.' % share(prev)))
    ext = None
    if d.flags is not None:
        ex = d.flags[d.flags['flag'] == 'excess_extremes']
        if len(ex):
            m = re.search(r'of (\d+) genera; a reference person has a median of ([\d.]+)', ex.iloc[0]['detail'])
            if m:
                ext = (int(m.group(1)), float(m.group(2)), [(d.no[r['sample']], int(num(r['value'])))
                                                             for _, r in ex.iterrows()])
    if ext:
        blocks.append(Para('**Everyone has some genera in the tails.** Of %d genera, a typical reference person has '
                           '%s in the outer %g%% on either side; your samples have %s. Extreme genera are therefore '
                           'expected, and one of them on its own says little.'
                           % (ext[0], ('%g' % ext[1]), d.a.flag_low_pct,
                              listing(['%d ([[%d]])' % (v, n) for n, v in ext[2]]))))
    stable = d.stable_extremes()
    if d.stability is not None:
        st = d.stability[d.stability['genus'] != 'Unassigned']
        share_stable = (st['stable'] == 'yes').mean() if len(st) else np.nan
        blocks.append(Para('Each position was computed twice: against the reference panel as it is, and after '
                           'batch correction removed the differences between the reference studies. **%s of all '
                           'positions fall in the same band both times** (low: %sth percentile or below; high: %sth '
                           'or above; otherwise middle). Only positions that hold both times are listed below.'
                           % (share(share_stable), int(lo_t), int(hi_t))))
    order = sorted(stable, key=lambda g: (-len(stable[g]), 0 if 'high' in stable[g].values() else 1, g))
    show = order[:16]
    if show:
        pct_u = d.stability.set_index(['genus', 'sample'])
        values = {g: [(d.no[s], num(pct_u.loc[(g, s), 'percentile_uncorrected'])) for s in d.personal
                      if (g, s) in pct_u.index] for g in show}
        blocks.append(Figure(percentile_dots(show, values, lo_t, hi_t), 0.27 * len(show) + 0.75,
                             'Genera whose position is extreme and stable in the most of your samples (up to 16, '
                             'most consistent first). Each dot is one of your samples, coloured as its marker.'))
        rows = []
        for g in show:
            cells = ['*%s*' % g]
            for s in d.personal:
                if (g, s) in pct_u.index:
                    r = pct_u.loc[(g, s)]
                    flag = stable[g].get(s)
                    v = '%s → %s' % (ordinal(num(r['percentile_uncorrected'])), ordinal(num(r['percentile_corrected'])))
                    cells.append('**%s**' % v if flag else v)
                else:
                    cells.append('–')
            n_hi = sum(v == 'high' for v in stable[g].values())
            n_lo = sum(v == 'low' for v in stable[g].values())
            cells.append(('high in %d' % n_hi if n_hi else '') + (', ' if n_hi and n_lo else '') +
                         ('low in %d' % n_lo if n_lo else ''))
            rows.append(cells)
        w = 0.62 / max(len(d.personal), 1)
        blocks.append(Table(['Genus'] + ['[[%d]]' % d.no[s] for s in d.personal] + ['Stable extreme'], rows,
                            [0.22] + [w] * len(d.personal) + [0.16],
                            ['l'] + ['r'] * len(d.personal) + ['l'],
                            'Percentile before → after batch correction; bold where the position is extreme both '
                            'times. The 0th percentile means lower than every reference person, the 100th at least '
                            'as high as every reference person.'))
    else:
        blocks.append(Para('No genus is extreme both before and after batch correction.'))
    allrows = []
    piv = d.pct.pivot(index='genus', columns='sample', values='percentile')
    for g in sorted(piv.index, key=str.lower):
        if g == 'Unassigned':
            continue
        allrows.append(['*%s*' % g] + [ordinal(num(piv.at[g, s])) if s in piv.columns else '–' for s in d.personal])
    w = 0.7 / max(len(d.personal), 1)
    blocks.append(Details('All %d genera, percentile before batch correction' % len(allrows),
                          [Table(['Genus'] + ['[[%d]]' % d.no[s] for s in d.personal], allrows,
                                 [0.3] + [w] * len(d.personal), ['l'] + ['r'] * len(d.personal))],
                          'The full list of genera is in the HTML report and in percentiles/*.taxon_percentiles.tsv.'))
    blocks.append(Callout('note', 'What it does not mean', [
        'A genus at the 2nd percentile is not "missing" or "deficient", and one at the 98th is not "overgrown". '
        'Healthy people differ widely in their genera, and for almost none is a healthy range known at this '
        'level. A position becomes interesting when it is extreme, stable, repeated in every sample, and beyond '
        'the laboratory spread; even then, it describes, it does not diagnose.']))
    return blocks


def section_checks(d):
    blocks = [Heading('Checks for gross abnormalities', 2, 'checks', ('screening', 'one comparison arm; provisional'))]
    if d.flags is None:
        return blocks + [Para('Not available.')]
    a = d.a
    blocks.append(Para('These checks look only for the extremes that 16S can see reliably: a collapsed community, '
                       'one genus taking over, a bloom of bacteria that are normally scarce, the absence of '
                       'bacteria nearly everyone has. A value counts as extreme only if it is extreme against the '
                       'whole reference panel **and** against every reference study on its own, so that one unusual '
                       'study cannot raise a flag.'))
    blocks.append(Para('**How often each check fires for ordinary people.** Every reference person was checked in '
                       'turn, as if they were you, against all the others. A check that fires for 1 in 5 ordinary '
                       'people is weak evidence on its own; one that fires for nobody is strong.'))
    fmt = dict(low=int(a.flag_low_pct), high=int(a.flag_high_pct), high_rest=int(a.flag_high_pct),
               dom=a.dominance, core=a.core_prevalence)
    rows = []
    for f in [x for x in FLAG_TEXT if x in set(d.flags['flag'])]:
        sub = d.flags[d.flags['flag'] == f].set_index('sample')
        cells = ['**%s**' % FLAG_TEXT[f][0], FLAG_TEXT[f][1].format(**fmt)]
        for s in d.personal:
            if s not in sub.index:
                cells.append('–')
                continue
            r = sub.loc[s]
            cells.append({'yes': '**raised**', 'no': 'no', 'not_evaluable': 'not measurable'}.get(r['raised'],
                                                                                                r['raised']))
        rate = num(sub.iloc[0]['reference_fire_rate'])
        cells.append(share(rate) if np.isfinite(rate) else '–')
        rows.append(cells)
    w = 0.24 / max(len(d.personal), 1)
    blocks.append(Table(['Check', 'Looks for'] + ['[[%d]]' % d.no[s] for s in d.personal] + ['Fires for ordinary '
                                                                                            'people'],
                        rows, [0.20, 0.42] + [w] * len(d.personal) + [0.14],
                        ['l', 'l'] + ['c'] * len(d.personal) + ['r']))
    raised = d.flags[d.flags['raised'] == 'yes']
    if len(raised):
        items = []
        for _, r in raised.iterrows():
            f = r['flag']
            what = r['detail'] if f == 'core_genus_missing' else '%s (threshold %s)' % (r['value'], r['threshold'])
            items.append('%s %s: %s.' % (d.mark(r['sample']), flag_name(f).lower(), what))
        blocks.append(Callout('info', 'Raised checks in detail', items, True))
    else:
        blocks.append(Para('**No check was raised for any of your samples.**'))
    blocks.append(Callout('note', 'What it does not mean', [
        'A raised check is a reason to look closer, not a diagnosis; read it next to how often it fires for '
        'ordinary people and next to when the sample was taken. No raised check does not mean "healthy": it '
        'means nothing grossly abnormal was visible to this method.']))
    return blocks


def section_changes(d):
    if len(d.personal) < 2:
        return []
    have_ref = d.wp_ref is not None and len(d.wp_ref) > 0 and d.changes is not None and \
        str(d.changes.iloc[0]['within_person_percentile']) != ''
    blocks = [Heading('Changes between your samples', 2, 'changes',
                      ('indicative', 'against people sampled twice in the reference studies') if have_ref else
                      ('na', 'no reference people sampled twice'))]
    blocks.append(Para('No two samples of one person are identical: diet, stool consistency, sampling and the '
                       'laboratory move a profile from day to day. So a change between your samples only means '
                       'something next to how much one person\'s samples ordinarily differ. That is measured on '
                       'reference people whose studies sampled them twice%s.'
                       % (' (%s)' % plural(len(d.wp_ref), 'pair') if have_ref else '')))
    t, _ = identity_answer(d)
    blocks.append(Para('**Are they one person?** %s' % t))
    if d.wp_sum is not None and have_ref:
        s = d.wp_sum[(d.wp_sum['metric'] == 'aitchison') & (d.wp_sum['study'] == 'all')].set_index('reference')

        def stats(scope):
            r = s.loc[scope]
            return {k: num(r[k]) for k in ('q05', 'median', 'q95', 'max')} | {'q50': num(r['median'])}
        marks = [('%d–%d' % (d.no[r['sample_a']], d.no[r['sample_b']]), num(r['aitchison']), '#24292f')
                 for _, r in d.changes.iterrows()]
        n_w, n_b = int(num(s.loc['within_person', 'n_pairs'])), int(num(s.loc['between_person', 'n_pairs']))
        blocks.append(Figure(range_rows([
            {'label': 'Two samples of one person (%s reference pairs)' % thousands(n_w),
             'values': col(d.wp_ref, 'aitchison').to_numpy(), 'marks': marks},
            {'label': 'Two different people (%s reference pairs)' % thousands(n_b), 'stats': stats('between_person'),
             'marks': marks}], xlabel='Aitchison distance between two samples'), 2.3,
            'Grey bars: the middle 90% of each kind of pair; dark tick: the median; markers: your pairs of samples '
            '(1–2 is sample 1 against sample 2).'))
    rows = []
    for _, r in d.changes.iterrows():
        pair = '%s→%s' % (d.mark(r['sample_a']), d.mark(r['sample_b']))
        if have_ref:
            dist = '%.1f (%s percentile)' % (num(r['aitchison']), ordinal(num(r['within_person_percentile'])))
            gen = '%s (about %s by chance)' % (r['genera_changed'], ('%g' % num(r['genera_changed_expected'])))
        else:
            dist, gen = r['aitchison'], 'not assessable'
        dsh = num(r['delta_shannon'])
        shan = '–' if not np.isfinite(dsh) else ('%+.2f' % dsh) + (
            ' (%s percentile)' % ordinal(num(r['delta_shannon_percentile'])) if str(r['delta_shannon_percentile'])
            else '')
        rows.append([pair, dist, shan, gen])
    blocks.append(Table(['Pair', 'Distance (percentile among one person\'s two samples)', 'Shannon change',
                         'Genera that changed more than usual'], rows, [0.12, 0.34, 0.22, 0.32],
                        ['l', 'r', 'r', 'r'],
                        'A percentile of 80 means the pair differs more than 80%% of the reference people\'s two '
                        'samples do. A genus "changed more than usual" when its change was more than %g times the '
                        'ordinary change of that genus between two samples of one person.' % d.a.change_z))
    if have_ref and d.gchanges is not None and len(d.gchanges):
        for _, r in d.changes.iterrows():
            sub = d.gchanges[(d.gchanges['sample_a'] == r['sample_a']) & (d.gchanges['sample_b'] == r['sample_b'])
                             & (d.gchanges['changed'] == 'yes')]
            if not len(sub):
                continue
            sub = sub.assign(absz=col(sub, 'z').abs()).sort_values('absz', ascending=False).head(10)
            grows = [['*%s*' % g['genus'], 'up' if g['direction'] == 'up' else 'down',
                      '%s → %s' % (share(num(g['rel_abundance_a'])), share(num(g['rel_abundance_b']))),
                      '%.1f×' % abs(num(g['z']))] for _, g in sub.iterrows()]
            blocks.append(Para('%s→%s **the largest changes**' % (d.mark(r['sample_a']), d.mark(r['sample_b'])),
                               'subhead'))
            blocks.append(Table(['Genus', 'Direction', 'Share of reads', 'Size (× its usual change)'], grows,
                                [0.36, 0.14, 0.26, 0.24], ['l', 'l', 'r', 'r']))
    notes = []
    if d.changes is not None:
        cav = ' '.join(d.changes['caveat'])
        if 'different runs' in cav:
            notes.append('Your samples were sequenced in different runs, while the reference pairs were not; part '
                         'of the difference may be run-to-run technical variation, so read the counts as an upper '
                         'bound.')
        if 'fewer reads' in cav:
            notes.append('Some of your samples have fewer reads than nearly all reference samples. A rare genus '
                         'absent from a sample was only counted as a change where the sample had enough reads to '
                         'find it.')
    notes.append('The reference pairs were taken at the intervals their studies chose. Over a longer interval, '
                 'larger changes are ordinary.')
    blocks.append(Callout('info', 'Read these changes with', notes, True))
    blocks.append(Callout('note', 'What it does not mean', [
        'A change larger than usual says that something differed between the two samples: diet, travel, '
        'illness, medication, stool consistency, or the laboratory run. It does not say which, and not whether it '
        'was good or bad.']))
    return blocks


def qc_value(metric, v):
    if not np.isfinite(v):
        return '–'
    return thousands(v) if metric in ('input_reads', 'analysed_reads') else share(v)


def section_technical(d):
    blocks = [Heading('Technical comparability', 2, 'technical', ('solid', 'from the read counts at each step'))]
    if d.qc is None:
        return blocks + [Para('Not available.')]
    blocks.append(Para('Before any comparison, the reads pass several cleaning steps. If one of your samples lost '
                       'far more reads at a step than the reference samples did, its composition may partly '
                       'reflect that. Each step is compared with the reference samples; a value is marked '
                       '**atypical** when it is far outside their usual range (a modified z-score beyond %g) against '
                       'the whole panel and against every reference study. Differences of less than 5 percentage '
                       'points are never marked.' % d.a.qc_max_z))
    rows = []
    for metric in QC_TEXT:
        sub = d.qc[d.qc['metric'] == metric].set_index('sample')
        if not len(sub):
            continue
        cells = ['**%s**' % QC_TEXT[metric][0], QC_TEXT[metric][1], qc_value(metric, num(sub.iloc[0]['panel_median']))]
        for s in d.personal:
            if s not in sub.index:
                cells.append('–')
                continue
            r = sub.loc[s]
            v = qc_value(metric, num(r['value']))
            cells.append('**%s (atypical)**' % v if r['atypical'] == 'yes' else v)
        rows.append(cells)
    w = 0.36 / max(len(d.personal), 1)
    blocks.append(Table(['Step', 'Share of …', 'Reference median'] + ['[[%d]]' % d.no[s] for s in d.personal],
                        rows, [0.18, 0.32, 0.14] + [w] * len(d.personal), ['l', 'l', 'r'] + ['r'] * len(d.personal)))
    atyp = d.qc_atypical()
    bad = sorted({m for v in atyp.values() for m in v}, key=list(QC_TEXT).index)
    common = set.intersection(*[set(v) for v in atyp.values()]) if atyp else set()
    if common and len(d.personal) > 1:
        names = [QC_TEXT[m][0].lower() for m in QC_TEXT if m in common]
        extra = ['%s also %s' % (d.mark(s), ', '.join(QC_TEXT[m][0].lower() for m in v if m not in common))
                 for s, v in atyp.items() if set(v) - common]
        blocks.append(Para('All your samples are atypical at the same steps (%s). That points to your laboratory\'s '
                           'protocol, for example its read length or sequencing chemistry, rather than to one bad '
                           'sample.%s' % (listing(names), (' Atypical in addition: %s.' % '; '.join(extra)) if extra
                                          else '')))
    elif bad:
        blocks.append(Para('Atypical steps: %s.' % '; '.join('%s %s' % (d.mark(s), ', '.join(QC_TEXT[m][0].lower()
                                                                                          for m in v))
                                                             for s, v in atyp.items() if v)))
    else:
        blocks.append(Para('**Every step is within the usual range for all your samples.**'))
    blocks.append(Callout('note', 'What it does not mean', [
        'Atypical is not an error, and it does not make the results wrong. It says the sample was handled '
        'differently, so differences between you and the reference people are more likely to be partly '
        'technical.']))
    return blocks


def section_trust(d):
    blocks = [Heading('How far comparisons with the panel can be trusted', 2, 'trust',
                      ('solid', 'the pipeline\'s own gates'))]
    if d.pc is not None and len(d.pc):
        r = d.pc.iloc[0]
        if r['verdict'] == 'PASS':
            blocks.append(Para('**Positive control: passed.** To check that batch correction does not erase real '
                               'differences, the genus *%s* was doubled on purpose in half of the reference people of '
                               'every study, and the correction was run again. %s of the planted difference survived '
                               '(%s); at least %s had to.' % (r['genus'], share(num(r['retained'])), pvalue(num(r['p'])),
                                                               share(num(r['retain_required'])))))
        else:
            blocks.append(Para('**Positive control: %s.** %s' % (r['verdict'].lower().replace('_', ' '), r['note'])))
    up, bg = d.banner()
    if bg is not None and bg['status'] != 'not_assessable':
        blocks.append(Para('**Differences between studies.** %.1f%% of the variation among the reference people is '
                           'explained by which study they come from (PERMANOVA R² = %s). %s'
                           % (100 * num(bg['value']), bg['value'],
                              'That is above the alarm level of %.0f%%, so every comparison with the panel is '
                              'orientation only.' % (100 * num(bg['threshold'])) if up else
                              'That is below the alarm level of %.0f%%: the studies resemble each other enough to be '
                              'pooled.' % (100 * num(bg['threshold'])))))
    if d.gates is not None:
        cs = d.gates[d.gates['gate'].str.startswith('contrast:')]
        if len(cs):
            rows = []
            for _, r in cs.iterrows():
                v = r['gate'].split(':', 1)[1]
                plain = {'age_cat': 'Age band', 'sex': 'Sex', 'bmi_cat': 'BMI band', 'country': 'Country'}.get(v, v)
                verdict = {'finding': 'a finding', 'orientation_only': 'orientation only',
                           'not_testable': 'not testable'}.get(r['status'], r['status'])
                why = r['detail']
                if 'confounded' in why or 'collinear' in why:
                    why = 'cannot be separated from which study a person comes from'
                elif r['status'] == 'not_testable':
                    why = 'too few people with this information, or only one value'
                elif r['status'] in ('orientation_only', 'finding'):
                    m = re.search(r'R2\(\w+\) = ([\d.]+) \(p = ([^)]+)\) vs R2\(study\) = ([\d.]+)', why)
                    if m:
                        why = ('explains %.1f%% of the variation, %s the study does (%.1f%%)'
                               % (100 * float(m.group(1)), 'more than' if r['status'] == 'finding' else 'less than',
                                  100 * float(m.group(3))))
                rows.append([plain, verdict, why])
            blocks.append(Para('**Can the reference panel say anything about lifestyle or person factors?** A factor '
                               'counts as a finding only if it explains more of the variation than the study does.'))
            blocks.append(Table(['Factor', 'Verdict', 'Why'], rows, [0.18, 0.18, 0.64]))
    return blocks


def section_glossary():
    return [Heading('Glossary', 2, 'glossary'),
            Table(['Term', 'Meaning'], [['**%s**' % t, m] for t, m in GLOSSARY], [0.24, 0.76])]


def section_methods(d):
    a = d.a
    rows = [['Reference panel', '; '.join('%s (%d people)' % (s, n) for s, n in d.studies.items())],
            ['People sampled twice (noise floor)', thousands(len(d.repeats))],
            ['Comparison', ARMS.get(a.arm, a.arm)]]
    if d.pct_params is not None:
        p = d.pct_params.iloc[0]
        rows += [['Genera compared', '%s of %s (found in at least %s of the reference people)'
                  % (p['n_genera_tested'], p['n_genera_panel_total'], share(num(p['min_prevalence'])))],
                 ['Pseudocount', p['pseudocount']]]
    if d.alpha_params is not None:
        p = d.alpha_params.iloc[0]
        rows.append(['Rarefaction', '%s reads, %s repeats, seed %s' % (thousands(num(p['depth'])), p['iterations'],
                                                                        p['seed'])])
    rows += [['Extreme percentile (checks)', 'at or below %g, at or above %g' % (a.flag_low_pct, a.flag_high_pct)],
             ['Stable-position bands', 'low at or below %g, high at or above %g' % (a.stability_low_pct,
                                                                                 a.stability_high_pct)],
             ['Changed genus', 'more than %g× its ordinary within-person change' % a.change_z],
             ['Atypical processing step', 'modified z-score beyond %g' % a.qc_max_z]]
    return [Heading('Methods and settings', 2, 'methods'),
            Para('Your samples and the reference samples were processed together by nf-core/ampliseq (DADA2 '
                 'sequence variants, %s taxonomy), then compared by gutportrait: genus counts, centred '
                 'log-ratio percentiles, rarefied diversity, PERMANOVA and MMUPHin batch correction of the panel '
                 'only, the positive-control and batch-dominance gates, leave-one-out calibration of every check, '
                 'and a within-person reference from people sampled twice. Every number in this report comes '
                 'from the tables in the results directory, which hold more detail than shown here.'
                 % {'silva': 'SILVA', 'sbdi-gtdb': 'SBDI-GTDB', 'gtdb': 'GTDB'}.get(a.taxonomy_db, a.taxonomy_db)),
            Table(['Setting', 'Value'], rows, [0.34, 0.66]),
            Para('Not medical advice. Nothing in this report is diagnostic.', 'small')]


def build(d):
    blocks = section_title(d) + section_first(d) + section_glance(d) + section_samples(d)
    blocks += section_diversity(d) + section_position(d) + section_genera(d) + section_checks(d)
    blocks += section_changes(d) + section_technical(d) + section_trust(d) + section_glossary()
    blocks += section_methods(d)
    return blocks


# ------------------------------------------------------------------------------ HTML

CSS = """
:root { --ink: #1f2328; --muted: #57606a; --line: #d0d7de; --soft: #f6f8fa; }
* { box-sizing: border-box; }
body { margin: 0; background: #ffffff; color: var(--ink);
  font: 16px/1.6 -apple-system, BlinkMacSystemFont, "Segoe UI", "Helvetica Neue", Arial, sans-serif; }
main { max-width: 860px; margin: 0 auto; padding: 32px 20px 64px; }
h1 { font-size: 2rem; line-height: 1.2; margin: 0 0 8px; }
h2 { font-size: 1.35rem; margin: 48px 0 12px; padding-top: 20px; border-top: 1px solid var(--line); }
h2 .badge { margin-left: 10px; vertical-align: middle; }
p { margin: 10px 0; }
p.lead { font-size: 1.15rem; color: var(--ink); }
p.subhead { margin-top: 18px; }
p.meta, p.small, .note, .caption { color: var(--muted); font-size: 0.85rem; }
.why { color: var(--muted); font-size: 0.85rem; margin: -6px 0 12px; }
nav { display: flex; flex-wrap: wrap; gap: 4px 16px; background: var(--soft); border: 1px solid var(--line);
  border-radius: 8px; padding: 10px 16px; margin: 20px 0; font-size: 0.92rem; }
nav a { color: #0969da; text-decoration: none; white-space: nowrap; }
.scroll { overflow-x: auto; }
a { color: #0969da; }
table { width: 100%; border-collapse: collapse; margin: 12px 0; font-size: 0.9rem; }
th { text-align: left; font-weight: 600; background: var(--soft); border-bottom: 2px solid var(--line); }
th, td { padding: 6px 8px; vertical-align: top; }
td { border-bottom: 1px solid #eaeef2; }
td.r, th.r { text-align: right; } td.c, th.c { text-align: center; }
.callout { border-left: 4px solid; border-radius: 6px; padding: 10px 14px; margin: 16px 0; }
.callout h4 { margin: 0 0 6px; font-size: 0.95rem; }
.callout ul { margin: 4px 0; padding-left: 20px; } .callout p { margin: 4px 0; }
.callout.warn { border-color: #bf8700; background: #fff8c5; }
.callout.info { border-color: #0969da; background: #ddf4ff; }
.callout.note { border-color: #8c959f; background: var(--soft); }
.badge { display: inline-block; color: #fff; border-radius: 999px; padding: 1px 9px; font-size: 0.75rem;
  font-weight: 600; white-space: nowrap; }
.dot { display: inline-block; min-width: 1.35em; height: 1.35em; line-height: 1.35em; border-radius: 50%;
  color: #fff; text-align: center; font-size: 0.75em; font-weight: 700; vertical-align: 0.1em; }
figure { margin: 16px 0; } figure svg { width: 100%; height: auto; display: block; }
figcaption { color: var(--muted); font-size: 0.82rem; margin-top: 4px; }
details { margin: 12px 0; } summary { cursor: pointer; color: #0969da; }
@media (max-width: 600px) { body { font-size: 15px; } main { padding: 20px 16px 48px; }
  table { font-size: 0.8rem; } th, td { padding: 5px 4px; } }
@media print { nav, details { display: none; } h2 { break-after: avoid; } figure, tr, .callout { break-inside: avoid; }
  body { font-size: 11pt; } }
"""


def html_runs(text):
    out = []
    for t, st in runs(text):
        if st.get('marker'):
            out.append('<span class="dot" style="background:%s">%s</span>' % (colour(st['marker']), t))
        elif st.get('badge'):
            out.append('<span class="badge" style="background:%s">%s</span>' % (LEVELS[st['badge']][1], html.escape(t)))
        else:
            e = html.escape(t)
            e = re.sub(r'\[([^\]]+)\]\(#(\w+)\)', r'<a href="#\2">\1</a>', e)
            e = re.sub(r'\*(\S[^*]*?)\*', r'<em>\1</em>', e)
            out.append('<strong>%s</strong>' % e if st.get('bold') else e)
    return ''.join(out)


def svg_of(fig_block, width=7.0):
    fig = plt.figure(figsize=(width, fig_block.height))
    fig_block.draw(fig, [0, 0, 1, 1])
    buf = io.StringIO()
    fig.savefig(buf, format='svg', metadata={'Date': None})
    plt.close(fig)
    svg = buf.getvalue()
    svg = svg[svg.index('<svg'):]
    svg = re.sub(r'<svg([^>]*?) width="[^"]+" height="[^"]+"', r'<svg\1', svg, count=1)
    return svg


def to_html(blocks, title):
    body, toc = [], [(b.anchor, b.text) for b in blocks if isinstance(b, Heading) and b.level == 2]

    def render(b):
        if isinstance(b, Heading):
            badge = ''
            if b.badge:
                lvl, why = b.badge
                badge = ' <span class="badge" style="background:%s">%s</span>' % (LEVELS[lvl][1], LEVELS[lvl][0])
            out = '<h%d id="%s">%s%s</h%d>' % (b.level, b.anchor, html.escape(b.text), badge, b.level)
            if b.badge and b.badge[1]:
                out += '<p class="why">Confidence: %s, %s.</p>' % (LEVELS[b.badge[0]][0].lower(), html.escape(b.badge[1]))
            if b.level == 1:
                return out
            return out
        if isinstance(b, Para):
            return '<p class="%s">%s</p>' % (b.style, html_runs(b.text))
        if isinstance(b, Bullets):
            return '<ul>%s</ul>' % ''.join('<li>%s</li>' % html_runs(i) for i in b.items)
        if isinstance(b, Table):
            al = b.align or ['l'] * len(b.header)
            head = ''.join('<th class="%s">%s</th>' % (al[i], html_runs(h)) for i, h in enumerate(b.header))
            rows = ''.join('<tr>%s</tr>' % ''.join('<td class="%s">%s</td>' % (al[i], html_runs(str(c)))
                                                   for i, c in enumerate(r)) for r in b.rows)
            note = '<p class="note">%s</p>' % html_runs(b.note) if b.note else ''
            return '<table><thead><tr>%s</tr></thead><tbody>%s</tbody></table>%s' % (head, rows, note)
        if isinstance(b, Figure):
            return '<figure>%s<figcaption>%s</figcaption></figure>' % (svg_of(b), html_runs(b.caption))
        if isinstance(b, Callout):
            inner = ('<ul>%s</ul>' % ''.join('<li>%s</li>' % html_runs(i) for i in b.items)) if b.bullets else \
                ''.join('<p>%s</p>' % html_runs(i) for i in b.items)
            return '<div class="callout %s"><h4>%s</h4>%s</div>' % (b.kind, html.escape(b.title), inner)
        if isinstance(b, Details):
            return '<details><summary>%s</summary>%s</details>' % (html.escape(b.summary),
                                                                   ''.join(render(x) for x in b.blocks))
        return ''

    for i, b in enumerate(blocks):
        body.append(render(b))
        if isinstance(b, Para) and b.style == 'meta':
            body.append('<nav>%s</nav>' % ''.join('<a href="#%s">%s</a>' % (a, html.escape(t)) for a, t in toc))
    return ('<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
            '<meta name="viewport" content="width=device-width, initial-scale=1">\n<title>%s</title>\n'
            '<style>%s</style>\n</head>\n<body>\n<main>\n%s\n</main>\n</body>\n</html>\n'
            % (html.escape(title), CSS, '\n'.join(body)))


# ------------------------------------------------------------------------------ PDF

class Pdf:
    """A small typesetter on matplotlib: A4 pages, wrapped text, tables, charts."""
    W, H = 8.27, 11.69
    ML, MR, MT, MB = 0.8, 0.8, 0.75, 0.8

    def __init__(self, path, title):
        self.pdf = PdfPages(path, metadata={'Title': title, 'Creator': 'gutportrait', 'CreationDate': None})
        reg = font_manager.findfont(font_manager.FontProperties(family='DejaVu Sans'))
        bold = font_manager.findfont(font_manager.FontProperties(family='DejaVu Sans', weight='bold'))
        self.fonts = {False: FT2Font(reg), True: FT2Font(bold)}
        self.cw = self.W - self.ML - self.MR
        self.fig, self.page, self.y = None, 0, 0

    # metrics
    def width(self, s, size, bold=False):
        f = self.fonts[bold]
        f.set_size(size, 72)
        f.set_text(s, 0.0)
        return f.get_width_height()[0] / 64.0 / 72.0

    def tok_width(self, tok, st, size):
        if st.get('marker'):
            return self.width(tok, size * 0.8, True) + size * 0.75 / 72
        if st.get('badge'):
            return self.width(tok, size * 0.8, True) + size * 0.9 / 72
        return self.width(tok, size, st.get('bold', False))

    def space(self, size):
        return self.width('a a', size) - self.width('aa', size)

    def wrap(self, text, size, width):
        tokens, space = [], False
        for t, st in runs(re.sub(r'\[([^\]]+)\]\(#\w+\)', r'\1', text)):
            if st.get('atomic'):
                tokens.append((t, st, space))
                space = False
                continue
            t = re.sub(r'\*(\S[^*]*?)\*', r'\1', t)
            parts = t.split(' ')
            for i, p in enumerate(parts):
                if i > 0:
                    space = True
                if p:
                    tokens.append((p, st, space))
                    space = False
        lines, cur, curw = [], [], 0.0
        sp_w = self.space(size)
        for tok, st, sp in tokens:
            w = self.tok_width(tok, st, size)
            add = w + (sp_w if sp and cur else 0)
            if cur and curw + add > width:
                lines.append(cur)
                cur, curw, sp, add = [], 0.0, False, w
            cur.append((tok, st, sp and bool(cur)))
            curw += add
        if cur:
            lines.append(cur)
        return lines

    # pages
    def new_page(self):
        self.finish_page()
        self.fig = plt.figure(figsize=(self.W, self.H))
        self.page += 1
        self.y = self.H - self.MT

    def finish_page(self):
        if self.fig is None:
            return
        self.fig.text(self.ML / self.W, 0.45 / self.H, 'gutportrait · not medical advice', fontsize=7, color=MUTED)
        self.fig.text((self.W - self.MR) / self.W, 0.45 / self.H, str(self.page), fontsize=7, color=MUTED, ha='right')
        self.pdf.savefig(self.fig)
        plt.close(self.fig)
        self.fig = None

    def ensure(self, h):
        if self.fig is None or self.y - h < self.MB:
            self.new_page()

    def text_at(self, x, y, s, size, bold=False, color=INK, **kw):
        self.fig.text(x / self.W, y / self.H, s, fontsize=size, fontweight='bold' if bold else 'normal',
                      color=color, va='baseline', ha='left', **kw)

    def draw_line(self, line, x, y, size, color=INK):
        sp_w = self.space(size)
        for tok, st, sp in line:
            if sp:
                x += sp_w
            w = self.tok_width(tok, st, size)
            if st.get('marker'):
                self.fig.text((x + w / 2) / self.W, (y + size * 0.3 / 72) / self.H, tok, fontsize=size * 0.75,
                              fontweight='bold', color='white', ha='center', va='center',
                              bbox=dict(boxstyle='circle,pad=0.18', fc=colour(st['marker']), ec='none'))
            elif st.get('badge'):
                self.fig.text((x + w / 2) / self.W, (y + size * 0.3 / 72) / self.H, tok, fontsize=size * 0.75,
                              fontweight='bold', color='white', ha='center', va='center',
                              bbox=dict(boxstyle='round,pad=0.25', fc=LEVELS[st['badge']][1], ec='none'))
            else:
                self.text_at(x, y, tok, size, st.get('bold', False), color)
            x += w

    def para(self, text, size=9.5, x=None, width=None, color=INK, gap=0.06, lh=1.42):
        x = self.ML if x is None else x
        width = self.cw if width is None else width
        lines = self.wrap(text, size, width)
        step = size * lh / 72
        for line in lines:
            self.ensure(step)
            self.y -= step
            self.draw_line(line, x, self.y + step * 0.25, size, color)
        self.y -= gap

    # blocks
    def heading(self, b):
        if b.level == 1:
            self.ensure(0.6)
            self.y -= 0.32
            self.text_at(self.ML, self.y, b.text, 20, True)
            self.y -= 0.12
            return
        self.ensure(1.3)
        self.y -= 0.22
        self.fig.add_artist(plt.Line2D([self.ML / self.W, (self.W - self.MR) / self.W],
                                       [(self.y + 0.06) / self.H] * 2, color='#d0d7de', lw=0.8))
        self.y -= 0.24
        self.text_at(self.ML, self.y, b.text, 13.5, True)
        if b.badge:
            lvl, why = b.badge
            x = self.ML + self.width(b.text, 13.5, True) + 0.15
            self.draw_line([(LEVELS[lvl][0], {'badge': lvl, 'atomic': True}, False)], x, self.y + 0.01, 9)
            self.y -= 0.06
            self.para('Confidence: %s, %s.' % (LEVELS[lvl][0].lower(), why), 8, color=MUTED, gap=0.02)
        self.y -= 0.06

    def table(self, b):
        size, pad = 8.2, 0.05
        widths = [w * self.cw for w in b.widths]
        al = b.align or ['l'] * len(b.header)

        def cell_lines(row, bold=False):
            return [self.wrap(('**%s**' % c) if bold and c and '[[' not in c and '<<' not in c else str(c),
                              size, widths[i] - 2 * pad) for i, c in enumerate(row)]
        step = size * 1.35 / 72
        head = cell_lines(b.header, True)
        head_h = max(len(c) for c in head) * step + 2 * pad

        def draw_row(cells, h, fill=None):
            top = self.y
            if fill:
                self.fig.add_artist(Rectangle((self.ML / self.W, (top - h) / self.H), self.cw / self.W, h / self.H,
                                              transform=self.fig.transFigure, color=fill, lw=0, zorder=0))
            x = self.ML
            for i, lines in enumerate(cells):
                yy = top - pad
                for line in lines:
                    yy -= step
                    lw = sum(self.tok_width(t, s, size) for t, s, _ in line) + \
                        sum(self.space(size) for _, _, sp in line if sp)
                    xx = x + pad if al[i] == 'l' else (x + widths[i] - pad - lw if al[i] == 'r' else
                                                       x + (widths[i] - lw) / 2)
                    self.draw_line(line, xx, yy + step * 0.25, size)
                x += widths[i]
            self.y = top - h
            self.fig.add_artist(plt.Line2D([self.ML / self.W, (self.W - self.MR) / self.W], [self.y / self.H] * 2,
                                           color='#d0d7de' if fill else '#eaeef2', lw=0.7))

        self.ensure(head_h + step * 3)
        self.y -= 0.04
        draw_row(head, head_h, '#f6f8fa')
        for row in b.rows:
            cells = cell_lines(row)
            h = max(len(c) for c in cells) * step + 2 * pad
            if self.y - h < self.MB:
                self.new_page()
                draw_row(head, head_h, '#f6f8fa')
            draw_row(cells, h)
        self.y -= 0.06
        if b.note:
            self.para(b.note, 7.8, color=MUTED)
        self.y -= 0.04

    def figure(self, b):
        cap_lines = len(self.wrap(b.caption, 7.8, self.cw)) if b.caption else 0
        h = b.height * (self.cw / 7.0)
        self.ensure(h + cap_lines * 7.8 * 1.4 / 72 + 0.1)
        self.y -= 0.05
        rect = [self.ML / self.W, (self.y - h) / self.H, self.cw / self.W, h / self.H]
        b.draw(self.fig, rect)
        self.y -= h
        if b.caption:
            self.para(b.caption, 7.8, color=MUTED)
        self.y -= 0.05

    def callout(self, b):
        size, pad, bar = 8.8, 0.1, 0.05
        inner = self.cw - 2 * pad - bar - (0.16 if b.bullets else 0)
        items = [self.wrap(i, size, inner) for i in b.items]
        step = size * 1.4 / 72
        h = 2 * pad + 0.2 + sum(len(i) * step + 0.03 for i in items)
        if h > self.y - self.MB:
            self.new_page()
        top = self.y - 0.04
        colours = {'warn': ('#fff8c5', '#bf8700'), 'info': ('#ddf4ff', '#0969da'), 'note': ('#f6f8fa', '#8c959f')}
        bg, edge = colours[b.kind]
        self.fig.add_artist(Rectangle((self.ML / self.W, (top - h) / self.H), self.cw / self.W, h / self.H,
                                      transform=self.fig.transFigure, color=bg, lw=0, zorder=0))
        self.fig.add_artist(Rectangle((self.ML / self.W, (top - h) / self.H), bar / self.W, h / self.H,
                                      transform=self.fig.transFigure, color=edge, lw=0, zorder=0))
        y = top - pad - 0.14
        x0 = self.ML + bar + pad
        self.text_at(x0, y, b.title, 9.2, True)
        y -= 0.06
        for lines in items:
            for k, line in enumerate(lines):
                y -= step
                if b.bullets and k == 0:
                    self.text_at(x0, y + step * 0.25, '•', size)
                self.draw_line(line, x0 + (0.16 if b.bullets else 0), y + step * 0.25, size)
            y -= 0.03
        self.y = top - h - 0.08

    def render(self, blocks):
        for b in blocks:
            if self.fig is None:
                self.new_page()
            if isinstance(b, Heading):
                self.heading(b)
            elif isinstance(b, Para):
                if b.style == 'subhead':
                    self.ensure(1.2)          # keep a sub-heading with the start of its table
                size = {'lead': 11.5, 'small': 8, 'meta': 8.2}.get(b.style, 9.5)
                self.para(b.text, size, color=MUTED if b.style in ('meta', 'small') else INK,
                          gap=0.12 if b.style == 'lead' else 0.06)
            elif isinstance(b, Bullets):
                for item in b.items:
                    self.ensure(9.5 * 1.42 / 72)
                    y0 = self.y
                    self.para(item, 9.5, x=self.ML + 0.2, width=self.cw - 0.2, gap=0.02)
                    self.text_at(self.ML + 0.05, y0 - 9.5 * 1.42 / 72 * 0.75, '•', 9.5)
                self.y -= 0.05
            elif isinstance(b, Table):
                self.table(b)
            elif isinstance(b, Figure):
                self.figure(b)
            elif isinstance(b, Callout):
                self.callout(b)
            elif isinstance(b, Details) and b.pdf_note:
                self.para(b.pdf_note, 8, color=MUTED)
        self.finish_page()
        self.pdf.close()


# ------------------------------------------------------------------------------ main

def main(argv=None):
    a = parse_args(argv)
    d = Data(a)
    if not d.personal:
        sys.exit('FATAL: the samplesheet has no personal samples')
    blocks = build(d)
    title = 'gutportrait report'
    with open('%s.report.html' % a.prefix, 'w', encoding='utf-8', newline='\n') as fh:
        fh.write(to_html(blocks, title))
    if not a.no_pdf:
        Pdf('%s.report.pdf' % a.prefix, title).render(blocks)
    print('report: %d sections for %d personal samples against %d panel samples'
          % (sum(isinstance(b, Heading) and b.level == 2 for b in blocks), len(d.personal), len(d.panel)))


if __name__ == '__main__':
    main()
