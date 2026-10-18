#!/usr/bin/env python3
"""
Generate the synthetic test fixture used by `-profile test`.

WHY SYNTHETIC
ampliseq's own test data is too small for this pipeline: three samples cannot exercise a
PERMANOVA, a percentile or a matched cohort, and ANCOM-BC already fails on it. Real panel
data cannot ship in a public repo. So the fixture is generated, deterministically, with
known signals planted in it. Tests can then check that the pipeline finds what was put
there, and only that.

WHAT IS PLANTED  (every later stage has a known answer to recover)
  1. AGE EFFECT          Older age bands get a more even composition, so Shannon diversity
                         rises with age_cat. This is the real-effect positive control (D12 b).
  2. STUDY EFFECT        STUDY_C carries a multiplicative bias on three genera (Bacteroides
                         x3, Segatella x0.3, Bifidobacterium x2), a batch effect for the
                         batch-dominance rule to detect.
  3. DEVIANT SAMPLE      Personal sample self_P3 has an Escherichia-Shigella bloom (~30%), for
                         the Proteobacteria deviation flag.
  4. MISSING METADATA    STUDY_C has no per-sample age/sex/BMI, mirroring the ENA panels,
                         which publish none (see TODO.md section 1b).
  5. ONE PERSON          The three personal samples are one person sampled three times
                         (within-person noise 0.15 against 0.45 between people), for the
                         identity check. self_P3's bloom is a within-person change that must
                         stand out from the borrowed noise floor.
  6. REPEAT SAMPLES      16 panel people have a second sample (role `repeat`, `subject`
                         pointing at the first), the reference for within-person variation.
  7. TECHNICAL OUTLIER   self_P1 loses half its reads at quality filtering, in
                         overall_summary.tsv, for the technical-comparability check.

The personal rows are invented. Their metadata bands are deliberately not anyone's real
ones: this repository is public, and a fixture is no place for a real person's profile.

OUTPUT LAYOUT  mirrors an ampliseq 2.18.0 results directory, so the pipeline reads it the
same way it reads real output:
  synthetic/samplesheet.csv  (and variants: _collinear, _norepeats, _mixup)
  synthetic/ampliseq_results/overall_summary.tsv               (ampliseq's read tracking)
  synthetic/ampliseq_results/dada2/ASV_table.tsv
  synthetic/ampliseq_results/dada2/ASV_tax_species.silva.tsv   (SILVA column layout)

Re-generate with:  python tests/fixtures/make_synthetic_fixture.py
The seed is fixed; the output is byte-identical on every run.
"""

import csv
import hashlib
import math
import os
import random

SEED = 20261015
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'synthetic')

# genus -> (phylum, class, order, family, baseline log-abundance)
GENERA = {
    'Bacteroides':                   ('Bacteroidota', 'Bacteroidia', 'Bacteroidales', 'Bacteroidaceae', 3.0),
    'Segatella':                     ('Bacteroidota', 'Bacteroidia', 'Bacteroidales', 'Prevotellaceae', 2.2),
    'Alistipes':                     ('Bacteroidota', 'Bacteroidia', 'Bacteroidales', 'Rikenellaceae', 1.6),
    'Parabacteroides':               ('Bacteroidota', 'Bacteroidia', 'Bacteroidales', 'Tannerellaceae', 1.2),
    'Faecalibacterium':              ('Bacillota', 'Clostridia', 'Oscillospirales', 'Ruminococcaceae', 2.6),
    'Blautia':                       ('Bacillota', 'Clostridia', 'Lachnospirales', 'Lachnospiraceae', 2.0),
    'Roseburia':                     ('Bacillota', 'Clostridia', 'Lachnospirales', 'Lachnospiraceae', 1.5),
    'Coprococcus':                   ('Bacillota', 'Clostridia', 'Lachnospirales', 'Lachnospiraceae', 1.1),
    'Dorea':                         ('Bacillota', 'Clostridia', 'Lachnospirales', 'Lachnospiraceae', 0.8),
    'Lachnospira':                   ('Bacillota', 'Clostridia', 'Lachnospirales', 'Lachnospiraceae', 0.9),
    'Subdoligranulum':               ('Bacillota', 'Clostridia', 'Oscillospirales', 'Ruminococcaceae', 1.3),
    'Ruminococcus':                  ('Bacillota', 'Clostridia', 'Oscillospirales', 'Ruminococcaceae', 1.0),
    'Christensenellaceae R-7 group': ('Bacillota', 'Clostridia', 'Christensenellales', 'Christensenellaceae', 0.9),
    'UCG-002':                       ('Bacillota', 'Clostridia', 'Oscillospirales', 'Oscillospiraceae', 0.8),
    'Dialister':                     ('Bacillota', 'Negativicutes', 'Veillonellales-Selenomonadales', 'Veillonellaceae', 0.9),
    'Veillonella':                   ('Bacillota', 'Negativicutes', 'Veillonellales-Selenomonadales', 'Veillonellaceae', 0.2),
    'Streptococcus':                 ('Bacillota', 'Bacilli', 'Lactobacillales', 'Streptococcaceae', 0.4),
    'Bifidobacterium':               ('Actinomycetota', 'Actinobacteria', 'Bifidobacteriales', 'Bifidobacteriaceae', 1.4),
    'Collinsella':                   ('Actinomycetota', 'Coriobacteriia', 'Coriobacteriales', 'Coriobacteriaceae', 0.9),
    'Akkermansia':                   ('Verrucomicrobiota', 'Verrucomicrobiae', 'Verrucomicrobiales', 'Akkermansiaceae', 0.7),
    'Escherichia-Shigella':          ('Pseudomonadota', 'Gammaproteobacteria', 'Enterobacterales', 'Enterobacteriaceae', 0.3),
    'Sutterella':                    ('Pseudomonadota', 'Gammaproteobacteria', 'Burkholderiales', 'Sutterellaceae', 0.4),
}
# ASVs with no genus call (the 'Unassigned' row must survive collapsing)
UNASSIGNED = [('Bacillota', 'Clostridia', 'Oscillospirales', '', 1.0),
              ('Bacillota', 'Clostridia', '', '', 0.6)]

AGE_CATS = ['20s', '30s', '40s', '50s', '60s']
STUDY_C_BIAS = {'Bacteroides': 3.0, 'Segatella': 0.3, 'Bifidobacterium': 2.0}
WITHIN_SD = 0.15            # within-person noise, against 0.45 between people
REPEATS = ['studya_%02d' % i for i in range(1, 7)] + ['studyb_%02d' % i for i in range(1, 7)] + \
          ['studyc_%02d' % i for i in range(1, 5)]
TECH_OUTLIER = 'self_P1'    # loses half its reads at quality filtering


def make_asvs(rng):
    """2-3 ASVs per genus; ids are md5 of the sequence, as ampliseq does."""
    asvs = []
    rows = [(g, v[:4], v[4]) for g, v in GENERA.items()] + \
           [('', u[:4], u[4]) for u in UNASSIGNED]
    for genus, ranks, base in rows:
        for k in range(rng.choice([2, 3])):
            seq = ''.join(rng.choice('ACGT') for _ in range(400))
            asvs.append({'id': hashlib.md5(seq.encode()).hexdigest(), 'seq': seq,
                         'genus': genus, 'ranks': ranks, 'base': base - 0.5 * k})
    return asvs


def latent(rng, asvs, age_idx, study=None):
    """A person's log-abundance profile: noise around the baseline, then the planted effects.

    Age: the log-abundance spread is compressed toward its mean by 12% per age band, so
    older samples are more even and have higher Shannon diversity. Study C: a fixed
    multiplicative bias on three genera.
    """
    x = [a['base'] + rng.gauss(0, 0.45) for a in asvs]
    mean = sum(x) / len(x)
    squeeze = 1.0 - 0.12 * age_idx
    x = [mean + (v - mean) * squeeze for v in x]
    if study == 'STUDY_C':
        x = [v + math.log(STUDY_C_BIAS[a['genus']]) if a['genus'] in STUDY_C_BIAS else v
             for v, a in zip(x, asvs)]
    return x


def resample_person(rng, x, sd=WITHIN_SD):
    """Another sample of the same person: their profile plus within-person noise."""
    return [v + rng.gauss(0, sd) for v in x]


def composition(rng, asvs, age_idx, study=None, bloom=False, x=None):
    """Relative abundances for one sample; Bloom raises Escherichia-Shigella to ~30%."""
    if x is None:
        x = latent(rng, asvs, age_idx, study)
    w = [math.exp(v) for v in x]
    if bloom:
        tot = sum(w)
        ec = [i for i, a in enumerate(asvs) if a['genus'] == 'Escherichia-Shigella']
        for i in ec:                                   # ~30% of the community after renorm
            w[i] += 0.43 * tot / len(ec)
    s = sum(w)
    return [v / s for v in w]


def multinomial(rng, n, p):
    """Plain-python multinomial via sequential binomial-by-normal approximation (no numpy)."""
    counts, left, rem = [], n, 1.0
    for i, pi in enumerate(p):
        if i == len(p) - 1:
            counts.append(left)
            break
        q = 0 if rem <= 0 else min(1.0, pi / rem)
        mu, sd = left * q, math.sqrt(max(left * q * (1 - q), 0))
        c = int(round(min(left, max(0, rng.gauss(mu, sd)))))
        counts.append(c)
        left -= c
        rem -= pi
    return counts


def main():
    rng = random.Random(SEED)
    asvs = make_asvs(rng)

    samples = []        # (sample, role, study, run, age_cat, sex, bmi_cat, country)
    for study, run, with_meta in (('STUDY_A', 'FCA0001', True),
                                  ('STUDY_B', 'FCB0001', True),
                                  ('STUDY_C', 'FCC0001', False)):
        for i in range(12):
            age = AGE_CATS[i % len(AGE_CATS)]
            sex = 'male' if i % 2 else 'female'
            bmi = ['Normal', 'Normal', 'Overweight', 'Obese'][i % 4]
            sid = '%s_%02d' % (study.lower().replace('_', ''), i + 1)
            samples.append((sid, 'panel', study, run,
                            age if with_meta else '', sex if with_meta else '',
                            bmi if with_meta else '', 'Germany' if study != 'STUDY_C' else ''))
    for i, run in enumerate(('FCP0001', 'FCP0002', 'FCP0003')):
        samples.append(('self_P%d' % (i + 1), 'personal', 'self', run,
                        '40s', 'female', 'Overweight', 'Netherlands'))

    # Panel samples are drawn exactly as before, so every signal planted in them is
    # unchanged. Each panel person's latent profile is kept for their repeat sample.
    table = {a['id']: [] for a in asvs}
    latents, person = {}, None
    rng_personal = random.Random(SEED + 3)
    for (sid, role, study, run, age, sex, bmi, country) in samples:
        age_idx = AGE_CATS.index(age) if age else 2     # unknown age: middle band
        if role == 'panel':
            x = latent(rng, asvs, age_idx, study=study)
            latents[sid] = x
        else:
            # The personal samples are ONE person, sampled three times.
            if person is None:
                person = latent(rng, asvs, age_idx)
            x = resample_person(rng_personal, person)
        p = composition(rng, asvs, age_idx, bloom=(sid == 'self_P3'), x=x)
        depth = rng.randint(6000, 20000)
        for a, c in zip(asvs, multinomial(rng, depth, p)):
            table[a['id']].append(c)

    # Repeat samples: a second sample of some panel people, from their own latent profile
    # plus within-person noise. Separate random stream, so nothing above changes.
    rng_rep = random.Random(SEED + 7)
    repeats = []
    info = {s[0]: s for s in samples}
    for subject in REPEATS:
        s = info[subject]
        # The panel sample IS the latent profile, without noise of its own, so the repeat
        # carries the noise of both samples: two personal samples differ by two noise draws.
        p = composition(rng_rep, asvs, 0, x=resample_person(rng_rep, latents[subject],
                                                            sd=WITHIN_SD * math.sqrt(2)))
        for a, c in zip(asvs, multinomial(rng_rep, rng_rep.randint(6000, 20000), p)):
            table[a['id']].append(c)
        repeats.append((subject + '_r', 'repeat', s[2], s[3], s[4], s[5], s[6], s[7], subject))

    columns = [s[0] for s in samples] + [r[0] for r in repeats]
    dada2 = os.path.join(OUT, 'ampliseq_results', 'dada2')
    os.makedirs(dada2, exist_ok=True)

    with open(os.path.join(dada2, 'ASV_table.tsv'), 'w', newline='\n') as fh:
        fh.write('ASV_ID\t' + '\t'.join(columns) + '\n')
        for a in asvs:
            fh.write(a['id'] + '\t' + '\t'.join(str(c) for c in table[a['id']]) + '\n')

    # ampliseq's read tracking, in its own format (thousands separators, a percent column).
    # Retention per step is ordinary for everyone except TECH_OUTLIER, which loses half its
    # reads at quality filtering: a technical outlier for TECHNICAL_QC to find.
    rng_tech = random.Random(SEED + 11)
    nonchim = {c: sum(table[a['id']][i] for a in asvs) for i, c in enumerate(columns)}
    with open(os.path.join(OUT, 'ampliseq_results', 'overall_summary.tsv'), 'w', newline='\n') as fh:
        fh.write('\t'.join(['sample', 'cutadapt_total_processed', 'cutadapt_reverse_complemented',
                            'cutadapt_passing_filters', 'cutadapt_passing_filters_percent', 'DADA2_input',
                            'filtered', 'denoisedF', 'denoisedR', 'denoisedPairs', 'merged', 'nonchim']) + '\n')
        for c in columns:
            chim = rng_tech.uniform(0.93, 0.99)
            merge = rng_tech.uniform(0.90, 0.97)
            denoise = rng_tech.uniform(0.98, 0.995)
            filt = 0.45 if c == TECH_OUTLIER else rng_tech.uniform(0.78, 0.88)
            primer = rng_tech.uniform(0.90, 0.97)
            merged = round(nonchim[c] / chim)
            pairs = round(merged / merge)
            filtered = round(pairs / denoise)
            dada2_in = round(filtered / filt)
            total = round(dada2_in / primer)
            fh.write('\t'.join([c, '{:,}'.format(total), 'NA', '{:,}'.format(dada2_in),
                                '%.1f%%' % (100.0 * dada2_in / total), str(dada2_in), str(filtered),
                                str(round(pairs * 1.004)), str(round(pairs * 1.002)), str(pairs),
                                str(merged), str(nonchim[c])]) + '\n')

    with open(os.path.join(dada2, 'ASV_tax_species.silva.tsv'), 'w', newline='\n') as fh:
        fh.write('ASV_ID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\tconfidence\tsequence\n')
        for a in asvs:
            ph, cl, od, fa = a['ranks']
            fh.write('\t'.join([a['id'], 'Bacteria', ph, cl, od, fa, a['genus'], '',
                                '1.0' if a['genus'] else '0.9', a['seq']]) + '\n')

    header = ['sample', 'role', 'study', 'run', 'age_cat', 'sex', 'bmi_cat', 'country', 'context', 'subject']

    def rows(country_of=None):
        out = []
        for s in samples:
            row = list(s) + ['synthetic timepoint %s' % s[0][-1] if s[1] == 'personal' else '', '']
            if country_of is not None and s[1] == 'panel':
                row[7] = country_of.get(s[2], '')
            out.append(row)
        for r in repeats:
            row = list(r[:8]) + ['', r[8]]
            if country_of is not None:
                row[7] = country_of.get(r[2], '')
            out.append(row)
        return out

    def write(name, sheet_rows):
        with open(os.path.join(OUT, name), 'w', newline='\n') as fh:
            w = csv.writer(fh, lineterminator='\n')
            w.writerow(header)
            w.writerows(sheet_rows)

    write('samplesheet.csv', rows())

    # The same sheet with country collinear with study (A = Germany, B = Austria), the shape
    # of a real panel where every study comes from one country. Tests that the batch models
    # report such a variable as confounded instead of failing on an undefined test.
    write('samplesheet_collinear.csv', rows({'STUDY_A': 'Germany', 'STUDY_B': 'Austria'}))

    # Without repeat samples there is no within-person reference: the identity check falls
    # back to its between-person half and the change statistics are not assessable.
    write('samplesheet_norepeats.csv', [r for r in rows() if r[1] != 'repeat'])

    # A lab mix-up: the third "personal" sample is in fact someone else (panel person
    # studya_07, who has no repeat), which the identity check must point out. studya_07 is an
    # ordinary stranger; with only 22 genera the within- and between-person ranges overlap more
    # than in real data, and the few strangers who happen to be close would pass as one person.
    mixup = []
    for r in rows():
        if r[0] == 'self_P3':
            continue
        if r[0] == 'studya_07':
            r = ['studya_07', 'personal', 'self', 'FCP0003', '40s', 'female', 'Overweight', 'Netherlands',
                 'synthetic timepoint 3', '']
        mixup.append(r)
    write('samplesheet_mixup.csv', mixup)

    print('wrote %d ASVs x %d samples (%d of them repeats) to %s'
          % (len(asvs), len(columns), len(repeats), OUT))


if __name__ == '__main__':
    main()
