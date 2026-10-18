# Biocentric/gutportrait: Output

## Introduction

This document describes the output produced by the pipeline. All paths are relative to the
top-level results directory (`--outdir`). File names start with the comparison arm
(`--arm`: `native` or `v4projected`), so the outputs of both arms can share one results
directory.

> [!IMPORTANT]
> Read [the gates](#gates) first. They say which numbers may be read as biology at all. A
> position that is not stable under batch correction, or that sits within the spread
> between reference studies, is a position, not a finding.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Genus tables](#genus-tables) - ASV counts collapsed to genus, with unassigned reads accounted for
- [Taxon percentiles](#taxon-percentiles) - where each personal sample sits in the panel, per genus, in CLR space
- [Batch diagnosis](#batch-diagnosis) - how much of the panel's variation is study, and how much each metadata variable
- [Batch correction](#batch-correction) - the panel's studies harmonised with MMUPHin; personal samples never corrected
- [Gates](#gates) - the injected-signal positive control and the batch-dominance rules
- [Percentile stability](#percentile-stability) - which genus positions survive batch correction
- [Alpha diversity](#alpha-diversity) - diversity at one fixed depth, and where each personal sample sits
- [Deviation flags](#deviation-flags) - extremes-only checks for each personal sample, each with how often it fires for an ordinary person
- [Within-person change and identity](#within-person-change-and-identity) - how much your samples differ, against how much one person's samples ordinarily differ, and whether they are one person
- [Technical comparability](#technical-comparability) - whether your samples were processed like the panel's, step by step
- [MultiQC](#multiqc) - software versions, parameters and the methods text
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### Genus tables

<details markdown="1">
<summary>Output files</summary>

- `genus/`
  - `<arm>.genus_counts.tsv`: genus × sample read counts, for every sample in the ASV table.
  - `<arm>.collapse_stats.tsv`: one row per sample, with `total_reads`, `unassigned_reads`,
    `unassigned_fraction`, and `n_taxa` (genus rows with at least one read, `Unassigned`
    included).

</details>

ASVs are summed by the `Genus` column of ampliseq's named-rank taxonomy. ASVs without a
genus call are not dropped: they become one `Unassigned` row. Dropping them would quietly
raise the share of every other genus. A high `unassigned_fraction` therefore shows up as a
large `Unassigned` part, instead of disappearing.

### Taxon percentiles

<details markdown="1">
<summary>Output files</summary>

- `percentiles/`
  - `<arm>.taxon_percentiles.tsv`: one row per genus per personal sample.
  - `<arm>.percentile_params.tsv`: the settings and counts the percentiles were computed
    with: `pseudocount`, `min_prevalence`, `drop_unassigned`, `n_genera_tested`,
    `n_genera_panel_total`, `panel_samples` and `personal_samples`.

</details>

The calculation:

1. Genera are kept if they are present in at least `--min_prevalence` of the **panel**
   samples. Personal samples play no part in this choice.
2. Personal and panel samples are aligned to that same set of genera, and every count gets
   `--pseudocount` added.
3. Each sample is transformed by the centred log-ratio (CLR): the log of each part minus
   the mean log of all parts of that sample.
4. For each genus, the percentile is the share of panel samples whose CLR is **at or
   below** the personal sample's.

| Column             | Meaning                                                                                                   |
| ------------------ | --------------------------------------------------------------------------------------------------------- |
| `genus`            | genus, or `Unassigned`                                                                                    |
| `sample`           | the personal sample                                                                                       |
| `context`          | the sample's `context` from the samplesheet, if given                                                     |
| `clr`              | its centred-log-ratio value for this genus                                                                |
| `percentile`       | the share of panel samples (%) at or below that value                                                     |
| `panel_median_clr` | the panel median, for scale                                                                               |
| `panel_n`          | number of panel samples compared against                                                                  |
| `panel_prevalence` | fraction of panel samples in which the genus is present                                                   |
| `personal_counts`  | raw reads for this genus in the personal sample, so a percentile resting on a handful of reads is visible |

Rows are sorted by sample and then by percentile, highest first.

How to read it:

- CLR rather than relative abundance, because a relative abundance is a ratio to an
  unmeasured total: one genus's share can move only because another genus changed. In
  CLR, a value moves only when that genus moves relative to the rest of the sample.
- A percentile's resolution is limited by `panel_n`. Against 36 panel samples, each sample
  is a step of about 2.8 points, and 100 means "above every panel sample", not "extreme
  on some absolute scale".
- With dozens of genera per sample, some sit at the edges by chance alone. A single high or
  low percentile is not a deviation. The [deviation flags](#deviation-flags) apply stricter,
  pre-registered criteria.

### Percentile stability

<details markdown="1">
<summary>Output files</summary>

- `percentiles/`
  - `<arm>.corrected.taxon_percentiles.tsv`: the same percentiles against the
    batch-corrected panel.
  - `<arm>.percentile_stability.tsv`: per genus and personal sample, both percentiles,
    their categories (`low` ≤ `--stability_low_percentile`, `high` ≥
    `--stability_high_percentile`, `middle`), whether they agree (`stable`) and the shift.

</details>

A position is stable when its category is the same against the uncorrected and the
corrected panel. An extreme that appears in only one of the two was created, or hidden, by
differences between studies, and is not a finding.

### Batch diagnosis

<details markdown="1">
<summary>Output files</summary>

- `batch/`
  - `<arm>.permanova.tsv`: one row per model term: `model`, `term`, `df`, `R2`, `F`, `p`,
    `n`, and a `note` when a model was not testable or is confounded.
  - `<arm>.batch_summary.tsv`: `r2_study`, whether it exceeds `--banner_r2_study`
    (`study_banner`), and every setting used.

</details>

PERMANOVA (vegan `adonis2`) on Aitchison distance, the Euclidean distance between CLR
profiles, computed on **panel** samples only with the same genus filter and pseudocount
as the percentiles. Models:

- `study`: R2(study), how much of the panel's variation is between studies. Above
  `--banner_r2_study` (default 0.10), every population-facing claim is orientation only.
- `study + <variable>` for each metadata column with at least two levels: marginal R2 for
  both terms. A variable explains more than study only if its R2 is the larger one, which
  is the batch-dominance rule of DECISIONS.md D12. If every sample with that variable comes
  from one study, the two cannot be separated; the model is fitted without study and the
  `note` says `confounded`.

With a single panel study, R2(study) cannot be estimated and is reported as missing.

### Batch correction

<details markdown="1">
<summary>Output files</summary>

- `batch/`
  - `<arm>.genus_counts_corrected.tsv`: the genus table with the **panel** studies
    harmonised by MMUPHin (`adjust_batch`, batch = study) and the personal samples
    unchanged. Same layout as `genus/<arm>.genus_counts.tsv`.
  - `<arm>.correction_summary.tsv`: method, studies, samples, the covariates protected
    (`--correction_covariates`, only those filled for every panel sample) and a note.

</details>

The personal samples form their own batch, and three samples cannot separate their lab
from their biology. Correcting that batch would pull them toward the panel and erase the
deviation the pipeline measures, so they are never corrected. With a single panel study
there is nothing to harmonise, and the table is passed through with a note.

### Gates

<details markdown="1">
<summary>Output files</summary>

- `gates/`
  - `<arm>.positive_control.tsv`: verdict (`PASS`, `FAIL` or `NOT_APPLICABLE`), the planted
    genus, `n_spiked`, the expected CLR shift, the median surviving shift, the `retained`
    share and the paired Wilcoxon `p`.
  - `<arm>.gates.tsv`: one row per gate, with `status`, `value`, `threshold` and `detail`.

</details>

**Positive control, injected signal.** A common, mid-abundance genus is multiplied by
`--pc_fold` in a random half of every panel study, and the same correction is applied to
the spiked and the unspiked table. For each spiked sample, the shift that survives is its
CLR in the spiked, corrected table minus its CLR in the unspiked, corrected one. The control
passes when the median surviving shift has the planted sign, a paired Wilcoxon test gives
p < `--pc_alpha`, and at least `--pc_retain` of the planted shift survives. **If it fails,
the run stops**: a correction that erases an effect planted on purpose would erase real ones.

**Batch dominance** (`gates.tsv`):

| Gate                            | Meaning                                                                                                                                                                               |
| ------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `study_banner`                  | `BANNER` when R²(study) exceeds `--banner_r2_study`: every population-facing claim is orientation only                                                                                |
| `contrast:<variable>`           | `finding` only if the variable's R² exceeds R²(study) in the same model and p < `--dominance_alpha`; otherwise `orientation_only` or `not_testable`                                   |
| `personal_positioning:<sample>` | `interpretable` only if the sample lies further from the pooled panel centroid than q95 of within-study distances plus the largest study shift (see below); otherwise `within_spread` |
| `personal_positioning`          | the summary over the personal samples: `interpretable`, `within_spread` or `mixed`                                                                                                    |

Why that bound: by the triangle inequality, an ordinary person (within the 95th percentile of
their own study's spread) in a study whose centroid is shifted by at most the largest shift
seen in the panel lies at most q95 + that shift from the pooled centroid. A personal sample
beyond the bound cannot be explained as an ordinary person measured in a shifted lab or
cohort; one within it can, so its deviation is not read as biology. The gate needs at least
two panel studies with `--min_study_n` samples each; with fewer, it is `not_assessable`.

### Alpha diversity

<details markdown="1">
<summary>Output files</summary>

- `alpha/`
  - `<arm>.alpha_diversity.tsv`: per sample: `reads`, `rarefied`, `shannon`, `observed`
    (ASVs) and `pielou`. Samples below `--rarefaction_depth` keep a row with empty metrics.
  - `<arm>.alpha_position.tsv`: per personal sample and metric: `panel_percentile`,
    `panel_median`, the range of study medians, and whether the value lies below the low
    (or above the high) percentile of **every** panel study.
  - `<arm>.alpha_params.tsv`: depth, iterations, seed and sample counts.

</details>

Every sample is subsampled to the same depth (`--rarefaction_depth`, default 5,000 reads),
`--rarefaction_iterations` times with fixed seeds, and the metrics are averaged. Without a
common depth, observed richness would rank samples by how deeply they were sequenced.
Panel studies with fewer than `--min_study_n` samples are left out of the per-study bounds.

### Deviation flags

<details markdown="1">
<summary>Output files</summary>

- `flags/`
  - `<arm>.deviation_flags.tsv`: every check for every personal sample, raised or not, with
    `value`, `threshold`, `detail`, `reference_fire_rate`, `arm` and `agreement`.
  - `<arm>.flag_calibration.tsv`: per check, how often it fires for panel people
    (`reference_n`, `reference_fired`, `reference_fire_rate`), the `nominal_rate` it was
    designed for, and a note.

</details>

| Flag                         | Raised when                                                                                                            |
| ---------------------------- | ---------------------------------------------------------------------------------------------------------------------- |
| `low_alpha`                  | Shannon below the low percentile of the pooled panel **and** of every panel study                                      |
| `single_taxon_dominance`     | one assigned genus above `--dominance_threshold` (default 50%) of the sample                                           |
| `enterobacteriaceae_bloom`   | family _Enterobacteriaceae_ above the high percentile of the pooled panel and of every study                           |
| `proteobacteria_bloom`       | phylum _Proteobacteria_ (SILVA) or _Pseudomonadota_ (GTDB), same rule                                                  |
| `faecalibacterium_depletion` | genus _Faecalibacterium_ below the low percentile, same rule                                                           |
| `core_genus_missing`         | a genus in ≥ `--core_prevalence` of the panel has no reads, although the panel median predicts ≥ `--core_min_expected` |
| `beta_outlier`               | Aitchison distance to the panel centroid above the high percentile of the panel's own (leave-one-out) distances        |
| `excess_extremes`            | more genera in the outer percentile tails than the high percentile of what panel people have                           |

**Calibration.** How often would a check fire for an ordinary person? Every panel sample is
treated as if it were a personal sample and checked against the panel without it, with the
same code. `reference_fire_rate` is the share of panel people for whom the check fires. A
check that fires for one ordinary person in five is weak evidence on its own; one that
fires for nobody is strong. Missing core genera, for example, are common: with dozens of
genera each present in 95% of people, many ordinary people lack one or two. The
leave-one-out counts of extreme genera are also the reference for `excess_extremes`:
everybody has some genera in the outer tails, roughly one in ten by chance alone, so only
an excess over what panel people show is a flag.

How to read it:

- Every check is reported, so "nothing raised" is a visible result, and the strongest
  negative result 16S can give. It is not a certificate of health.
- Read a raised flag next to its `reference_fire_rate`.
- Abundance checks use relative abundance, so a flag reads as "31% against a panel 95th
  percentile of 2.6%". The distance check uses CLR.
- Requiring a value to be extreme against **every** study, not just the pooled panel, is
  the between-study-spread rule: a value is unusual only if no single reference study would
  call it ordinary.
- A flag from one arm is provisional. The design requires the native and V4-projected arms
  to agree before a flag reaches the report, and the `agreement` column says which.
- Read each flag next to its sample's `context`. A flag that appears only in a sample taken
  just after an unusual event, such as travel or an illness, is a candidate for that event, not evidence
  about the usual state.

### Within-person change and identity

<details markdown="1">
<summary>Output files</summary>

- `within_person/`
  - `<arm>.within_person_reference.tsv`: one row per reference pair (two samples of one
    panel person, from the `repeat` rows): `aitchison`, `trimmed_distance`, `delta_shannon`,
    `same_run`, and `genera_changed_loo`, the number of genera that changed more than usual,
    each pair measured against the others.
  - `<arm>.within_person_summary.tsv`: the within-person and between-person distributions of
    both distances, pooled and per study.
  - `<arm>.personal_changes.tsv`: one row per pair of your samples: the Aitchison distance
    with its `within_person_percentile` and `between_person_percentile`, the change in
    Shannon diversity (`delta_shannon`, with a z-score and a percentile against the
    reference pairs), `genera_changed` against `genera_changed_expected`, and a `caveat`.
  - `<arm>.personal_genus_changes.tsv`: per pair and genus: `z`, the change in CLR abundance
    in units of that genus's ordinary within-person change (`ordinary_change_sd`); `changed`
    when |z| exceeds `--change_z`; the direction and both relative abundances.
  - `<arm>.identity.tsv`: per pair of your samples, a `verdict` with both distances and their
    percentiles. The last row (`all`) is the overall verdict.

</details>

Two samples of one person are never identical. Whether a difference between your samples
means anything depends on how much one person's samples ordinarily differ, and three
samples of one person cannot measure that. The reference is borrowed from panel studies
that sampled the same people twice (`repeat` rows, see
[repeat samples](usage.md#repeat-samples)). Both references, two samples of one person and
two different people, come from within one study, processed in one lab. Your samples are
compared only with each other.

How to read `personal_changes.tsv`:

- `within_person_percentile` 80 means your two samples differ more than 80% of the
  reference pairs do: an ordinary change. Near 100 is a change larger than nearly anyone's.
- `genera_changed` counts genera that moved more than `--change_z` times their ordinary
  change. Some always do, by chance: `genera_changed_expected` is the median for a
  reference person, and `genera_changed_percentile` places your count among theirs.
- The `caveat` notes when your two samples were sequenced in different runs while most
  reference pairs were not: part of the difference may then be technical.
- The reference pairs were taken at the intervals their studies chose. Over a longer
  interval, larger changes are expected.

Identity, from `identity.tsv`. You sent in your own samples, so one person is the starting
assumption, and the check looks for evidence against it, as a guard against a mix-up:

| Verdict                     | Meaning                                                                                                                                                                                                     |
| --------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `same_person`               | the Aitchison distance is ordinary for one person: at or below the `--identity_within_percentile` of the within-person pairs, or below the `--identity_between_percentile` of the pairs of different people |
| `same_person_large_change`  | it is not, but the trimmed distance is: the difference sits in a few genera, as after a bloom or a crash of one group                                                                                       |
| `possibly_different_people` | neither is: most genera differ, as between two people. Check for a sample mix-up                                                                                                                            |
| `unclear`                   | without a within-person reference, a pair that is not closer than nearly all pairs of different people                                                                                                      |

The trimmed distance is the Aitchison distance after centring the CLR differences on their
median and leaving out the 5% of genera that changed most, so that one genus cannot make a
sample look like someone else's. The overall row reads `one_person`, `one_sample_differs`
(naming the sample), `possibly_different_people` or `unclear`.

### Technical comparability

<details markdown="1">
<summary>Output files</summary>

- `qc/`
  - `<arm>.technical_metrics.tsv`: per sample: `input_reads`, `primer_retention`,
    `filter_retention`, `denoise_retention`, `merge_rate`, `chimera_retention`,
    `end_to_end_retention`, `analysed_reads` and `unassigned_fraction`, from ampliseq's
    `overall_summary.tsv` and the genus collapse.
  - `<arm>.technical_qc.tsv`: per personal sample and metric: `panel_percentile`,
    `panel_median`, `modified_z` and `atypical`; one `overall` row per sample.

</details>

A profile can look unusual for technical reasons: a sample that lost most of its reads at
quality filtering, merged poorly or came out mostly chimeric has a different composition for
reasons that have nothing to do with the gut. Each step of ampliseq's read tracking is turned
into a retention rate and compared with the panel. A metric is `atypical` when, on its worse
side, its modified z-score (Iglewicz and Hoaglin) is beyond `--qc_max_z` against the pooled
panel and against every panel study. Rates within 5 percentage points of the panel median
are never atypical, and read counts are compared on a log scale.

Atypical is not an error. It means the sample was processed differently, and its other
results should be read with that in mind. When all your samples are atypical in the same
metrics, your lab's protocol differed from the panel's (read length, amplicon, chemistry),
which is worth knowing in itself.

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. For gutportrait it currently carries the software versions, the parameters of the run and a methods description with citations. The dedicated gutportrait report (HTML and PDF) is planned. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://docs.seqera.io/platform-cloud/reports/overview) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
