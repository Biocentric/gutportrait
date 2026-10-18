# Biocentric/gutportrait: Output

## Introduction

This document describes the output produced by the pipeline. All paths are relative to the
top-level results directory (`--outdir`). File names start with the comparison arm
(`--arm`: `native` or `v4projected`), so the outputs of both arms can share one results
directory.

> [!IMPORTANT]
> Batch correction and the two gates are not built yet. The batch diagnosis below measures
> how much of the panel's variation is study, but nothing is corrected for it. Until the
> gates exist, every number is computed against the **whole, uncorrected** panel. Read it
> as a position, not a finding.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Genus tables](#genus-tables) - ASV counts collapsed to genus, with unassigned reads accounted for
- [Taxon percentiles](#taxon-percentiles) - where each personal sample sits in the panel, per genus, in CLR space
- [Batch diagnosis](#batch-diagnosis) - how much of the panel's variation is study, and how much each metadata variable
- [Alpha diversity](#alpha-diversity) - diversity at one fixed depth, and where each personal sample sits
- [Deviation flags](#deviation-flags) - extremes-only checks for each personal sample
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
    `value`, `threshold`, `detail`, `arm` and `agreement`.

</details>

| Flag                         | Raised when                                                                                                            |
| ---------------------------- | ---------------------------------------------------------------------------------------------------------------------- |
| `low_alpha`                  | Shannon below the low percentile of the pooled panel **and** of every panel study                                      |
| `single_taxon_dominance`     | one assigned genus above `--dominance_threshold` (default 50%) of the sample                                           |
| `enterobacteriaceae_bloom`   | family _Enterobacteriaceae_ above the high percentile of the pooled panel and of every study                           |
| `proteobacteria_bloom`       | phylum _Proteobacteria_ (SILVA) or _Pseudomonadota_ (GTDB), same rule                                                  |
| `faecalibacterium_depletion` | genus _Faecalibacterium_ below the low percentile, same rule                                                           |
| `core_genus_missing`         | a genus in ≥ `--core_prevalence` of the panel has no reads, although the panel median predicts ≥ `--core_min_expected` |
| `beta_outlier`               | Aitchison distance to the panel centroid above the high percentile of the panel's own distances                        |

How to read it:

- Every check is reported, so "nothing raised" is a visible result, and the strongest
  negative result 16S can give. It is not a certificate of health.
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
