# Biocentric/gutportrait: Output

## Introduction

This document describes the output produced by the pipeline. All paths are relative to the
top-level results directory (`--outdir`). File names start with the comparison arm
(`--arm`: `native` or `v4projected`), so the outputs of both arms can share one results
directory.

> [!IMPORTANT]
> The integration stage and its gates are not built yet. Until they are, every number
> below is computed against the **whole, uncorrected** panel. Read it as a position, not
> a finding.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Genus tables](#genus-tables) - ASV counts collapsed to genus, with unassigned reads accounted for
- [Taxon percentiles](#taxon-percentiles) - where each personal sample sits in the panel, per genus, in CLR space
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

| Column             | Meaning                                                                                           |
| ------------------ | ------------------------------------------------------------------------------------------------- |
| `genus`            | genus, or `Unassigned`                                                                            |
| `sample`           | the personal sample                                                                               |
| `clr`              | its centred-log-ratio value for this genus                                                        |
| `percentile`       | the share of panel samples (%) at or below that value                                             |
| `panel_median_clr` | the panel median, for scale                                                                       |
| `panel_n`          | number of panel samples compared against                                                          |
| `panel_prevalence` | fraction of panel samples in which the genus is present                                           |
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
  low percentile is not a deviation. The deviation flags (planned) apply stricter,
  pre-registered criteria.

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
