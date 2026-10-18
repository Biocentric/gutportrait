# Biocentric/gutportrait: Usage

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

gutportrait is the second of two pipeline runs. nf-core/ampliseq turns reads into ASVs;
gutportrait reads ampliseq's results directory and positions your samples against a
reference panel that was processed **in the same ampliseq run**. It never touches a FASTQ.

A full analysis is therefore:

1. one nf-core/ampliseq run over your samples and the panel together (the
   [README](../README.md#step-1--profile-with-ampliseq) lists the pinned settings), then
2. one gutportrait run over that results directory, with a samplesheet that says which
   samples are yours and which are the panel.

## Samplesheet input

```bash
--input '[path to samplesheet file]'
```

A comma-separated file with a header row and one row per sample of the ampliseq run that
should take part in the comparison:

```csv title="samplesheet.csv"
sample,role,study,run,age_cat,sex,bmi_cat,country,context,subject
self_sample1,personal,self,RUN_A,40s,female,Normal,Netherlands,baseline,
self_sample2,personal,self,RUN_B,40s,female,Normal,Netherlands,after a diet change,
panel_ERR0000001,panel,PRJEB00000,RUN_C,,,,,,
panel_ERR0000002,repeat,PRJEB00000,RUN_C,,,,,,panel_ERR0000001
cohort_S0001,panel,COHORT_X,RUN_D,30s,male,Overweight,Germany,,
```

| Column    | Required | Description                                                                                                                                                                                                                   |
| --------- | -------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`  | yes      | Sample name exactly as it appears as a column of ampliseq's `dada2/ASV_table.tsv`. Must match `^[a-zA-Z][a-zA-Z0-9_]+$`, the same rule ampliseq applies.                                                                      |
| `role`    | yes      | `personal` for your samples, `panel` for reference samples, `repeat` for a second sample of a panel person (see [repeat samples](#repeat-samples)).                                                                           |
| `study`   | yes      | The batch unit, usually the ENA project. Batch diagnosis and correction work per study. Personal samples need a study of their own: a study that contains both roles is rejected, so your samples can never pass for a batch. |
| `run`     | yes      | The sequencing run (flowcell), the same value as the `run` column of the ampliseq samplesheet. It is the unit DADA2 learned its error model on.                                                                               |
| `age_cat` | no       | `20s`, `30s`, `40s`, `50s`, `60s` or `70+`. For panel rows, a set such as `30s                                                                                                                                                | 40s | 50s` records a study's inclusion criteria when the person's own band is not published (see [cohort matching](#cohort-matching)). |
| `sex`     | no       | `male`, `female` or `other`; for panel rows also a set, as for `age_cat`.                                                                                                                                                     |
| `bmi_cat` | no       | `Underweight`, `Normal`, `Overweight` or `Obese`; for panel rows also a set, as for `age_cat`.                                                                                                                                |
| `country` | no       | Free text.                                                                                                                                                                                                                    |
| `context` | no       | Free text, for personal samples: when or why the sample was taken, e.g. `baseline` or `after a diet change`. It is carried into the percentile, alpha and flag tables, so every number can be read next to it.                |
| `subject` | no       | For `repeat` rows, required: the panel person the sample belongs to, given as that panel row's `subject`, or its `sample` name when `subject` is empty. Your own samples are always taken to be one person.                   |

Leave optional cells empty where a value is unknown. For public ENA panels that is the
normal case: they publish study-level inclusion criteria, not per-sample metadata. The
metadata columns are for cohort matching and the real-effect positive control. Samples
without metadata are still used, but only as part of the panel as a whole.

The run stops before any work starts when:

- a sample name appears twice, or does not match the pattern above;
- there is no personal sample, or there are fewer panel samples than `--min_panel_samples`;
- a study contains both personal and panel samples;
- a `repeat` row's `subject` names no panel person of the same study;
- a personal sample gives a set of bands (with `|`) instead of single values.

It stops at the first module when a samplesheet sample has no column in the ASV table.
Samples that are in the ASV table but not in the samplesheet are ignored, with a note in
the task log; that is how you analyse a subset of one large ampliseq run.

An [example samplesheet](../assets/samplesheet.csv) is provided with the pipeline.

### Cohort matching

The comparisons use the whole panel by default. `MATCH_COHORT` adds a comparison with the
panel people who share your bands, the `age_cat`, `sex` and `bmi_cat` of your samples. The
match starts exact and is widened one step at a time, in a fixed order, until at least
`--min_matched` people (default 50) qualify:

1. same age band, sex and BMI band;
2. age within one band either side, same sex and BMI band;
3. age within one band either side, same sex;
4. age within one band either side.

If no step reaches the minimum, there is no matched comparison, and the report says how far
the panel's metadata went. A panel person counts only if their band is **known** to satisfy
the step: an empty cell never matches. Where a study publishes only its inclusion criteria,
give the possible bands as a set, for example `30s|40s|50s` for "adults 30 to 50"; such a
person matches only when every band in the set satisfies the step, so a match is never a
guess. Sets are treated as unknown by the batch models, which need one value per person.

Most public panels publish no per-person age or sex, so on an ENA-only panel a matched
comparison is often impossible. Studies that publish them (American Gut, or a cohort shared
on request) are what make it work.

### Repeat samples

Some panel studies sampled the same people twice, for example before and after an
intervention. Their second samples go in as `repeat` rows, with `subject` naming the
person. They never enter the reference distributions (percentiles, batch statistics,
flags); they serve one purpose: measuring how much two samples of one person ordinarily
differ. That is the noise floor against which the differences between your own samples
are read, and the reference for the identity check (see
[within-person change and identity](output.md#within-person-change-and-identity)).

- At least `--min_reference_pairs` pairs (default 10) are needed. Without them, the change
  statistics between your samples are not assessable, and the identity check uses only its
  comparison with pairs of different people.
- Use a person's own second sample, from the same study and sequencing run where the study
  allows it, so that the pair differs the way two of your samples differ.
- The interval matters: two samples taken days apart differ less than two taken months
  apart. Prefer studies whose interval resembles yours, and read the result with the
  interval in mind.

## The ampliseq results directory

```bash
--ampliseq_results '[path to the --outdir of an nf-core/ampliseq run]'
--taxonomy_db silva
```

gutportrait reads two files from it:

- `dada2/ASV_table.tsv`, the ASV counts for every sample of the run;
- `dada2/ASV_tax_species.<taxonomy_db>.tsv`, the taxonomy, with named rank columns.

`--taxonomy_db` must name the database ampliseq classified against. The default, `silva`,
matches ampliseq's `--dada_ref_taxonomy silva`. The rank-numbered
`qiime2/rel_abundance_tables/rel-table-N.tsv` files are not used: their numbering is
positional and shifts with the database, so under SBDI-GTDB `rel-table-6` is Family.

## Comparison arms

```bash
--arm native        # or: v4projected
```

Run gutportrait once per arm. The arm names the outputs, so both arms can share one
results directory.

- `native`: a V3–V4 comparison. Your reads and a panel amplified with the same primers,
  processed at full length.
- `v4projected`: a V4-only comparison against the much larger body of V4 data. It comes
  from a separate, single-end ampliseq run in which every forward read is cut at the 515F
  primer site and truncated to a fixed length, so V3–V4 and V4 datasets share one window.

## Analysis options

| Parameter                  | Default    | Meaning                                                                                                                                                                  |
| -------------------------- | ---------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `--min_panel_samples`      | `20`       | Refuse to run with fewer panel samples than this.                                                                                                                        |
| `--min_prevalence`         | `0.1`      | Keep genera present in at least this fraction of **panel** samples. The filter never looks at personal samples: selecting features for being present in you is circular. |
| `--pseudocount`            | `0.5`      | Added to every count before the centred-log-ratio transform.                                                                                                             |
| `--min_study_n`            | `10`       | Panel studies smaller than this are left out of per-study bounds and PERMANOVA models.                                                                                   |
| `--rarefaction_depth`      | `5000`     | Reads per sample for alpha diversity; samples below it are listed and get no alpha diversity.                                                                            |
| `--rarefaction_iterations` | `20`       | Rarefaction repeats averaged per sample.                                                                                                                                 |
| `--permutations`           | `999`      | Permutations for the PERMANOVA tests.                                                                                                                                    |
| `--banner_r2_study`        | `0.1`      | If study explains more of the panel's variance than this, every population claim is orientation only.                                                                    |
| `--min_matched`            | `50`       | The smallest matched group a comparison with people like you is reported for (see [cohort matching](#cohort-matching)).                                                  |
| `--seed`                   | `20261018` | Seed for rarefaction and permutations, so every run is reproducible.                                                                                                     |

Deviation flags (DECISIONS.md D12):

| Parameter                | Default | Meaning                                                                                                           |
| ------------------------ | ------- | ----------------------------------------------------------------------------------------------------------------- |
| `--flag_low_percentile`  | `5`     | A low value must fall below this percentile of the pooled panel and of every panel study.                         |
| `--flag_high_percentile` | `95`    | A high value must exceed this percentile of the pooled panel and of every panel study.                            |
| `--dominance_threshold`  | `0.5`   | Flag one genus above this fraction of the sample.                                                                 |
| `--core_prevalence`      | `0.95`  | Genera present in at least this fraction of the panel count as core.                                              |
| `--core_min_expected`    | `10`    | A missing core genus is flagged only if the panel median predicts at least this many reads at the sample's depth. |

Batch correction and gates (DECISIONS.md D12):

| Parameter                     | Default | Meaning                                                                                                        |
| ----------------------------- | ------- | -------------------------------------------------------------------------------------------------------------- |
| `--correction_covariates`     | none    | Samplesheet columns to protect during batch correction, e.g. `age_cat,sex`. Incomplete columns are dropped.    |
| `--pc_fold`                   | `2`     | Positive control: fold change planted in half of each panel study.                                             |
| `--pc_alpha`                  | `0.01`  | Positive control: the planted effect must survive correction at this significance.                             |
| `--pc_retain`                 | `0.5`   | Positive control: share of the planted shift that must survive; below it the run stops.                        |
| `--pc_min_prevalence`         | `0.9`   | Positive control: the planted genus is chosen among genera present in at least this share of the panel.        |
| `--dominance_alpha`           | `0.05`  | A contrast is a finding only if its R² exceeds R²(study) and p is below this.                                  |
| `--stability_low_percentile`  | `10`    | Genus positions at or below this percentile count as low.                                                      |
| `--stability_high_percentile` | `90`    | Genus positions at or above this percentile count as high.                                                     |
| `--qc_max_z`                  | `3.5`   | Technical comparability: a metric is atypical beyond this modified z-score, against the panel and every study. |

Within-person change and identity:

| Parameter                       | Default | Meaning                                                                                                                |
| ------------------------------- | ------- | ---------------------------------------------------------------------------------------------------------------------- |
| `--change_z`                    | `2`     | A genus changed between two of your samples when its change exceeds this many times its ordinary within-person change. |
| `--min_reference_pairs`         | `10`    | Fewer `repeat` pairs than this: no within-person reference is used.                                                    |
| `--identity_within_percentile`  | `95`    | Identity: a distance at or below this percentile of the within-person pairs is ordinary for one person.                |
| `--identity_between_percentile` | `5`     | Identity: a distance below this percentile of the pairs of different people is ordinary for one person.                |

The defaults are deliberately strict. The values used are written next to the results, in
the `*_params.tsv`, `*.batch_summary.tsv`, `*.correction_summary.tsv` and
`*.positive_control.tsv` files.

## Test profile

```bash
nextflow run Biocentric/gutportrait -profile test,docker --outdir results_test
```

`-profile test` needs no input of your own. It runs on a synthetic dataset in
`tests/fixtures/synthetic/`: 64 ASVs over 36 panel samples in three studies, 16 repeat
samples of panel people, and three personal samples of one invented person, laid out
exactly like an ampliseq 2.18.0 results directory, read tracking included. Known signals
are planted in it, so that every stage has an answer to recover:

| Planted signal                                                                          | Recovered by                                                                                                                                                                                                                                                                                                                                           |
| --------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| Personal sample `self_P3` carries an _Escherichia-Shigella_ bloom (about 31% of reads)  | `TAXON_PERCENTILES` (100th percentile) and `DEVIATION_FLAGS` (Enterobacteriaceae and Proteobacteria blooms, beta outlier and an excess of extreme genera, all raised; nothing raised for `self_P1`/`self_P2`)                                                                                                                                          |
| The three personal samples are one person, sampled three times                          | `WITHIN_PERSON`: identity `one_person`; the two pairs with `self_P3` are `same_person_large_change`, because the difference sits in one genus                                                                                                                                                                                                          |
| The bloom is a change within that person                                                | `WITHIN_PERSON`: _Escherichia-Shigella_ changes by z = 17 between `self_P1` and `self_P3`; 5 genera changed against 1 expected by chance, and 1 between `self_P1` and `self_P2`                                                                                                                                                                        |
| Personal sample `self_P1` lost half its reads at quality filtering                      | `TECHNICAL_QC`: filter and end-to-end retention atypical for `self_P1`; nothing atypical for `self_P2`/`self_P3`                                                                                                                                                                                                                                       |
| `samplesheet_mixup.csv`: the third "personal" sample is someone else                    | `WITHIN_PERSON`: identity `one_sample_differs`, naming that sample                                                                                                                                                                                                                                                                                     |
| `samplesheet_norepeats.csv`: no repeat samples                                          | `WITHIN_PERSON`: change statistics not assessable, identity from the between-person comparison only                                                                                                                                                                                                                                                    |
| Every result above                                                                      | `REPORT`: one self-contained HTML page and an A4 PDF that explain each of them in plain language, with the bloom, the technical outlier and the identity verdict in the "At a glance" table                                                                                                                                                            |
| Older age bands have a more even composition, so Shannon diversity rises with `age_cat` | `ALPHA_DIVERSITY` (mean Shannon 3.68 in the 20s rising to 4.04 in the 60s) and `BATCH_DIAGNOSE` (age R² 0.37, p = 0.001); the real-effect positive control is planned                                                                                                                                                                                  |
| Study C carries _Bacteroides_ ×3, _Bifidobacterium_ ×2 and _Segatella_ ×0.3             | `BATCH_DIAGNOSE`: R²(study) 0.34, p = 0.001, so the orientation-only banner is raised; between the unbiased studies A and B, R² 0.02. `BATCH_CORRECT` removes it (study medians of _Bacteroides_ 0.14 / 0.12 / 0.30 before, 0.17 / 0.15 / 0.16 after)                                                                                                  |
| A genus doubled on purpose in half of every study (made by `POSITIVE_CONTROL` itself)   | `POSITIVE_CONTROL`: 94% of the planted shift survives correction, p = 7.6 × 10⁻⁶, so the gate passes                                                                                                                                                                                                                                                   |
| Study C has no per-sample metadata, like the ENA panels                                 | `BATCH_DIAGNOSE` tests the metadata variables on studies A and B only. `MATCH_COHORT` never matches anyone from study C; with the default `--min_matched 50` no step of the ladder qualifies (2, 4, 6, 14 people), so there is no matched comparison; with a minimum of 10 the match widens to age bands only (14 people, 7 each from studies A and B) |
| Personal samples carry a `context` label                                                | every personal-sample table carries it through                                                                                                                                                                                                                                                                                                         |

The fixture is invented, including the personal samples and their metadata.
`python tests/fixtures/make_synthetic_fixture.py` regenerates it byte-identically.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run Biocentric/gutportrait \
    --input ./samplesheet.csv \
    --ampliseq_results ./ampliseq_results \
    --outdir ./results \
    -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/running/run-pipelines#configuring-pipelines), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run Biocentric/gutportrait -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
ampliseq_results: './ampliseq_results'
outdir: './results/'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull Biocentric/gutportrait
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

For gutportrait this matters twice over: pin the gutportrait version **and** the ampliseq version, primers, truncation lengths and taxonomy database. Results from two different ampliseq configurations are not comparable, whatever gutportrait does with them.

First, go to the [Biocentric/gutportrait releases page](https://github.com/Biocentric/gutportrait/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Runs on the in-repo synthetic fixture, so needs no other parameters
- `test_full`
  - The same synthetic fixture without the CI resource caps. There is no public full-size dataset: a real run includes one person's reads, which cannot be redistributed
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://charliecloud.io/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow `24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#set-max-resources) and [customise process resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#customize-process-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#update-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#modifying-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
