# Biocentric/gutportrait: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [unreleased<!-- TODO nf-core: replace with date on release -->]

Initial release of Biocentric/gutportrait, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- Samplesheet schema with `sample`, `role` (personal / panel), `study` and `run`, plus
  optional `age_cat`, `sex`, `bmi_cat` and `country` for cohort matching. Validation
  rejects a study that mixes personal and panel samples, and panels smaller than
  `--min_panel_samples`.
- `--ampliseq_results`, `--taxonomy_db` and `--arm`: the pipeline reads
  `dada2/ASV_table.tsv` and the named-rank `dada2/ASV_tax_species.<db>.tsv` from an
  nf-core/ampliseq results directory.
- `COLLAPSE_GENUS`: genus collapse that keeps unassigned reads as an explicit
  `Unassigned` row, with per-sample accounting.
- `TAXON_PERCENTILES`: per-genus centred-log-ratio percentile of each personal sample
  against the panel, with a panel-only prevalence filter (`--min_prevalence`) and an
  additive pseudocount (`--pseudocount`).
- Synthetic test fixture with planted signals (age effect, study bias, deviant sample,
  missing metadata), generated deterministically by
  `tests/fixtures/make_synthetic_fixture.py`, and used by `-profile test`.
- Tool citations in the MultiQC methods description.
- `BATCH_DIAGNOSE`: PERMANOVA (vegan `adonis2`) on Aitchison distance for the panel:
  R²(study) with an orientation-only banner above `--banner_r2_study`, and study +
  variable models with marginal R² per term (the batch-dominance rule of D12).
- `ALPHA_DIVERSITY`: Shannon, observed ASVs and Pielou evenness at a fixed rarefaction
  depth (`--rarefaction_depth`, `--rarefaction_iterations`, `--seed`), with each personal
  sample's position in the pooled panel and against every panel study.
- `DEVIATION_FLAGS`: the extremes-only checks of D12 (low alpha, single-taxon dominance,
  Enterobacteriaceae and Proteobacteria blooms, Faecalibacterium depletion, missing core
  genera, beta outlier), each reported raised or not, with value and threshold. All
  thresholds are parameters.
- Optional samplesheet column `context`, carried into every personal-sample table.
- `BATCH_CORRECT`: MMUPHin `adjust_batch` harmonises the panel studies, optionally
  protecting complete covariates (`--correction_covariates`). Personal samples are never
  corrected.
- `POSITIVE_CONTROL`: injected-signal control (D12a). A planted fold change must survive the
  same correction (paired Wilcoxon, `--pc_alpha`, `--pc_retain`), or the run stops.
- `BATCH_DOMINANCE`: the study banner, a finding / orientation-only verdict per contrast, and
  personal positioning against the between-study spread (D12).
- `PERCENTILE_STABILITY` and `TAXON_PERCENTILES_CORRECTED`: every genus position against
  the uncorrected and the corrected panel; positions that change category are unstable.
- Samplesheet role `repeat` and column `subject`: a second sample of a panel person, used
  only as the within-person reference. Validation rejects a repeat that names no panel
  person of its study.
- `WITHIN_PERSON`: differences between the personal samples against a noise floor borrowed
  from the panel's repeat pairs (distance percentiles, change in Shannon diversity, per-genus
  change in units of its ordinary within-person change, with a leave-one-out count of what
  changes by chance; an absence counts only where the shallower sample could have detected
  the genus), and a sample-identity check that tells a large change in a few genera from a
  different person (`--change_z`, `--min_reference_pairs`, `--identity_within_percentile`,
  `--identity_between_percentile`).
- `TECHNICAL_QC`: retention at every step of ampliseq's read tracking, analysed reads and the
  unassigned fraction, against the panel; atypical beyond a modified z-score of `--qc_max_z`
  against the pooled panel and every study.
- `DEVIATION_FLAGS` calibration: every check is run leave-one-out on the panel, and each flag
  carries its `reference_fire_rate`; `flag_calibration.tsv` lists the rates. New check
  `excess_extremes`: more genera in the outer percentile tails than panel people have.
- Test fixture: the three personal samples are now one invented person; 16 panel people have
  a repeat sample; read tracking (`overall_summary.tsv`) with one technical outlier; and two
  more samplesheets, one without repeats and one with a sample mix-up.

### `Changed`

- `beta_outlier` compares a sample with the panel's leave-one-out distances to its centroid,
  so that the reference distances are measured the way a personal sample's is.

### `Fixed`

- `bin/` scripts are now committed executable.

### `Dependencies`

- `r-vegan` 2.7.5 (R 4.5.3) for `BATCH_DIAGNOSE`.
- `bioconductor-mmuphin` 1.16.0 (R 4.3.2) for `BATCH_CORRECT` and `POSITIVE_CONTROL`.

### `Deprecated`
