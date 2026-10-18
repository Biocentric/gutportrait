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

### `Fixed`

- `bin/` scripts are now committed executable.

### `Dependencies`

- `r-vegan` 2.7.5 (R 4.5.3) for `BATCH_DIAGNOSE`.

### `Deprecated`
