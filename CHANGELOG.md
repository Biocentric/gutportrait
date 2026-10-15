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

### `Fixed`

- `bin/` scripts are now committed executable.

### `Dependencies`

### `Deprecated`
