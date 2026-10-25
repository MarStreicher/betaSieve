# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased] - 2026-XX-XX

### Added

### Changed

### Removed

## [0.1.0] - 2026-07-26

### Added

- Duplicate Analysis of EPIC V2 beta files
- Initial public release

### Changed
/

### Fixed
/

## [0.2.0] - 2026-07-26

### Added
/

### Changed

- Change version for PyPI support

### Fixed
/

## [1.0.0] - 2026-10-25

### Added

- In-memory analysis through `sieve_betas` and `SieveConfig`.
- Example pipeline script `example/usage_betaSieve.py`. Local paths are read from a gitignored `.env` file; `.env.example` lists the variables.

### Changed

- File-based runs use `PipelineConfig` and `run_sieve_pipeline`.
- The CLI is generated from the config dataclasses with
  [tyro](https://brentyi.github.io/tyro/). Option names changed
  (`--betas` → `--betas-path`, `--threshold` → `--analysis.threshold`).
- `SieveConfig` now defaults to an automatic threshold search.
- Report design.
- Report wording uses design replicates vs exact replicates.
- Difference and p̂ histograms support Absolute / Percentage toggles.
- Report subsections implement `_figures()` and `generate()` is shared in the base class.
- Report and README text uses "CpG site" instead of "CpG-site", and difference plots are labelled "Max-min β-value difference".
- Difference histograms extend to the largest observed difference instead of stopping at 0.3.
- Report figures download as SVG.
- Require `epicv2-io>=1.0.5` and `python-dotenv`. `matplotlib` is no longer a package dependency.
- The package version is taken from git tags with setuptools-scm.
- The release workflow runs tests on a version tag, publishes to PyPI, and creates a GitHub release.

### Removed

- Hand-written argparse CLI (`build_parser`, `from_namespace`).
