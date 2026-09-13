# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased] - 2026-09-13

### Added

- In-memory analysis through `sieve_betas` and `SieveConfig`.

### Changed

- File-based runs use `PipelineConfig` and `run_sieve_pipeline`.
- The CLI is generated from the config dataclasses with
  [tyro](https://brentyi.github.io/tyro/). Option names changed
  (`--betas` → `--betas-path`, `--threshold` → `--analysis.threshold`).
- `SieveConfig` now defaults to an automatic threshold search.

### Removed

- Hand-written argparse CLI (`build_parser`, `from_namespace`).

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
