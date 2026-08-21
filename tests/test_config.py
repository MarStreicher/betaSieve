from pathlib import Path
from dataclasses import replace

import pytest

from betasieve.config import (
    SieveConfig,
    ReportConfig,
    validate_sieve_config,
    validate_report_config,
)


def test_output_directory_properties_are_derived_from_root(
    sieve_args: ReportConfig,
) -> None:
    assert sieve_args.csv_dir == sieve_args.out_dir / "csv"
    assert sieve_args.figures_dir == sieve_args.out_dir / "figures"
    assert sieve_args.report_dir == sieve_args.out_dir / "report"
    assert sieve_args.pkl_dir == sieve_args.out_dir / "pkl"


def test_validate_accepts_fixed_threshold(sieve_args: ReportConfig) -> None:
    validate_report_config(sieve_args)


def test_validate_accepts_complete_threshold_sweep(
    sieve_args: ReportConfig,
) -> None:
    sieve_args.analysis = SieveConfig(
        threshold_min=0.01,
        threshold_max=0.1,
        threshold_step=0.01,
    )

    validate_report_config(sieve_args)


@pytest.mark.parametrize(
    ("updates", "message"),
    [
        ({"fdr": "invalid"}, "fdr 'invalid' is not supported"),
        ({"confidence": 0.0}, "confidence must be between 0 and 1"),
        ({"confidence": 1.0}, "confidence must be between 0 and 1"),
        ({"target_p0": 0.0}, "target_p0 must be between 0 and 1"),
        ({"threshold": 0.0}, r"threshold must be in \(0, 1\]"),
        ({"threshold": 1.1}, r"threshold must be in \(0, 1\]"),
    ],
)
def test_validate_rejects_invalid_individual_values(
    sieve_args: ReportConfig, updates: dict, message: str
) -> None:
    analysis = replace(sieve_args.analysis, **updates)

    with pytest.raises(ValueError, match=message):
        validate_sieve_config(analysis)


@pytest.mark.parametrize(
    ("betas_path", "message"),
    [
        ("betas.tsv", "betas_path must be a pathlib.Path"),
        (Path("missing.tsv"), "betas file does not exist"),
    ],
)
def test_validate_rejects_invalid_betas_path(
    sieve_args: ReportConfig, betas_path: object, message: str
) -> None:
    sieve_args.betas_path = betas_path  # type: ignore[assignment]

    with pytest.raises(ValueError, match=message):
        validate_report_config(sieve_args)


def test_validate_rejects_directory_as_betas_path(
    sieve_args: ReportConfig, tmp_path: Path
) -> None:
    sieve_args.betas_path = tmp_path

    with pytest.raises(ValueError, match="betas path is not a file"):
        validate_report_config(sieve_args)


def test_default_config_is_a_valid_threshold_sweep() -> None:
    validate_sieve_config(SieveConfig())


def test_validate_reports_missing_sweep_fields(sieve_args: ReportConfig) -> None:
    sieve_args.analysis = SieveConfig(
        threshold_min=0.01,
        threshold_max=None,
        threshold_step=None,
    )

    with pytest.raises(ValueError) as exc_info:
        validate_sieve_config(sieve_args.analysis)

    assert "automatic search requires: threshold_max, threshold_step" in str(
        exc_info.value
    )


def test_validate_requires_fixed_threshold_or_sweep(sieve_args: ReportConfig) -> None:
    sieve_args.analysis = SieveConfig(
        threshold_min=None,
        threshold_max=None,
        threshold_step=None,
    )

    with pytest.raises(ValueError, match="Either set threshold"):
        validate_sieve_config(sieve_args.analysis)


@pytest.mark.parametrize(
    ("minimum", "maximum", "step", "message"),
    [
        (0.1, 0.1, 0.01, "threshold_min \\(0.1\\) must be less"),
        (0.2, 0.1, 0.01, "threshold_min \\(0.2\\) must be less"),
        (0.01, 0.1, 0.0, "threshold_step must be in"),
        (0.01, 0.1, 1.1, "threshold_step must be in"),
    ],
)
def test_validate_rejects_invalid_sweep_ranges(
    sieve_args: ReportConfig,
    minimum: float,
    maximum: float,
    step: float,
    message: str,
) -> None:
    sieve_args.analysis = SieveConfig(
        threshold_min=minimum,
        threshold_max=maximum,
        threshold_step=step,
    )

    with pytest.raises(ValueError, match=message):
        validate_sieve_config(sieve_args.analysis)


def test_validation_collects_multiple_errors(sieve_args: ReportConfig) -> None:
    sieve_args = replace(
        sieve_args,
        analysis=replace(
            sieve_args.analysis,
            fdr="bad",  # type: ignore[arg-type]
            confidence=2.0,
            threshold=-1.0,
        ),
    )

    with pytest.raises(ValueError) as exc_info:
        validate_report_config(sieve_args)

    message = str(exc_info.value)
    assert message.startswith("Invalid arguments for betaSieve:")
    assert "fdr 'bad'" in message
    assert "confidence must be" in message
    assert "threshold must be" in message
