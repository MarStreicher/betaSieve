import argparse
from pathlib import Path

import pytest

from betasieve.config import SieveArgs, validate_sieve_args


def test_output_directory_properties_are_derived_from_root(
    sieve_args: SieveArgs,
) -> None:
    assert sieve_args.csv_dir == sieve_args.out_dir / "csv"
    assert sieve_args.figures_dir == sieve_args.out_dir / "figures"
    assert sieve_args.report_dir == sieve_args.out_dir / "report"
    assert sieve_args.pkl_dir == sieve_args.out_dir / "pkl"


def test_from_namespace_maps_cli_values_and_defaults_pkl() -> None:
    namespace = argparse.Namespace(
        betas=Path("betas.tsv"),
        threshold=0.2,
        fdr="holm",
        confidence=0.9,
        threshold_min=0.01,
        threshold_max=0.1,
        threshold_step=0.01,
        out_dir=Path("output"),
        report=False,
    )

    args = SieveArgs.from_namespace(namespace)

    assert args == SieveArgs(
        betas_path=Path("betas.tsv"),
        threshold=0.2,
        fdr="holm",
        confidence=0.9,
        threshold_min=0.01,
        threshold_max=0.1,
        threshold_step=0.01,
        out_dir=Path("output"),
        report=False,
        pkl=False,
    )


def test_validate_accepts_fixed_threshold(sieve_args: SieveArgs) -> None:
    validate_sieve_args(sieve_args)


def test_validate_accepts_complete_threshold_sweep(
    sieve_args: SieveArgs,
) -> None:
    sieve_args.threshold = None
    sieve_args.threshold_min = 0.01
    sieve_args.threshold_max = 0.1
    sieve_args.threshold_step = 0.01

    validate_sieve_args(sieve_args)


@pytest.mark.parametrize(
    ("updates", "message"),
    [
        ({"betas_path": "betas.tsv"}, "betas_path must be a pathlib.Path"),
        ({"betas_path": Path("missing.tsv")}, "betas file does not exist"),
        ({"fdr": "invalid"}, "fdr 'invalid' is not supported"),
        ({"confidence": 0.0}, "confidence must be between 0 and 1"),
        ({"confidence": 1.0}, "confidence must be between 0 and 1"),
        ({"target_p0": 0.0}, "target_p0 must be between 0 and 1"),
        ({"threshold": 0.0}, r"threshold must be in \(0, 1\]"),
        ({"threshold": 1.1}, r"threshold must be in \(0, 1\]"),
    ],
)
def test_validate_rejects_invalid_individual_values(
    sieve_args: SieveArgs, updates: dict, message: str
) -> None:
    for name, value in updates.items():
        setattr(sieve_args, name, value)

    with pytest.raises(ValueError, match=message):
        validate_sieve_args(sieve_args)


def test_validate_rejects_directory_as_betas_path(
    sieve_args: SieveArgs, tmp_path: Path
) -> None:
    sieve_args.betas_path = tmp_path

    with pytest.raises(ValueError, match="betas path is not a file"):
        validate_sieve_args(sieve_args)


def test_validate_reports_missing_sweep_fields(sieve_args: SieveArgs) -> None:
    sieve_args.threshold = None
    sieve_args.threshold_min = 0.01
    sieve_args.threshold_max = None
    sieve_args.threshold_step = None

    with pytest.raises(ValueError) as exc_info:
        validate_sieve_args(sieve_args)

    assert "automatic search requires: threshold_max, threshold_step" in str(
        exc_info.value
    )


def test_validate_requires_fixed_threshold_or_sweep(sieve_args: SieveArgs) -> None:
    sieve_args.threshold = None
    sieve_args.threshold_min = None
    sieve_args.threshold_max = None
    sieve_args.threshold_step = None

    with pytest.raises(ValueError, match="Either set threshold"):
        validate_sieve_args(sieve_args)


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
    sieve_args: SieveArgs,
    minimum: float,
    maximum: float,
    step: float,
    message: str,
) -> None:
    sieve_args.threshold = None
    sieve_args.threshold_min = minimum
    sieve_args.threshold_max = maximum
    sieve_args.threshold_step = step

    with pytest.raises(ValueError, match=message):
        validate_sieve_args(sieve_args)


def test_validation_collects_multiple_errors(sieve_args: SieveArgs) -> None:
    sieve_args.fdr = "bad"
    sieve_args.confidence = 2.0
    sieve_args.threshold = -1.0

    with pytest.raises(ValueError) as exc_info:
        validate_sieve_args(sieve_args)

    message = str(exc_info.value)
    assert message.startswith("Invalid arguments for betaSieve:")
    assert "fdr 'bad'" in message
    assert "confidence must be" in message
    assert "threshold must be" in message
