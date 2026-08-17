from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional

VALID_FDR_METHODS = (
    "bonferroni",
    "sidak",
    "holm-sidak",
    "holm",
    "simes-hochberg",
    "hommel",
    "fdr_bh",
    "fdr_by",
    "fdr_tsbh",
    "fdr_tsbky",
)


@dataclass
class SieveConfig:
    threshold: Optional[float] = None
    fdr: str = "fdr_bh"
    confidence: float = 0.95
    threshold_min: Optional[float] = None
    threshold_max: Optional[float] = None
    threshold_step: Optional[float] = None
    target_p0: float = 0.05


@dataclass
class ReportConfig:
    betas_path: Path
    analysis: SieveConfig = field(default_factory=SieveConfig)
    out_dir: Path = Path("results")
    report: bool = True
    pkl: bool = False
    csv_files: bool = True

    @property
    def csv_dir(self) -> Path:
        return self.out_dir / "csv"

    @property
    def figures_dir(self) -> Path:
        return self.out_dir / "figures"

    @property
    def report_dir(self) -> Path:
        return self.out_dir / "report"

    @property
    def pkl_dir(self) -> Path:
        return self.out_dir / "pkl"

    @classmethod
    def from_namespace(cls, namespace: argparse.Namespace) -> ReportConfig:
        return cls(
            betas_path=namespace.betas,
            analysis=SieveConfig(
                threshold=namespace.threshold,
                fdr=namespace.fdr,
                confidence=namespace.confidence,
                threshold_min=namespace.threshold_min,
                threshold_max=namespace.threshold_max,
                threshold_step=namespace.threshold_step,
                target_p0=getattr(namespace, "target_p0", 0.05),
            ),
            out_dir=namespace.out_dir,
            report=namespace.report,
            pkl=getattr(namespace, "pkl", False),
        )


def _sieve_config_errors(config: SieveConfig) -> List[str]:
    errors: List[str] = []

    if config.fdr not in VALID_FDR_METHODS:
        errors.append(
            f"fdr {config.fdr!r} is not supported. "
            f"Choose one of: {', '.join(VALID_FDR_METHODS)}."
        )

    if not (0.0 < config.confidence < 1.0):
        errors.append(
            "confidence must be between 0 and 1 (exclusive), "
            f"got {config.confidence}."
        )

    if not (0.0 < config.target_p0 < 1.0):
        errors.append(
            f"target_p0 must be between 0 and 1 (exclusive), got {config.target_p0}."
        )

    def check_threshold_value(name: str, value: Optional[float]) -> None:
        if value is None:
            return
        if value <= 0.0 or value > 1.0:
            errors.append(f"{name} must be in (0, 1], got {value}.")

    sweep_fields = (
        ("threshold_min", config.threshold_min),
        ("threshold_max", config.threshold_max),
        ("threshold_step", config.threshold_step),
    )
    has_threshold = config.threshold is not None
    has_all_sweep = all(v is not None for _, v in sweep_fields)
    has_any_sweep = any(v is not None for _, v in sweep_fields)

    if has_threshold:
        check_threshold_value("threshold", config.threshold)
    elif has_all_sweep:
        for name, value in sweep_fields:
            check_threshold_value(name, value)
        if config.threshold_step is not None and config.threshold_step <= 0.0:
            errors.append(f"threshold_step must be > 0, got {config.threshold_step}.")
        if (
            config.threshold_min is not None
            and config.threshold_max is not None
            and config.threshold_min >= config.threshold_max
        ):
            errors.append(
                f"threshold_min ({config.threshold_min}) must be less than "
                f"threshold_max ({config.threshold_max})."
            )
    else:
        if has_any_sweep:
            missing = [name for name, value in sweep_fields if value is None]
            errors.append(
                "threshold is not set; automatic search requires: "
                + ", ".join(missing)
                + "."
            )
        else:
            errors.append(
                "Either set threshold, or set threshold_min, threshold_max, "
                "and threshold_step for automatic threshold search."
            )

    return errors


def _raise_validation_errors(title: str, errors: List[str]) -> None:
    if errors:
        message = f"{title}:\n" + "\n".join(f"  • {err}" for err in errors)
        raise ValueError(message)


def validate_sieve_config(config: SieveConfig) -> None:
    _raise_validation_errors(
        "Invalid betaSieve analysis configuration",
        _sieve_config_errors(config),
    )


def validate_report_config(args: ReportConfig) -> None:
    errors: List[str] = []
    betas_path = args.betas_path
    if not isinstance(betas_path, Path):
        errors.append(
            f"betas_path must be a pathlib.Path, got {type(betas_path).__name__}."
        )
    elif not betas_path.exists():
        errors.append(f"betas file does not exist: {betas_path}")
    elif not betas_path.is_file():
        errors.append(f"betas path is not a file: {betas_path}")

    errors.extend(_sieve_config_errors(args.analysis))
    _raise_validation_errors("Invalid arguments for betaSieve", errors)
