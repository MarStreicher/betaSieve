from __future__ import annotations

import argparse
from dataclasses import dataclass
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
class SieveArgs:
    betas_path: Path
    threshold: Optional[float] = None
    fdr: str = "fdr_bh"
    confidence: float = 0.95
    threshold_min: Optional[float] = None
    threshold_max: Optional[float] = None
    threshold_step: Optional[float] = None
    target_p0: float = 0.05
    out_dir: Path = Path("results")
    report: bool = True
    pkl: bool = False

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
    def from_namespace(cls, namespace: argparse.Namespace) -> SieveArgs:
        return cls(
            betas_path=namespace.betas,
            threshold=namespace.threshold,
            fdr=namespace.fdr,
            confidence=namespace.confidence,
            threshold_min=namespace.threshold_min,
            threshold_max=namespace.threshold_max,
            threshold_step=namespace.threshold_step,
            out_dir=namespace.out_dir,
            report=namespace.report,
            pkl=getattr(namespace, "pkl", False),
        )


def validate_sieve_args(args: SieveArgs) -> None:
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

    if args.fdr not in VALID_FDR_METHODS:
        errors.append(
            f"fdr {args.fdr!r} is not supported. "
            f"Choose one of: {', '.join(VALID_FDR_METHODS)}."
        )

    if not (0.0 < args.confidence < 1.0):
        errors.append(
            f"confidence must be between 0 and 1 (exclusive), got {args.confidence}."
        )

    if not (0.0 < args.target_p0 < 1.0):
        errors.append(
            f"target_p0 must be between 0 and 1 (exclusive), got {args.target_p0}."
        )

    def check_threshold_value(name: str, value: Optional[float]) -> None:
        if value is None:
            return
        if value <= 0.0 or value > 1.0:
            errors.append(f"{name} must be in (0, 1], got {value}.")

    sweep_fields = (
        ("threshold_min", args.threshold_min),
        ("threshold_max", args.threshold_max),
        ("threshold_step", args.threshold_step),
    )
    has_threshold = args.threshold is not None
    has_all_sweep = all(v is not None for _, v in sweep_fields)
    has_any_sweep = any(v is not None for _, v in sweep_fields)

    if has_threshold:
        check_threshold_value("threshold", args.threshold)
    elif has_all_sweep:
        for name, value in sweep_fields:
            check_threshold_value(name, value)
        if args.threshold_step is not None and args.threshold_step <= 0.0:
            errors.append(f"threshold_step must be > 0, got {args.threshold_step}.")
        if (
            args.threshold_min is not None
            and args.threshold_max is not None
            and args.threshold_min >= args.threshold_max
        ):
            errors.append(
                f"threshold_min ({args.threshold_min}) must be less than "
                f"threshold_max ({args.threshold_max})."
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

    if errors:
        message = "Invalid arguments for betaSieve:\n" + "\n".join(
            f"  • {err}" for err in errors
        )
        raise ValueError(message)
