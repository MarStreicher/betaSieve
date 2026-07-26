from __future__ import annotations

import argparse
from pathlib import Path

from .config import VALID_FDR_METHODS, SieveArgs
from .pipeline import run_beta_sieve


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Compute max-min beta value differences for duplicate CpG probe groups "
            "and flag sites with statistically significant variability."
        )
    )
    parser.add_argument(
        "--betas",
        type=Path,
        required=True,
        help="Path to betas CSV (IlmnID × samples).",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=None,
        help="Fixed β-values max–min threshold. If omitted, sweep threshold_min…max.",
    )
    parser.add_argument(
        "--confidence",
        type=float,
        default=0.95,
        help="Confidence level for intervals and p-value flags (default: 0.95).",
    )
    parser.add_argument(
        "--fdr",
        type=str,
        default="fdr_bh",
        choices=VALID_FDR_METHODS,
        help="Multiple-testing method (statsmodels multipletests).",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("results"),
        dest="out_dir",
        help="Root output directory (default: ./results; creates csv/ and report/ subdirs).",
    )
    parser.add_argument(
        "--report",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Generate the HTML analysis report (default: enabled).",
    )
    parser.add_argument(
        "--pkl",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Write SieveArgs and SieveResults pickles to out-dir/pkl/.",
    )
    parser.add_argument(
        "--threshold-min",
        type=float,
        default=0.01,
        dest="threshold_min",
        help="Lower bound for automatic threshold search (inclusive).",
    )
    parser.add_argument(
        "--threshold-max",
        type=float,
        default=0.1,
        dest="threshold_max",
        help="Upper bound for automatic threshold search (inclusive).",
    )
    parser.add_argument(
        "--threshold-step",
        type=float,
        default=0.01,
        dest="threshold_step",
        help="Step size for automatic threshold search.",
    )
    return parser


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    namespace = parser.parse_args(argv)
    run_beta_sieve(SieveArgs.from_namespace(namespace))


if __name__ == "__main__":
    main()
