from __future__ import annotations

from typing import List, Optional

import tyro

from .config import ReportConfig
from .pipeline import run_beta_sieve


def main(argv: Optional[List[str]] = None) -> None:
    run_beta_sieve(tyro.cli(ReportConfig, args=argv))


if __name__ == "__main__":
    main()
