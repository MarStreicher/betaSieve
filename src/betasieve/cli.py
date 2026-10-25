from __future__ import annotations

from typing import List, Optional

import tyro

from .config import PipelineConfig
from .pipeline import run_sieve_pipeline


def main(argv: Optional[List[str]] = None) -> None:
    run_sieve_pipeline(tyro.cli(PipelineConfig, args=argv))


if __name__ == "__main__":
    main()
