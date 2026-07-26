from pathlib import Path

import pandas as pd
import pytest

from betasieve.analysis import Col, SieveResults, _add_flags, _add_statistics
from betasieve.cg_probe_table import DesignGroup
from betasieve.config import SieveArgs


@pytest.fixture
def betas_path(tmp_path: Path) -> Path:
    path = tmp_path / "betas.tsv"
    path.write_text("IlmnID\tSample_A\ncg00000001_TC11\t0.1\n")
    return path


@pytest.fixture
def sieve_args(betas_path: Path, tmp_path: Path) -> SieveArgs:
    return SieveArgs(
        betas_path=betas_path,
        threshold=0.1,
        out_dir=tmp_path / "results",
        report=False,
    )


@pytest.fixture
def diff_frame() -> pd.DataFrame:
    return pd.DataFrame(
        {
            Col.GROUP: [
                DesignGroup.EXACT_REPLICATES.value,
                DesignGroup.EXACT_REPLICATES.value,
                DesignGroup.PAIR_TYPE.value,
                DesignGroup.PAIR_DESIGN.value,
            ],
            "Sample_A": [0.0, 0.0, 0.2, 0.0],
            "Sample_B": [0.2, 0.0, 0.2, 0.0],
            "Sample_C": [0.0, 0.0, 0.2, 0.0],
            "Sample_D": [0.2, 0.2, 0.2, 0.0],
        },
        index=pd.Index(["cg_er1", "cg_er2", "cg_pair1", "cg_pair2"], name=Col.SITE),
    )


@pytest.fixture
def sieve_results(diff_frame: pd.DataFrame) -> SieveResults:
    statistics = _add_statistics(diff_frame, 0.1, "fdr_bh", 0.95)
    flagged = _add_flags(statistics)
    return SieveResults(
        diff_frame=diff_frame,
        threshold=0.1,
        statistics_frame=statistics.copy(),
        flagged_frame=flagged,
        candidate_cpgs=pd.Series(
            ["cg_pair1_TC11", "cg_pair1_TC21"], name="IlmnID"
        ),
    )
