import numpy as np
import pandas as pd
import pytest

from betasieve.cg_probe_table import DesignGroup
from betasieve.columns import Col
from betasieve.site_statistics import _add_flags, _add_statistics


def test_add_statistics_uses_exact_replicates_as_empirical_null(
    diff_frame: pd.DataFrame,
) -> None:
    result, _ = _add_statistics(diff_frame, 0.1, "fdr_bh", 0.95)

    assert (result[Col.N] == 4).all()
    assert (result[Col.THRESHOLD] == 0.1).all()
    np.testing.assert_allclose(result[Col.P0], 3 / 8)
    assert result.loc["cg_pair1", Col.ABOVE] == 4
    assert result.loc["cg_pair1", Col.P_HAT] == pytest.approx(1.0)
    assert result.loc["cg_pair1", Col.P_EMPIR] == pytest.approx(1 / 3)
    assert result.loc["cg_pair2", Col.P_EMPIR] == pytest.approx(1.0)
    assert result.loc["cg_er1", Col.P_ADJUSTED] == result.loc["cg_er1", Col.P_VALUE]
    assert np.isfinite(result[Col.CI_LOWER]).all()
    assert np.isfinite(result[Col.CI_UPPER]).all()


def test_add_flags_uses_strict_alpha_and_ci_comparisons() -> None:
    frame = pd.DataFrame(
        {
            Col.CONFIDENCE: [0.95, 0.95],
            Col.GROUP: [DesignGroup.PAIR_TYPE.value] * 2,
            Col.P_VALUE: [0.049, 0.051],
            Col.P_ADJUSTED: [0.049, 0.051],
            Col.P_EMPIR: [0.049, 0.051],
            Col.P_EMPIR_ADJUSTED: [0.049, 0.051],
            Col.P_BETA: [0.049, 0.051],
            Col.P_BETA_ADJUSTED: [0.049, 0.051],
            Col.CI_LOWER: [0.2, 0.1],
            Col.P0: [0.1, 0.1],
        }
    )

    result = _add_flags(frame)

    assert result[Col.P_FLAG].tolist() == [True, False]
    assert result[Col.P_ADJ_FLAG].tolist() == [True, False]
    assert result[Col.P_EMPIR_FLAG].tolist() == [True, False]
    assert result[Col.P_EMPIR_ADJ_FLAG].tolist() == [True, False]
    assert result[Col.P_BETA_FLAG].tolist() == [True, False]
    assert result[Col.P_BETA_ADJ_FLAG].tolist() == [True, False]
    assert result[Col.CI_FLAG].tolist() == [True, False]
