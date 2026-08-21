import pandas as pd
import pytest

from betasieve.cg_probe_table import DesignGroup
from betasieve.columns import Col
from betasieve.threshold import (
    _find_threshold,
    _flag_rates_by_group,
    _select_threshold_for_replicates,
    _sweep_thresholds,
)


def test_flag_rates_by_group_calculates_counts_and_percentages() -> None:
    frame = pd.DataFrame(
        {
            Col.GROUP: ["g1", "g1", "g2"],
            Col.CI_FLAG: [True, False, True],
            Col.P_ADJ_FLAG: [False, False, True],
            Col.P_EMPIR_ADJ_FLAG: [True, True, False],
        }
    )

    result = _flag_rates_by_group(frame).set_index(Col.GROUP)

    assert result.loc["g1", "n_sites"] == 2
    assert result.loc["g1", Col.PCT_CI_FLAGGED] == pytest.approx(50.0)
    assert result.loc["g1", Col.PCT_P_ADJ_FLAGGED] == pytest.approx(0.0)
    assert result.loc["g1", Col.PCT_EMPIR_ADJ_FLAGGED] == pytest.approx(100.0)


def test_sweep_thresholds_includes_end_point(
    diff_frame: pd.DataFrame,
) -> None:
    result = _sweep_thresholds(
        diff_frame,
        threshold_min=0.05,
        threshold_max=0.15,
        threshold_step=0.05,
        fdr="fdr_bh",
        confidence=0.95,
    )

    assert result[Col.THRESHOLD].unique().tolist() == pytest.approx([0.05, 0.1, 0.15])
    assert DesignGroup.EXACT_REPLICATES.value in result[Col.GROUP].values


def test_select_threshold_returns_first_qualifying_replicate_threshold() -> None:
    sweep = pd.DataFrame(
        {
            Col.GROUP: [DesignGroup.EXACT_REPLICATES.value] * 3,
            Col.THRESHOLD: [0.2, 0.1, 0.3],
            Col.P0: [0.04, 0.08, 0.01],
        }
    )

    assert _select_threshold_for_replicates(sweep, target_p0=0.05) == 0.2


def test_select_threshold_falls_back_to_lowest_p0(
    capsys: pytest.CaptureFixture[str],
) -> None:
    sweep = pd.DataFrame(
        {
            Col.GROUP: [DesignGroup.EXACT_REPLICATES.value] * 2,
            Col.THRESHOLD: [0.1, 0.2],
            Col.P0: [0.2, 0.1],
        }
    )

    chosen = _select_threshold_for_replicates(sweep, target_p0=0.05)

    assert chosen == 0.2
    assert "target_p0=0.05 was not reached" in capsys.readouterr().out


def test_select_threshold_requires_exact_replicate_results() -> None:
    sweep = pd.DataFrame({Col.GROUP: ["other"], Col.THRESHOLD: [0.1], Col.P0: [0.1]})

    with pytest.raises(ValueError, match="No sweep results"):
        _select_threshold_for_replicates(sweep, target_p0=0.05)


def test_find_threshold_returns_sweep_results(diff_frame: pd.DataFrame) -> None:
    chosen, sweep = _find_threshold(
        diff_frame,
        threshold_min=0.05,
        threshold_max=0.15,
        threshold_step=0.05,
        fdr="fdr_bh",
        confidence=0.95,
        target_p0=0.4,
    )

    assert chosen == pytest.approx(0.05)
    assert not sweep.empty
