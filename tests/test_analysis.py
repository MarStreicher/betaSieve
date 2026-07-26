from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from betasieve import analysis
from betasieve.analysis import (
    Col,
    _add_flags,
    _add_statistics,
    _collect_max_min_differences,
    _create_cpg_list,
    _diff_value_columns,
    _find_threshold,
    _flag_rates_by_group,
    _select_threshold_for_replicates,
    _sweep_thresholds,
    run_duplicate_analysis,
)
from betasieve.cg_probe_table import CgProbeTable, DesignGroup
from betasieve.config import SieveArgs


def test_diff_value_columns_returns_only_numeric_sample_columns(
    diff_frame: pd.DataFrame,
) -> None:
    frame = diff_frame.copy()
    frame[Col.N] = 4
    frame["note"] = "text"

    assert _diff_value_columns(frame) == [
        "Sample_A",
        "Sample_B",
        "Sample_C",
        "Sample_D",
    ]


def test_collect_max_min_differences_for_designs_and_replicates() -> None:
    cg_by_sample = pd.DataFrame(
        {
            "A": [0.1, 0.4, 0.2, 0.5],
            "B": [0.2, 0.3, 0.6, 0.1],
        },
        index=["cg1_TC11", "cg1_TC21", "cg2_BC11", "cg2_BC12"],
    )
    groups = CgProbeTable.from_probe_ids(pd.Series(cg_by_sample.index))

    result = _collect_max_min_differences(cg_by_sample, groups)

    assert result.index.tolist() == ["cg1", "cg2"]
    assert result[Col.GROUP].tolist() == [
        DesignGroup.PAIR_TYPE.value,
        DesignGroup.EXACT_REPLICATES.value,
    ]
    np.testing.assert_allclose(
        result[["A", "B"]].to_numpy(),
        [[0.3, 0.1], [0.3, 0.5]],
    )


def test_add_statistics_uses_exact_replicates_as_empirical_null(
    diff_frame: pd.DataFrame,
) -> None:
    result = _add_statistics(diff_frame, 0.1, "fdr_bh", 0.95)

    assert (result[Col.N] == 4).all()
    assert (result[Col.THRESHOLD] == 0.1).all()
    np.testing.assert_allclose(result[Col.P0], 3 / 8)
    assert result.loc["cg_pair1", Col.ABOVE] == 4
    assert result.loc["cg_pair1", Col.P_HAT] == pytest.approx(1.0)
    assert result.loc["cg_pair1", Col.P_BETA] == pytest.approx(1 / 3)
    assert result.loc["cg_pair2", Col.P_BETA] == pytest.approx(1.0)
    assert result.loc["cg_er1", Col.P_ADJUSTED] == result.loc["cg_er1", Col.P_VALUE]
    assert np.isfinite(result[Col.CI_LOWER]).all()
    assert np.isfinite(result[Col.CI_UPPER]).all()


def test_add_flags_uses_strict_alpha_and_ci_comparisons() -> None:
    frame = pd.DataFrame(
        {
            Col.CONFIDENCE: [0.95, 0.95],
            Col.P_VALUE: [0.049, 0.051],
            Col.P_ADJUSTED: [0.049, 0.051],
            Col.P_BETA: [0.049, 0.051],
            Col.P_BETA_ADJUSTED: [0.049, 0.051],
            Col.CI_LOWER: [0.2, 0.1],
            Col.P0: [0.1, 0.1],
        }
    )

    result = _add_flags(frame)

    assert result[Col.P_FLAG].tolist() == [True, False]
    assert result[Col.P_ADJ_FLAG].tolist() == [True, False]
    assert result[Col.P_BETA_FLAG].tolist() == [True, False]
    assert result[Col.P_BETA_ADJ_FLAG].tolist() == [True, False]
    assert result[Col.CI_FLAG].tolist() == [True, False]


def test_create_cpg_list_returns_all_probe_instances_at_flagged_sites() -> None:
    flagged = pd.DataFrame(
        {Col.P_BETA_ADJ_FLAG: [True, False]},
        index=["cg1", "cg2"],
    )
    cg_by_sample = pd.DataFrame(
        {"A": [0.1, 0.2, 0.3, 0.4]},
        index=["cg1_TC11", "cg1_TC21", "cg2_TC11", "invalid"],
    )

    result = _create_cpg_list(flagged, cg_by_sample)

    assert result.tolist() == ["cg1_TC11", "cg1_TC21"]
    assert result.name == "IlmnID"


def test_create_cpg_list_raises_when_no_candidates_match() -> None:
    flagged = pd.DataFrame(
        {Col.P_BETA_ADJ_FLAG: [True]},
        index=["cg_missing"],
    )
    cg_by_sample = pd.DataFrame({"A": [0.1]}, index=["cg1_TC11"])

    with pytest.raises(ValueError, match="No candidate CpGs"):
        _create_cpg_list(flagged, cg_by_sample)


def test_flag_rates_by_group_calculates_counts_and_percentages() -> None:
    frame = pd.DataFrame(
        {
            Col.GROUP: ["g1", "g1", "g2"],
            Col.CI_FLAG: [True, False, True],
            Col.P_ADJ_FLAG: [False, False, True],
            Col.P_BETA_ADJ_FLAG: [True, True, False],
        }
    )

    result = _flag_rates_by_group(frame).set_index(Col.GROUP)

    assert result.loc["g1", "n_sites"] == 2
    assert result.loc["g1", Col.PCT_CI_FLAGGED] == pytest.approx(50.0)
    assert result.loc["g1", Col.PCT_P_ADJ_FLAGGED] == pytest.approx(0.0)
    assert result.loc["g1", Col.PCT_BETA_ADJ_FLAGGED] == pytest.approx(100.0)


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

    assert result[Col.THRESHOLD].unique().tolist() == pytest.approx(
        [0.05, 0.1, 0.15]
    )
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
    sweep = pd.DataFrame(
        {Col.GROUP: ["other"], Col.THRESHOLD: [0.1], Col.P0: [0.1]}
    )

    with pytest.raises(ValueError, match="No sweep results"):
        _select_threshold_for_replicates(sweep, target_p0=0.05)


def test_find_threshold_writes_sweep_csv(
    diff_frame: pd.DataFrame, tmp_path: Path
) -> None:
    chosen, sweep = _find_threshold(
        diff_frame,
        threshold_min=0.05,
        threshold_max=0.15,
        threshold_step=0.05,
        fdr="fdr_bh",
        confidence=0.95,
        target_p0=0.4,
        out_dir=tmp_path / "csv",
    )

    assert chosen == pytest.approx(0.05)
    written = pd.read_csv(tmp_path / "csv" / "threshold_sweep_summary.csv")
    pd.testing.assert_frame_equal(written, sweep, check_dtype=False)


def test_run_duplicate_analysis_orchestrates_fixed_threshold(
    sieve_args: SieveArgs,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    cg_data = pd.DataFrame({"A": [0.1]}, index=["cg1_TC11"])
    diff = pd.DataFrame({Col.GROUP: ["group"], "A": [0.2]}, index=["cg1"])
    statistics = pd.DataFrame({"stat": [1]}, index=["cg1"])
    flagged = pd.DataFrame({Col.P_BETA_ADJ_FLAG: [True]}, index=["cg1"])
    candidates = pd.Series(["cg1_TC11"], name="IlmnID")

    class FakeLoader:
        def __init__(self, path):
            assert path == sieve_args.betas_path

        def load_data(self):
            return cg_data

    monkeypatch.setattr(analysis, "BetasLoader", FakeLoader)
    monkeypatch.setattr(
        analysis.CgProbeTable, "from_probe_ids", lambda ids: pd.DataFrame()
    )
    monkeypatch.setattr(
        analysis, "_collect_max_min_differences", lambda data, groups: diff
    )
    monkeypatch.setattr(
        analysis, "_add_statistics", lambda frame, threshold, fdr, confidence: statistics
    )
    monkeypatch.setattr(analysis, "_add_flags", lambda frame: flagged)
    monkeypatch.setattr(
        analysis, "_create_cpg_list", lambda frame, data: candidates
    )

    result = run_duplicate_analysis(sieve_args)

    assert result.diff_frame is diff
    assert result.threshold == 0.1
    assert result.statistics_frame is statistics
    assert result.flagged_frame is flagged
    assert result.candidate_cpgs is candidates
    assert result.sweep_df is None
