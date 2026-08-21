import numpy as np
import pandas as pd
import pytest

from betasieve import analysis
from betasieve.analysis import (
    _collect_max_min_differences,
    _create_cpg_list,
    sieve_betas,
    validate_betas_frame,
)
from betasieve.cg_probe_table import CgProbeTable, DesignGroup
from betasieve.columns import Col
from betasieve.config import SieveConfig


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


def test_create_cpg_list_returns_all_probe_instances_at_flagged_sites() -> None:
    flagged = pd.DataFrame(
        {Col.P_EMPIR_ADJ_FLAG: [True, False]},
        index=["cg1", "cg2"],
    )
    cg_by_sample = pd.DataFrame(
        {"A": [0.1, 0.2, 0.3, 0.4]},
        index=["cg1_TC11", "cg1_TC21", "cg2_TC11", "invalid"],
    )

    result = _create_cpg_list(flagged, cg_by_sample)

    assert result.tolist() == ["cg1_TC11", "cg1_TC21"]
    assert result.name == "IlmnID"


def test_create_cpg_list_warns_when_no_candidates_match() -> None:
    flagged = pd.DataFrame(
        {Col.P_EMPIR_ADJ_FLAG: [True]},
        index=["cg_missing"],
    )
    cg_by_sample = pd.DataFrame({"A": [0.1]}, index=["cg1_TC11"])

    with pytest.warns(UserWarning, match="No candidate CpGs"):
        result = _create_cpg_list(flagged, cg_by_sample)

    assert result.empty
    assert result.name == "IlmnID"


@pytest.mark.parametrize(
    ("frame", "message"),
    [
        (pd.DataFrame(), "at least one probe and one sample"),
        (
            pd.DataFrame({"A": ["not-a-number"]}, index=["cg00000001_TC11"]),
            "all sample columns must be numeric",
        ),
        (
            pd.DataFrame({"A": [1.1]}, index=["cg00000001_TC11"]),
            "between 0 and 1",
        ),
        (
            pd.DataFrame(
                {"A": pd.Series([pd.NA], dtype="Float64").array},
                index=["cg00000001_TC11"],
            ),
            "must not contain missing",
        ),
        (
            pd.DataFrame({"A": [0.1]}, index=["invalid"]),
            "valid EPICv2 IlmnID",
        ),
    ],
)
def test_validate_betas_frame_rejects_invalid_frames(
    frame: pd.DataFrame, message: str
) -> None:
    with pytest.raises(ValueError, match=message):
        validate_betas_frame(frame)


def test_sieve_betas_does_not_load_from_disk(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    frame = pd.DataFrame({"A": [0.1]}, index=["cg00000001_TC11"])
    expected = object()

    monkeypatch.setattr(analysis, "validate_sieve_config", lambda config: None)
    monkeypatch.setattr(analysis, "validate_betas_frame", lambda data: None)
    monkeypatch.setattr(
        analysis.CgProbeTable, "from_probe_ids", lambda ids: pd.DataFrame()
    )
    monkeypatch.setattr(
        analysis,
        "_collect_max_min_differences",
        lambda data, groups: pd.DataFrame(
            {
                Col.GROUP: [
                    DesignGroup.EXACT_REPLICATES.value,
                    DesignGroup.PAIR_TYPE.value,
                ],
                "A": [0.0, 0.2],
            },
            index=["cg0", "cg1"],
        ),
    )
    monkeypatch.setattr(
        analysis,
        "_add_statistics",
        lambda *args: (pd.DataFrame(), expected),
    )
    monkeypatch.setattr(
        analysis,
        "_add_flags",
        lambda frame: pd.DataFrame({Col.P_EMPIR_ADJ_FLAG: [False]}, index=["cg1"]),
    )

    result = sieve_betas(frame, SieveConfig(threshold=0.1))

    assert result.null_models is expected
    pd.testing.assert_frame_equal(result.sieved_betas, frame)
