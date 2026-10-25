import pandas as pd

from betasieve.cg_probe_table import (
    CgProbeTable,
    DesignGroup,
    ProbeTableCol,
)


def test_parse_from_probe_ids_drops_invalid_ids_and_preserves_values() -> None:
    probe_ids = pd.Series(
        ["cg1_TC11", "invalid", "cg2_BC210"],
        index=[10, 11, 12],
    )

    result = CgProbeTable.parse_from_probe_ids(probe_ids)

    assert result.to_dict("records") == [
        {
            "probe": "cg1_TC11",
            "site_id": "cg1",
            "design_id": "cg1_TC1",
            "design_type": "TC1",
            "replicate_id": "1",
        },
        {
            "probe": "cg2_BC210",
            "site_id": "cg2",
            "design_id": "cg2_BC2",
            "design_type": "BC2",
            "replicate_id": "10",
        },
    ]


def test_parse_from_probe_ids_returns_empty_typed_frame_for_no_matches() -> None:
    result = CgProbeTable.parse_from_probe_ids(pd.Series(["ch1", "rs2"]))

    assert result.empty
    assert result.columns.tolist() == [
        column.value
        for column in (
            ProbeTableCol.PROBE,
            ProbeTableCol.SITE_ID,
            ProbeTableCol.DESIGN_ID,
            ProbeTableCol.DESIGN_TYPE,
            ProbeTableCol.REPLICATE_ID,
        )
    ]


def test_from_probe_ids_classifies_design_groups_and_exact_replicates() -> None:
    probe_ids = pd.Series(
        [
            "cg1_TC11",
            "cg1_TC21",
            "cg2_TC11",
            "cg2_BC11",
            "cg3_TC11",
            "cg3_BC11",
            "cg3_TO11",
            "cg4_TC11",
            "cg4_TC21",
            "cg4_BC11",
            "cg4_BC21",
            "cg5_TC11",
            "cg5_TC12",
            "cg6_TC11",
        ]
    )

    result = CgProbeTable.from_probe_ids(probe_ids)

    assert set(result.loc[result["site_id"] == "cg1", "group"]) == {
        DesignGroup.PAIR_TYPE.value
    }
    assert set(result.loc[result["site_id"] == "cg2", "group"]) == {
        DesignGroup.PAIR_DESIGN.value
    }
    assert set(result.loc[result["site_id"] == "cg3", "group"]) == {
        DesignGroup.TRIPLET.value
    }
    assert set(result.loc[result["site_id"] == "cg4", "group"]) == {
        DesignGroup.QUADRUPLET.value
    }
    assert result.loc[["cg5_TC11", "cg5_TC12"], "replicate"].tolist() == [
        DesignGroup.EXACT_REPLICATES.value,
        DesignGroup.EXACT_REPLICATES.value,
    ]
    assert pd.isna(result.loc["cg6_TC11", "group"])
    assert pd.isna(result.loc["cg6_TC11", "replicate"])
    assert result.index.name == ProbeTableCol.PROBE.value
