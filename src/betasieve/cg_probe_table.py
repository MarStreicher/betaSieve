from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Optional

import pandas as pd

from epicv2io import CgProbeId
from plotly.graph_objs import Frame


class ProbeTableCol(str, Enum):
    PROBE = "probe"
    SITE_ID = "site_id"
    DESIGN_ID = "design_id"
    DESIGN_TYPE = "design_type"
    REPLICATE_ID = "replicate_id"
    GROUP_COL = "group"
    EXACT_REPLICATE_COL = "replicate"


PARSED_PROBE_COLUMNS = [
    ProbeTableCol.SITE_ID,
    ProbeTableCol.DESIGN_ID,
    ProbeTableCol.DESIGN_TYPE,
    ProbeTableCol.REPLICATE_ID,
]


# Pairs where both designs belong to the same chemistry family (TC or BC);
# any other 2-design site is a "design" difference.
_WITHIN_FAMILY_PAIRS = (frozenset({"TC1", "TC2"}), frozenset({"BC1", "BC2"}))


class DesignGroup(str, Enum):
    PAIR_TYPE = "Pair type"
    PAIR_DESIGN = "Pair design"
    TRIPLET = "Triplet"
    QUADRUPLET = "Quadruplet"
    EXACT_REPLICATES = "Exact replicates"


class CgProbeTable:
    @classmethod
    def parse_from_probe_ids(cls, probe_ids: pd.Series) -> pd.DataFrame:
        parsed = probe_ids.map(CgProbeId.parse)
        keep = parsed.notna()

        frame = pd.DataFrame(
            parsed[keep].tolist(),
            columns=[column.value for column in PARSED_PROBE_COLUMNS],
        )
        frame.insert(0, ProbeTableCol.PROBE.value, probe_ids[keep].to_numpy())
        return frame

    @classmethod
    def from_probe_ids(cls, probe_ids: pd.Series) -> pd.DataFrame:
        frame = cls.parse_from_probe_ids(probe_ids)

        def _add_group_col(frame: pd.DataFrame) -> None:
            site_design_sets = frame.groupby(ProbeTableCol.SITE_ID.value)[
                ProbeTableCol.DESIGN_TYPE.value
            ].agg(lambda x: frozenset(str(d) for d in x if d is not None))

            def _classify_site(designs: frozenset) -> Optional[str]:
                n = len(designs)
                if n < 2:
                    return None
                if n == 2:
                    return (
                        DesignGroup.PAIR_TYPE.value
                        if designs in _WITHIN_FAMILY_PAIRS
                        else DesignGroup.PAIR_DESIGN.value
                    )
                return (
                    DesignGroup.TRIPLET.value
                    if n == 3
                    else DesignGroup.QUADRUPLET.value
                )

            site_group = site_design_sets.map(_classify_site).dropna()
            frame[ProbeTableCol.GROUP_COL] = frame[ProbeTableCol.SITE_ID].map(
                site_group
            )

        def _add_exact_replicate_col(frame: pd.DataFrame) -> None:
            design_id_probes = frame.groupby(ProbeTableCol.DESIGN_ID.value)[
                ProbeTableCol.PROBE.value
            ].agg(list)
            replicate_designs = design_id_probes[design_id_probes.map(len) > 1]
            replicate_labels = pd.Series(
                DesignGroup.EXACT_REPLICATES.value,
                index=replicate_designs.index,
            )
            frame[ProbeTableCol.EXACT_REPLICATE_COL] = frame[
                ProbeTableCol.DESIGN_ID
            ].map(replicate_labels)

        _add_group_col(frame)
        _add_exact_replicate_col(frame)

        frame.index = frame[ProbeTableCol.PROBE]
        return frame
