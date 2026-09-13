from enum import Enum
from typing import List

import pandas as pd

from .cg_probe_table import ProbeTableCol


class Col(str, Enum):
    SITE = "site"
    SAMPLE = "sample"
    DIFF = "diff"
    GROUP = ProbeTableCol.GROUP_COL.value
    N = "n"
    THRESHOLD = "threshold"
    CONFIDENCE = "confidence"
    ABOVE = "above_threshold"
    P_HAT = "p_hat"
    P0 = "p0"
    Z = "z"
    Z_OBS = "Z_obs"
    CI_LOWER = "ci_lower"
    CI_UPPER = "ci_upper"
    P_VALUE = "p_value"
    P_ADJUSTED = "p_adjusted"
    P_FLAG = "p_flag"
    P_ADJ_FLAG = "p_adjusted_flag"
    CI_FLAG = "ci_flag"
    PCT_CI_FLAGGED = "pct_ci_flagged"
    PCT_P_ADJ_FLAGGED = "pct_p_adj_flagged"
    P_BETA = "p_beta"
    P_BETA_FLAG = "p_beta_flagged"
    P_BETA_ADJUSTED = "p_beta_adj"
    P_BETA_ADJ_FLAG = "p_beta_adj_flagged"
    PCT_BETA_ADJ_FLAGGED = "pct_p_beta_adj_flagged"
    P_EMPIR = "p_empir"
    P_EMPIR_FLAG = "p_empir_flagged"
    P_EMPIR_ADJUSTED = "p_empir_adj"
    P_EMPIR_ADJ_FLAG = "p_empir_adj_flagged"
    PCT_EMPIR_ADJ_FLAGGED = "pct_p_empir_adj_flagged"


_STAT_META_COLUMNS = frozenset(
    {
        Col.GROUP,
        Col.N,
        Col.THRESHOLD,
        Col.CONFIDENCE,
        Col.ABOVE,
        Col.P_HAT,
        Col.P0,
        Col.Z,
        Col.Z_OBS,
        Col.CI_LOWER,
        Col.CI_UPPER,
        Col.P_VALUE,
        Col.P_ADJUSTED,
        Col.P_FLAG,
        Col.P_ADJ_FLAG,
        Col.CI_FLAG,
        Col.P_BETA,
        Col.P_BETA_FLAG,
        Col.P_BETA_ADJUSTED,
        Col.P_BETA_ADJ_FLAG,
        Col.PCT_BETA_ADJ_FLAGGED,
        Col.P_EMPIR,
        Col.P_EMPIR_FLAG,
        Col.P_EMPIR_ADJUSTED,
        Col.P_EMPIR_ADJ_FLAG,
        Col.PCT_EMPIR_ADJ_FLAGGED,
    }
)


def _diff_value_columns(frame: pd.DataFrame) -> List[str]:
    return [
        col
        for col in frame.columns
        if col not in _STAT_META_COLUMNS and pd.api.types.is_numeric_dtype(frame[col])
    ]


__all__ = ["Col"]
