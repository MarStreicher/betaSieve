from typing import List, Tuple

import numpy as np
import pandas as pd

from .cg_probe_table import DesignGroup
from .columns import Col
from .site_statistics import _add_flags, _add_statistics


def _flag_rates_by_group(frame: pd.DataFrame) -> pd.DataFrame:
    return (
        frame.groupby(Col.GROUP, sort=False)
        .agg(
            n_sites=(Col.CI_FLAG, "count"),
            pct_ci_flagged=(Col.CI_FLAG, lambda s: 100.0 * s.mean()),
            pct_p_adj_flagged=(Col.P_ADJ_FLAG, lambda s: 100.0 * s.mean()),
            pct_p_empir_adj_flagged=(
                Col.P_EMPIR_ADJ_FLAG,
                lambda s: 100.0 * s.mean(),
            ),
        )
        .reset_index()
    )


def _sweep_thresholds(
    diff_frame: pd.DataFrame,
    *,
    threshold_min: float,
    threshold_max: float,
    threshold_step: float,
    fdr: str,
    confidence: float,
) -> pd.DataFrame:
    rows: List[pd.DataFrame] = []
    thresholds = np.arange(
        threshold_min, threshold_max + threshold_step / 2, threshold_step
    )

    for threshold in thresholds:
        statistics_frame, _ = _add_statistics(
            diff_frame, float(threshold), fdr, confidence
        )
        P0_value = statistics_frame[Col.P0].unique()[0]
        flagged_frame = _add_flags(statistics_frame)
        rates = _flag_rates_by_group(flagged_frame)
        rates[Col.THRESHOLD] = float(threshold)
        rates[Col.P0] = float(P0_value)
        rows.append(rates)

    return pd.concat(rows, ignore_index=True)


def _select_threshold_for_replicates(
    sweep_df: pd.DataFrame,
    target_p0: float,
) -> float:
    replicate_rates = (
        sweep_df.loc[sweep_df[Col.GROUP] == DesignGroup.EXACT_REPLICATES.value]
        .sort_values(Col.THRESHOLD)
        .reset_index(drop=True)
    )
    if replicate_rates.empty:
        raise ValueError("No sweep results found!")

    qualifying = replicate_rates.loc[
        replicate_rates[Col.P0] <= target_p0, Col.THRESHOLD
    ]
    if len(qualifying) > 0:
        return float(qualifying.iloc[0])

    best_idx = replicate_rates[Col.P0].idxmin()
    chosen = float(replicate_rates.loc[best_idx, Col.THRESHOLD])
    best_p0 = float(replicate_rates.loc[best_idx, Col.P0])
    print(
        f"Warning: target_p0={target_p0} was not reached in the sweep "
        f"(best p0={best_p0:.4f} at threshold={chosen}). "
        "Consider increasing target threshold."
    )
    return chosen


def _find_threshold(
    diff_frame: pd.DataFrame,
    *,
    threshold_min: float,
    threshold_max: float,
    threshold_step: float,
    fdr: str,
    confidence: float,
    target_p0: float,
) -> Tuple[float, pd.DataFrame]:
    print(
        f"Sweeping thresholds from {threshold_min} to {threshold_max} "
        f"(step {threshold_step})..."
    )
    sweep_df = _sweep_thresholds(
        diff_frame,
        threshold_min=threshold_min,
        threshold_max=threshold_max,
        threshold_step=threshold_step,
        fdr=fdr,
        confidence=confidence,
    )

    chosen = _select_threshold_for_replicates(sweep_df, target_p0)
    at_chosen = sweep_df.loc[
        (sweep_df[Col.GROUP] == DesignGroup.EXACT_REPLICATES.value)
        & (sweep_df[Col.THRESHOLD] == chosen)
    ]
    ci_pct = (
        float(at_chosen[Col.PCT_CI_FLAGGED].iloc[0]) if len(at_chosen) else float("nan")
    )
    empir_pct = (
        float(at_chosen[Col.PCT_EMPIR_ADJ_FLAGGED].iloc[0])
        if len(at_chosen)
        else float("nan")
    )
    print(
        f"Selected threshold {chosen} "
        f"(exact replicates: CI flag rate {ci_pct:.2g}%, "
        f"P_EMPIR_ADJ flag rate {empir_pct:.2g}%)"
    )
    return chosen, sweep_df
