from typing import Optional, Tuple

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from .cg_probe_table import DesignGroup
from .columns import Col, _diff_value_columns
from .null_models import NullModels


def _add_statistics(
    frame: pd.DataFrame,
    threshold: float,
    fdr: str,
    confidence: float,
) -> Tuple[pd.DataFrame, NullModels]:

    def calculate_empirical_background_rate(
        frame: pd.DataFrame, threshold: float
    ) -> float:
        value_cols = _diff_value_columns(frame)
        replicate_diffs = frame.loc[
            frame[Col.GROUP] == DesignGroup.EXACT_REPLICATES.value, value_cols
        ]
        n_above_threshold = (replicate_diffs > threshold).sum().sum()
        return float(n_above_threshold / replicate_diffs.size)

    def apply_fdr(
        column: str, confidence: float, method: str, mask: Optional[pd.Series] = None
    ) -> pd.Series:
        if mask is None:
            mask = pd.Series(True, index=samples_frame.index)

        p_raw = samples_frame.loc[mask, column]
        _, p_adj, _, _ = multipletests(
            p_raw,
            alpha=1 - confidence,
            method=method,
        )
        result = samples_frame.loc[:, column].copy()
        result[mask] = p_adj
        return result

    value_cols = _diff_value_columns(frame)
    diffs = frame[value_cols]
    n_samples = len(value_cols)

    samples_frame = pd.DataFrame(index=frame.index)
    samples_frame[Col.N] = n_samples
    samples_frame[Col.THRESHOLD] = threshold
    samples_frame[Col.CONFIDENCE] = confidence
    samples_frame[Col.ABOVE] = (diffs > threshold).sum(axis=1)
    samples_frame[Col.Z] = stats.norm.ppf(confidence)
    samples_frame[Col.P_HAT] = samples_frame[Col.ABOVE] / n_samples

    samples_frame[Col.P0] = calculate_empirical_background_rate(frame, threshold)

    # Wilson score confidence interval
    samples_frame[Col.Z_OBS] = (
        samples_frame[Col.P_HAT] - samples_frame[Col.P0]
    ) / np.sqrt(
        samples_frame[Col.P0] * (1 - samples_frame[Col.P0]) / samples_frame[Col.N]
    )

    _z = samples_frame[Col.Z]
    _p = samples_frame[Col.P_HAT]
    _n = samples_frame[Col.N]
    _z2 = _z**2
    _center = _p + _z2 / (2 * _n)
    _margin = _z * np.sqrt(_p * (1 - _p) / _n + _z2 / (4 * _n**2))
    _denom = 1 + _z2 / _n
    samples_frame[Col.CI_LOWER] = (_center - _margin) / _denom
    samples_frame[Col.CI_UPPER] = (_center + _margin) / _denom

    samples_frame[Col.GROUP] = frame[Col.GROUP]
    test_mask = samples_frame[Col.GROUP] != DesignGroup.EXACT_REPLICATES.value

    # empirical, beta-binomial, binomial
    null_models = NullModels.fit(
        n=n_samples,
        p0=float(samples_frame[Col.P0].iloc[0]),
        p_hat_reference=samples_frame.loc[~test_mask, Col.P_HAT].to_numpy(),
        threshold=threshold,
        confidence=confidence,
    )

    samples_frame[Col.P_EMPIR] = null_models.empirical.sf(
        samples_frame[Col.P_HAT].to_numpy()
    )
    samples_frame[Col.P_EMPIR_ADJUSTED] = apply_fdr(
        Col.P_EMPIR, confidence, fdr, mask=test_mask
    )

    samples_frame[Col.P_BETA] = null_models.beta_binomial.sf(
        samples_frame[Col.ABOVE].to_numpy()
    )
    samples_frame[Col.P_BETA_ADJUSTED] = apply_fdr(
        Col.P_BETA, confidence, fdr, mask=test_mask
    )

    samples_frame[Col.P_VALUE] = null_models.binomial.sf(
        samples_frame[Col.ABOVE].to_numpy()
    )
    samples_frame[Col.P_ADJUSTED] = apply_fdr(
        Col.P_VALUE, confidence, fdr, mask=test_mask
    )

    return samples_frame, null_models


def _add_flags(frame: pd.DataFrame) -> pd.DataFrame:
    alpha = 1 - frame[Col.CONFIDENCE]
    test_mask = frame[Col.GROUP] != DesignGroup.EXACT_REPLICATES.value
    frame[Col.P_FLAG] = frame[Col.P_VALUE] < alpha
    frame[Col.P_ADJ_FLAG] = frame[Col.P_ADJUSTED] < alpha
    frame[Col.P_EMPIR_FLAG] = (frame[Col.P_EMPIR] < alpha) & test_mask
    frame[Col.P_EMPIR_ADJ_FLAG] = (frame[Col.P_EMPIR_ADJUSTED] < alpha) & test_mask
    frame[Col.P_BETA_FLAG] = (frame[Col.P_BETA] < alpha) & test_mask
    frame[Col.P_BETA_ADJ_FLAG] = (frame[Col.P_BETA_ADJUSTED] < alpha) & test_mask
    frame[Col.CI_FLAG] = frame[Col.CI_LOWER] > frame[Col.P0]
    return frame
