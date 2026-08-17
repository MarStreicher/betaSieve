from enum import Enum
from typing import List, Optional, Tuple
import warnings

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
from dataclasses import dataclass

from .cg_probe_table import ProbeTableCol, DesignGroup, CgProbeTable
from .config import SieveConfig, validate_sieve_config
from .null_models import NullModels


@dataclass
class SieveResults:
    diff_frame: pd.DataFrame
    threshold: float
    statistics_frame: pd.DataFrame
    flagged_frame: pd.DataFrame
    sieved_betas: pd.DataFrame
    candidate_cpgs: pd.Series
    null_models: NullModels
    sweep_df: Optional[pd.DataFrame] = None


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


def _validate_betas_frame(cg_by_sample: pd.DataFrame) -> None:
    errors: List[str] = []

    if not isinstance(cg_by_sample, pd.DataFrame):
        raise TypeError(
            "cg_by_sample must be a pandas.DataFrame, "
            f"got {type(cg_by_sample).__name__}."
        )

    if cg_by_sample.empty:
        errors.append("cg_by_sample must contain at least one probe and one sample.")
    if not cg_by_sample.index.is_unique:
        errors.append("probe IDs in the index must be unique.")
    if not cg_by_sample.columns.is_unique:
        errors.append("sample names in the columns must be unique.")
    if cg_by_sample.index.hasnans:
        errors.append("probe IDs in the index must not be missing.")

    non_numeric = [
        str(column)
        for column in cg_by_sample.columns
        if not pd.api.types.is_numeric_dtype(cg_by_sample[column])
    ]
    if non_numeric:
        errors.append(
            "all sample columns must be numeric; non-numeric columns: "
            + ", ".join(non_numeric)
            + "."
        )
    elif not cg_by_sample.empty:
        if cg_by_sample.isna().to_numpy().any():
            errors.append("beta values must not contain missing or infinite values.")
        else:
            values = cg_by_sample.to_numpy(dtype=float)
            if not np.isfinite(values).all():
                errors.append(
                    "beta values must not contain missing or infinite values."
                )
            elif ((values < 0.0) | (values > 1.0)).any():
                errors.append("beta values must be between 0 and 1 (inclusive).")

    if len(cg_by_sample.index) > 0 and not cg_by_sample.index.hasnans:
        probe_ids = pd.Series(cg_by_sample.index.astype(str))
        try:
            parsed = CgProbeTable.parse_from_probe_ids(probe_ids)
        except (TypeError, ValueError):
            parsed = pd.DataFrame()
        if len(parsed) != len(probe_ids):
            errors.append(
                "every index value must be a valid EPICv2 IlmnID "
                "(for example, cg00000001_TC11)."
            )

    if errors:
        message = "Invalid beta-value DataFrame:\n" + "\n".join(
            f"  • {error}" for error in errors
        )
        raise ValueError(message)


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


def _create_cpg_list(
    frame: pd.DataFrame,
    cg_by_sample: pd.DataFrame,
) -> pd.Series:
    """
    Return all probes in ``cg_by_sample`` whose site belongs to a flagged site.
    """
    flagged_sites = set(frame.index[frame[Col.P_EMPIR_ADJ_FLAG]])
    print(f"Number of flagged sites: {len(flagged_sites)}")

    probe_index = pd.Series(cg_by_sample.index.astype(str))
    site_prefix = probe_index.str.extract(r"^(cg\d+)_", expand=False)
    final_cpgs = probe_index.loc[site_prefix.isin(flagged_sites)].reset_index(drop=True)
    final_cpgs.name = "IlmnID"

    if final_cpgs.empty:
        warnings.warn("No candidate CpGs found for flagged sites.")
        return final_cpgs

    print(f"Number of candidate probe instances: {len(final_cpgs)}")
    return final_cpgs


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


def _collect_max_min_differences(
    cg_by_sample: pd.DataFrame, cg_by_group: pd.DataFrame
) -> pd.DataFrame:
    sample_cols = cg_by_sample.columns.tolist()
    merged = cg_by_sample.join(cg_by_group, how="inner")

    diffs_groups = merged.groupby([ProbeTableCol.SITE_ID, ProbeTableCol.GROUP_COL])[
        sample_cols
    ].agg(lambda frame: frame.max() - frame.min())

    diffs_exact_replicates = merged.groupby(
        [ProbeTableCol.SITE_ID, ProbeTableCol.EXACT_REPLICATE_COL]
    )[sample_cols].agg(lambda frame: frame.max() - frame.min())

    index_names = [ProbeTableCol.SITE_ID.value, Col.GROUP.value]
    diffs_groups.index = diffs_groups.index.set_names(index_names)
    diffs_exact_replicates.index = diffs_exact_replicates.index.set_names(index_names)

    diffs_groups = diffs_groups.loc[
        diffs_groups.index.get_level_values(Col.GROUP).notna()
    ]
    diffs_exact_replicates = diffs_exact_replicates.loc[
        diffs_exact_replicates.index.get_level_values(Col.GROUP).notna()
    ]

    result = pd.concat([diffs_groups, diffs_exact_replicates], axis=0)
    return result.reset_index(level=Col.GROUP)


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


def sieve_betas(
    cg_by_sample: pd.DataFrame,
    config: SieveConfig,
) -> SieveResults:
    """Run betaSieve on an in-memory IlmnID-by-sample beta-value matrix."""
    validate_sieve_config(config)
    _validate_betas_frame(cg_by_sample)
    sweep_df: Optional[pd.DataFrame] = None

    print("Building (design) groups based on the IlmnID...")
    cg_by_group = CgProbeTable.from_probe_ids(pd.Series(cg_by_sample.index.tolist()))

    print("Computing max-min ranges per site per sample...")
    diff_frame = _collect_max_min_differences(cg_by_sample, cg_by_group)
    groups = set(diff_frame[Col.GROUP])
    if DesignGroup.EXACT_REPLICATES.value not in groups:
        raise ValueError(
            "The beta-value DataFrame must contain at least one exact-replicate "
            "probe group to estimate the empirical background rate."
        )
    if groups == {DesignGroup.EXACT_REPLICATES.value}:
        raise ValueError(
            "The beta-value DataFrame must contain at least one duplicate-design "
            "probe group to test."
        )

    if config.threshold is None:
        assert config.threshold_min is not None
        assert config.threshold_max is not None
        assert config.threshold_step is not None

        threshold, sweep_df = _find_threshold(
            diff_frame,
            threshold_min=config.threshold_min,
            threshold_max=config.threshold_max,
            threshold_step=config.threshold_step,
            fdr=config.fdr,
            confidence=config.confidence,
            target_p0=config.target_p0,
        )
    else:
        threshold = config.threshold

    print(f"Computing statistics at threshold {threshold}...")
    statistics_frame, null_models = _add_statistics(
        diff_frame, threshold, config.fdr, config.confidence
    )

    print("Adding flagged columns...")
    flagged_frame = _add_flags(statistics_frame)

    print("Creating list CpG pd.Series ...")
    cpg_serie = _create_cpg_list(flagged_frame, cg_by_sample)

    sieved_frame = cg_by_sample.loc[~cg_by_sample.index.isin(cpg_serie)]

    return SieveResults(
        diff_frame=diff_frame,
        threshold=threshold,
        statistics_frame=statistics_frame,
        flagged_frame=flagged_frame,
        sweep_df=sweep_df,
        candidate_cpgs=cpg_serie,
        null_models=null_models,
        sieved_betas=sieved_frame,
    )


__all__ = [
    "Col",
    "SieveResults",
    "sieve_betas",
]
