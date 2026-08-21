from typing import List, Optional
import warnings

import numpy as np
import pandas as pd
from dataclasses import dataclass

from .cg_probe_table import ProbeTableCol, DesignGroup, CgProbeTable
from .columns import Col
from .config import SieveConfig, validate_sieve_config
from .null_models import NullModels
from .site_statistics import _add_flags, _add_statistics
from .threshold import _find_threshold
from .validation import raise_validation_errors


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


def validate_betas_frame(cg_by_sample: pd.DataFrame) -> None:
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

    raise_validation_errors("Invalid beta-value DataFrame", errors)


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


def sieve_betas(
    cg_by_sample: pd.DataFrame,
    config: SieveConfig,
) -> SieveResults:
    """Run betaSieve on an in-memory IlmnID-by-sample beta-value matrix."""
    validate_sieve_config(config)
    validate_betas_frame(cg_by_sample)
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
    "SieveResults",
    "sieve_betas",
    "validate_betas_frame",
]
