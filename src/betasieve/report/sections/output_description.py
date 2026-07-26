from typing import List

from betasieve.analysis import Col
from betasieve.cg_probe_table import DesignGroup
from betasieve.report.report_section import ReportMainSection, ReportSubSection
from betasieve.report.tables import _data_dict_figure

_HEADERS = ["Column name", "Description"]
_COL_WIDTHS = [0.25, 0.75]


class OutputDescriptionSection(ReportMainSection):
    @property
    def title(self) -> str:
        return "Output Files"

    @property
    def description(self) -> str:
        return (
            "betaSieve writes up to three CSV files. "
            "Each subsection below lists the columns of one file."
        )

    @property
    def subsection_types(self) -> List[type[ReportSubSection]]:
        sections: List[type[ReportSubSection]] = [MinMaxDiffCsvSubSection]
        if self.results.sweep_df is not None:
            sections.append(SweepSummaryCsvSubSection)
        if self.results.candidate_cpgs is not None:
            sections.append(CandidateCpgsCsvSubSection)
        return sections


class MinMaxDiffCsvSubSection(ReportSubSection):
    @property
    def title(self) -> str:
        return f"min_max_difference_{{threshold}}.csv"

    @property
    def description(self) -> str:
        return (
            "Main results file with one row per CpG site. "
            f"The primary flag column is <code>{Col.P_BETA_ADJ_FLAG.value}</code>. "
            "Sites where the entry is <code>True</code> are considered discordant."
        )

    def generate(self) -> None:
        groups = ", ".join(f'"{g.value}"' for g in DesignGroup)
        rows = [
            (
                Col.SITE.value,
                "CpG site identifier (e.g. cg12345678). Forms the row index.",
            ),
            (Col.GROUP.value, f"Design group. One of: {groups}."),
            (Col.N.value, "Number of samples."),
            (Col.THRESHOLD.value, "Applied β-values max–min range threshold t."),
            (Col.CONFIDENCE.value, "Confidence level (1 − α)."),
            (Col.ABOVE.value, "Number of samples with β-values range > t."),
            (Col.P_HAT.value, "Observed exceedance rate: above_threshold / n."),
            (
                Col.P0.value,
                "Background exceedance rate estimated from exact-replicate sites.",
            ),
            (Col.Z.value, "Normal quantile Φ⁻¹(confidence)."),
            (Col.Z_OBS.value, "Observed z-score: (p̂ − p₀) / √(p₀(1−p₀)/n)."),
            (Col.CI_LOWER.value, "Lower bound of Wilson score CI for p̂."),
            (Col.CI_UPPER.value, "Upper bound of Wilson score CI for p̂."),
            (Col.P_VALUE.value, "One-sided z-test p-value (reference only)."),
            (Col.P_ADJUSTED.value, "FDR-adjusted z-test p-value (reference only)."),
            (Col.P_FLAG.value, "True if p_value < α (reference only)."),
            (Col.P_ADJ_FLAG.value, "True if p_adjusted < α (reference only)."),
            (Col.CI_FLAG.value, "True if ci_lower > p₀ (reference only)."),
            (
                Col.P_BETA.value,
                "Permutation p-value against exact-replicate null distribution.",
            ),
            (Col.P_BETA_FLAG.value, "True if p_beta < α."),
            (
                Col.P_BETA_ADJUSTED.value,
                "FDR-adjusted permutation p-value (primary criterion).",
            ),
            (
                Col.P_BETA_ADJ_FLAG.value,
                "True if p_beta_adj < α — site is discordant.",
            ),
            (
                "<sample_id>",
                "Per-sample β-values max–min range. One column per sample.",
            ),
        ]
        self.figures.append(
            _data_dict_figure(_HEADERS, rows, _COL_WIDTHS, row_height=30)
        )


class SweepSummaryCsvSubSection(ReportSubSection):
    @property
    def title(self) -> str:
        return "threshold_sweep_summary.csv"

    @property
    def description(self) -> str:
        return (
            "Written only when threshold selection is automatic. "
            "One row per (design group, candidate threshold) combination. "
            "Used to inspect how flagging rates change across the swept threshold range."
        )

    def generate(self) -> None:
        rows = [
            (Col.GROUP.value, "Design group."),
            ("n_sites", "Number of CpG sites in this group."),
            (Col.THRESHOLD.value, "Candidate threshold t."),
            (Col.P0.value, "Background exceedance rate at this threshold."),
            ("pct_p_beta_adj_flagged", "% of sites with p_beta_adj < α."),
            ("pct_p_adj_flagged", "% of sites with p_adjusted < α (reference)."),
            ("pct_ci_flagged", "% of sites with ci_lower > p₀ (reference)."),
        ]
        self.figures.append(_data_dict_figure(_HEADERS, rows, _COL_WIDTHS))


class CandidateCpgsCsvSubSection(ReportSubSection):
    @property
    def title(self) -> str:
        return "candidate_cpgs.csv"

    @property
    def description(self) -> str:
        return (
            "Single-column file listing all probe IDs at sites flagged as discordant. "
            "All probes at a flagged site are included, regardless of design type."
        )

    def generate(self) -> None:
        rows = [
            ("IlmnID", "Probe ID in EPICv2 format, e.g. cg12345678_TC11."),
        ]
        self.figures.append(_data_dict_figure(_HEADERS, rows, _COL_WIDTHS))
