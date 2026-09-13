from typing import List

from betasieve.columns import Col
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
            "<strong>betaSieve</strong> writes up to three CSV files. "
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
        return "min_max_difference_{threshold}.csv"

    @property
    def description(self) -> str:
        return (
            "Main results file with one row per CpG-site. "
            f"The primary flag column is <code>{Col.P_EMPIR_ADJ_FLAG.value}</code>. "
            "Sites where the entry is <code>True</code> are considered discordant."
        )

    def _figures(self):
        groups = ", ".join(f'"{g.value}"' for g in DesignGroup)
        rows = [
            (
                Col.SITE.value,
                "CpG-site identifier (e.g. cg12345678). Forms the row index.",
            ),
            (Col.GROUP.value, f"Replicate group. One of: {groups}."),
            (Col.N.value, "Number of samples."),
            (Col.THRESHOLD.value, "Applied β-values max-min difference threshold 𝑡."),
            (Col.CONFIDENCE.value, "Confidence level (1 - 𝛼)."),
            (
                Col.ABOVE.value,
                "Number of samples with β-values max-min difference > 𝑡.",
            ),
            (Col.P_HAT.value, "Observed exceedance rate: above_threshold / 𝑛."),
            (
                Col.P0.value,
                "Background exceedance rate estimated from exact-replicate sites.",
            ),
            (Col.Z.value, "Normal quantile Φ⁻¹(confidence)."),
            (Col.Z_OBS.value, "Observed z-score: (𝑝̂ - 𝑝₀) / √(𝑝₀(1-𝑝₀)/𝑛)."),
            (Col.CI_LOWER.value, "Lower bound of Wilson score CI for 𝑝̂."),
            (Col.CI_UPPER.value, "Upper bound of Wilson score CI for 𝑝̂."),
            (Col.P_VALUE.value, "One-sided z-test p-value (reference only)."),
            (Col.P_ADJUSTED.value, "FDR-adjusted z-test p-value (reference only)."),
            (Col.P_FLAG.value, "True if p_value < 𝛼 (reference only)."),
            (Col.P_ADJ_FLAG.value, "True if p_adjusted < 𝛼 (reference only)."),
            (Col.CI_FLAG.value, "True if ci_lower > 𝑝₀ (reference only)."),
            (
                Col.P_EMPIR.value,
                "Empirical upper-tail p-value against the exact-replicate reference distribution.",
            ),
            (Col.P_EMPIR_FLAG.value, "True if p_empir < 𝛼."),
            (
                Col.P_EMPIR_ADJUSTED.value,
                "Multiple-testing-adjusted empirical p-value (primary criterion).",
            ),
            (
                Col.P_EMPIR_ADJ_FLAG.value,
                "True if p_empir_adj < 𝛼 — site is discordant.",
            ),
            (
                "<sample_id>",
                "Per-sample β-values max-min difference. One column per sample.",
            ),
        ]
        return [_data_dict_figure(_HEADERS, rows, _COL_WIDTHS, row_height=30)]


class SweepSummaryCsvSubSection(ReportSubSection):
    @property
    def title(self) -> str:
        return "threshold_sweep_summary.csv"

    @property
    def description(self) -> str:
        return (
            "Written only when threshold selection is automatic. "
            "One row per (replicate group, candidate threshold) combination. "
            "Used to inspect how flagging rates change across the swept threshold range."
        )

    def _figures(self):
        rows = [
            (Col.GROUP.value, "Replicate group (design replicate or exact-replicate category)."),
            ("n_sites", "Number of CpG-sites in this replicate group."),
            (Col.THRESHOLD.value, "Candidate threshold 𝑡."),
            (Col.P0.value, "Background exceedance rate at this threshold."),
            ("pct_p_empir_adj_flagged", "% of sites with p_empir_adj < 𝛼."),
            ("pct_p_adj_flagged", "% of sites with p_adjusted < 𝛼 (reference)."),
            ("pct_ci_flagged", "% of sites with ci_lower > 𝑝₀ (reference)."),
        ]
        return [_data_dict_figure(_HEADERS, rows, _COL_WIDTHS)]


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

    def _figures(self):
        rows = [
            ("IlmnID", "Probe ID in EPICv2 format, e.g. cg12345678_TC11."),
        ]
        return [_data_dict_figure(_HEADERS, rows, _COL_WIDTHS)]
