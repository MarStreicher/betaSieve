from typing import List

import plotly.graph_objects as go
from plotly.graph_objects import Figure

from betasieve.analysis import Col
from betasieve.cg_probe_table import DesignGroup
from betasieve.report.domain.mappings import (
    BS_DARK,
    BS_GREEN,
    BS_LIGHTEST,
    BS_MID,
    DESIGN_GROUP_COLORS,
    REPORT_FONT_FAMILY,
)
from betasieve.report.plots import _layout_figure
from betasieve.report.report_section import ReportMainSection, ReportSubSection
from betasieve.report.tables import _summary_table_figure


def _add_selected_threshold_vline(fig: Figure, threshold: float) -> None:
    fig.add_vline(
        x=threshold,
        line_dash="dot",
        line_color=BS_DARK,
        annotation_text=f" {threshold:g}",
        annotation_position="top right",
        annotation_font_color=BS_DARK,
    )


def _fixed_threshold_placeholder(threshold: float) -> Figure:
    fig = go.Figure()
    fig.add_annotation(
        text=(
            f"No threshold sweep was performed — "
            f"a fixed threshold of {threshold:g} was used."
        ),
        xref="paper",
        yref="paper",
        x=0.5,
        y=0.5,
        xanchor="center",
        yanchor="middle",
        showarrow=False,
        font=dict(family=REPORT_FONT_FAMILY, size=14, color=BS_MID),
    )
    fig.update_layout(
        xaxis=dict(visible=False, fixedrange=True),
        yaxis=dict(visible=False, fixedrange=True),
        margin=dict(l=12, r=12, t=4, b=4),
        height=64,
        paper_bgcolor="white",
        plot_bgcolor=BS_LIGHTEST,
    )
    return fig


def _add_group_traces(
    fig: Figure,
    sweep_df,
    *,
    y_column: str,
) -> None:
    for group in DesignGroup:
        subset = sweep_df.loc[sweep_df[Col.GROUP] == group].sort_values(Col.THRESHOLD)
        if subset.empty:
            continue
        color = DESIGN_GROUP_COLORS.get(group, BS_GREEN)
        fig.add_trace(
            go.Scatter(
                x=subset[Col.THRESHOLD],
                y=subset[y_column],
                mode="lines+markers",
                name=group.value,
                line=dict(color=color),
                marker=dict(color=color),
            )
        )


class ThresholdSweepSection(ReportMainSection):
    @property
    def title(self) -> str:
        return "Threshold Sweep"

    @property
    def description(self) -> str:
        if self.results.sweep_df is None:
            return (
                "A user-specified fixed threshold was used for this run, therefore, "
                "threshold sensitivity curves are not available."
            )
        return (
            "The threshold t that separates assay noise from genuine probe disagreement "
            "was selected automatically by sweeping a range of candidate values. "
            "For each candidate t, the empirical background exceedance rate p₀ was "
            "computed from exact-replicate sites. "
            "The smallest t for which p₀ fell at or below the target proportion p₀ was chosen as "
            "the final threshold."
        )

    @property
    def subsection_types(self) -> List[type[ReportSubSection]]:
        if self.results.sweep_df is None:
            return [
                FixedThresholdPlaceholderSubSection,
                ThresholdSummaryTableSubSection,
            ]
        return [
            P0SubSection,
            BetaSection,
            ThresholdSummaryTableSubSection,
        ]


class FixedThresholdPlaceholderSubSection(ReportSubSection):
    @property
    def title(self) -> str:
        return "Fixed threshold"

    @property
    def description(self) -> str:
        return ""

    def _plot(self) -> Figure:
        return _fixed_threshold_placeholder(self.results.threshold)

    def generate(self) -> None:
        self.figures.append(self._plot())


class P0SubSection(ReportSubSection):
    @property
    def title(self) -> str:
        return "Background Exceedance Rate p₀ vs Threshold"

    @property
    def description(self) -> str:
        return (
            "The empirical background exceedance rate p₀ is the fraction of "
            "(exact-replicate site, sample) pairs for which the β-values max–min range "
            "exceeds the candidate threshold t. "
        )

    def _plot(self) -> Figure:
        sweep_df = self.results.sweep_df
        if sweep_df is None:
            return _fixed_threshold_placeholder(self.results.threshold)

        subset = sweep_df.drop_duplicates(subset=Col.THRESHOLD).sort_values(
            Col.THRESHOLD
        )
        threshold = self.results.threshold

        fig = go.Figure(
            go.Scatter(
                x=subset[Col.THRESHOLD],
                y=subset[Col.P0],
                mode="lines+markers",
                name="p₀ (exact replicates)",
                line=dict(color=BS_GREEN),
                marker=dict(color=BS_GREEN),
            )
        )
        _add_selected_threshold_vline(fig, threshold)
        return _layout_figure(
            fig,
            title="Empirical background exceedance rate p₀ vs threshold t",
            x_title="Threshold t",
            y_title="p₀ (background exceedance rate)",
            height=360,
            show_legend=False,
        )

    def generate(self) -> None:
        self.figures.append(self._plot())


class BetaSection(ReportSubSection):
    @property
    def title(self) -> str:
        return "Adjusted Empirical Flag Rate vs Threshold"

    @property
    def description(self) -> str:
        return (
            "Percentage of CpG sites per design group for which the FDR-adjusted "
            "empirical upper-tail p-value falls below α, plotted over the "
            "threshold candidates. "
            "The dotted vertical line marks the selected threshold."
        )

    def _plot(self) -> Figure:
        sweep_df = self.results.sweep_df
        if sweep_df is None:
            return _fixed_threshold_placeholder(self.results.threshold)

        fig = go.Figure()
        _add_group_traces(fig, sweep_df, y_column=Col.PCT_EMPIR_ADJ_FLAGGED)
        _add_selected_threshold_vline(fig, self.results.threshold)
        return _layout_figure(
            fig,
            title="Adjusted empirical flag rate vs threshold t",
            x_title="Threshold t",
            y_title="Sites with p-value < α (%)",
            height=400,
        )

    def generate(self) -> None:
        self.figures.append(self._plot())


class ThresholdSummaryTableSubSection(ReportSubSection):
    @property
    def title(self) -> str:
        return "Selected Threshold Summary"

    @property
    def description(self) -> str:
        return (
            "The final threshold t and corresponding background exceedance rate p₀ "
            "used for all downstream statistical tests and figures."
        )

    def _summary_table(self) -> Figure:
        flagged = self.results.flagged_frame
        rows = [
            (
                "Threshold t",
                round(float(flagged[Col.THRESHOLD].iloc[0]), 4),
            ),
            ("Background exceedance rate p₀", round(float(flagged[Col.P0].iloc[0]), 4)),
        ]
        return _summary_table_figure(rows)

    def generate(self) -> None:
        self.figures.append(self._summary_table())
