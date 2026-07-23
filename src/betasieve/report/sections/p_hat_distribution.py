from __future__ import annotations

from typing import List

import plotly.graph_objects as go
from plotly.graph_objects import Figure

from betasieve.analysis import Col
from betasieve.cg_probe_table import DesignGroup
from betasieve.report.domain.mappings import (
    BS_DARK,
    BS_GREEN,
    DESIGN_GROUP_COLORS,
    REPORT_FONT_FAMILY,
)
from betasieve.report.plots import _layout_figure
from betasieve.report.report_section import ReportMainSection, ReportSubSection


class PhatDistributionSection(ReportMainSection):
    @property
    def title(self) -> str:
        return "Observed Exceedance Rate"

    @property
    def description(self) -> str:
        return "The observed exceedance rate p̂ is the proportion of samples for which the β-value max–min range at a given CpG site exceeds the calibrated threshold t. "

    @property
    def subsection_types(self) -> List[type[ReportSubSection]]:
        return [PhatHistogram, PhatBoxplot]


class PhatHistogram(ReportSubSection):
    @property
    def title(self) -> str:
        return "p̂ Distribution by Design Group"

    @property
    def description(self) -> str:
        return (
            "Distribution of the observed exceedance rate p̂ across all CpG sites, "
            "stratified by design group. Please, click on a design group in the legend to include or exclude it."
        )

    def _plot(self) -> Figure:
        df = self.results.flagged_frame
        fig = go.Figure()

        for group in DesignGroup:
            mask = df[Col.GROUP] == group
            values = df.loc[mask, [Col.P_HAT]].dropna().to_numpy().ravel()
            color = DESIGN_GROUP_COLORS.get(group, BS_GREEN)
            r, g, b = int(color[1:3], 16), int(color[3:5], 16), int(color[5:7], 16)
            fill = f"rgba({r},{g},{b},0.4)"
            fig.add_trace(
                go.Histogram(
                    x=values,
                    name=group.value,
                    marker=dict(color=fill, line=dict(color=color, width=2)),
                    opacity=0.7,
                    nbinsx=40,
                )
            )
        fig.update_layout(barmode="overlay")
        return _layout_figure(
            fig,
            title="Distribution of observed exceedance rates p̂ by design group",
            x_title="p̂ (observed exceedance rate)",
            y_title="Number of CpG sites",
            height=400,
        )

    def generate(self) -> None:
        self.figures.append(self._plot())


class PhatBoxplot(ReportSubSection):
    @property
    def title(self) -> str:
        return "p̂ distribution by design group"

    @property
    def description(self) -> str:
        return (
            "Box plots of the observed exceedance rate p̂ per design group. "
            "The dashed horizontal line marks the empirical background exceedance "
            "rate p₀, estimated from exact-replicate sites. "
        )

    def _plot(self) -> Figure:
        flagged = self.results.flagged_frame
        p0 = float(flagged[Col.P0].iloc[0])

        fig = go.Figure()

        for group in DesignGroup:
            mask = flagged[Col.GROUP] == group
            values = flagged.loc[mask, Col.P_HAT].dropna().values
            if len(values) == 0:
                continue
            color = DESIGN_GROUP_COLORS.get(group, BS_GREEN)
            r = int(color[1:3], 16)
            g_val = int(color[3:5], 16)
            b = int(color[5:7], 16)
            fill = f"rgba({r},{g_val},{b},0.15)"
            fig.add_trace(
                go.Box(
                    y=values,
                    name=group.value,
                    boxpoints="outliers",
                    marker=dict(color=color, size=3, opacity=0.5),
                    line=dict(color=color),
                    fillcolor=fill,
                )
            )

        fig.add_hline(
            y=p0,
            line_dash="dash",
            line_color=BS_DARK,
            annotation_text=f"p\u2080 = {p0:.5f}",
            annotation_position="top right",
            annotation_font=dict(color=BS_DARK, family=REPORT_FONT_FAMILY),
        )

        return _layout_figure(
            fig,
            title="Observed exceedance rate p̂ by design group",
            x_title="Design group",
            y_title="p̂ (observed exceedance rate)",
            height=400,
        )

    def generate(self) -> None:
        self.figures.append(self._plot())
