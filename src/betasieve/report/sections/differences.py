from typing import List, Optional

import numpy as np
import plotly.graph_objects as go
from plotly.graph_objects import Figure

from betasieve.columns import Col
from betasieve.cg_probe_table import DesignGroup
from betasieve.report.domain.mappings import (
    BS_GREEN,
    BS_HEATMAP_COLORSCALE,
    DESIGN_GROUP_COLORS,
    REPORT_FONT_FAMILY,
)
from betasieve.report.plots import _layout_figure
from betasieve.report.report_section import ReportMainSection, ReportSubSection


class DifferencesSection(ReportMainSection):
    @property
    def title(self) -> str:
        return "Probe Range Analysis"

    @property
    def description(self) -> str:
        return (
            "For each duplicate CpG site and each sample, the max–min range of "
            "β-values across all probes in the group was computed. "
        )

    @property
    def subsection_types(self) -> List[type[ReportSubSection]]:
        heatmap_sections = [
            PairTypeHeatmapSubSection,
            PairDesignHeatmapSubSection,
            TripletHeatmapSubSection,
            QuadrupletHeatmapSubSection,
        ]
        available_heatmaps = [
            section
            for section in heatmap_sections
            if (self.results.diff_frame[Col.GROUP] == section.GROUP_KEY).any()
        ]
        return [
            OneDifferencesPercentageHistogram,
            DifferencesPercentageHistogram,
            DifferencesHistogram,
            DifferencesBoxplot,
            *available_heatmaps,
        ]


# TODO: Refine and fix
class OneDifferencesPercentageHistogram(ReportSubSection):
    @property
    def title(self) -> str:
        return "Relative β-values Range Distribution"

    @property
    def description(self) -> str:
        return (
            "Relative distribution of β-value max–min ranges, comparing exact "
            "replicates with all non-replicate design groups pooled together. "
            "Both distributions use the same bin edges."
        )

    def _plot(self) -> Figure:
        df = self.results.diff_frame
        fig = go.Figure()

        value_cols = [column for column in df.columns if column != Col.GROUP]
        x_max = 0.3
        shared_bins = dict(start=0, end=x_max, size=x_max / 100)

        exact_mask = df[Col.GROUP] == DesignGroup.EXACT_REPLICATES
        distributions = (
            ("Exact replicates", df.loc[exact_mask, value_cols], BS_GREEN),
            ("Other design groups", df.loc[~exact_mask, value_cols], "#F97316"),
        )

        for name, subset, color in distributions:
            values = subset.to_numpy().ravel()
            values = values[np.isfinite(values)]
            r, g, b = int(color[1:3], 16), int(color[3:5], 16), int(color[5:7], 16)
            fill = f"rgba({r},{g},{b},0.4)"
            fig.add_trace(
                go.Histogram(
                    x=values,
                    histnorm="percent",
                    name=name,
                    marker=dict(color=fill, line=dict(color=color, width=2)),
                    opacity=0.7,
                    xbins=shared_bins,
                )
            )
        fig.update_layout(barmode="overlay", bargap=0.05, width=700, autosize=False)
        fig.update_xaxes(range=[0, x_max])
        return _layout_figure(
            fig,
            title="",
            x_title="β-values max–min range",
            y_title="Percentage of observations",
            height=400,
        )

    def generate(self) -> None:
        self.figures.append(self._plot())


# TODO: Refine and fix (Also the test)
class DifferencesPercentageHistogram(ReportSubSection):
    @property
    def title(self) -> str:
        return "Relative β-values Range Distribution"

    @property
    def description(self) -> str:
        return "TODO"

    def _plot(self) -> Figure:
        df = self.results.diff_frame
        fig = go.Figure()

        for group in DesignGroup:
            sub_df = df[df[Col.GROUP] == group]
            values = (
                sub_df.drop(columns=[Col.GROUP], errors="ignore").to_numpy().ravel()
            )
            color = DESIGN_GROUP_COLORS.get(group, BS_GREEN)
            r, g, b = int(color[1:3], 16), int(color[3:5], 16), int(color[5:7], 16)
            fill = f"rgba({r},{g},{b},0.4)"
            fig.add_trace(
                go.Histogram(
                    x=values,
                    histnorm="percent",
                    name=group.value,
                    marker=dict(color=fill, line=dict(color=color, width=2)),
                    opacity=0.7,
                    nbinsx=40,
                )
            )
        fig.update_layout(barmode="overlay")
        return _layout_figure(
            fig,
            title="Distribution of β-values max–min ranges by design group",
            x_title="β-values max–min range",
            y_title="Number of observations",
            height=400,
        )

    def generate(self) -> None:
        self.figures.append(self._plot())


class DifferencesHistogram(ReportSubSection):
    @property
    def title(self) -> str:
        return "Absolute β-values Range Distribution"

    @property
    def description(self) -> str:
        return (
            "Distribution of the β-values max–min range across all (site, sample) "
            "observations, stratified by design group. Please, click on a design group in the legend to include or exclude it."
        )

    def _plot(self) -> Figure:
        df = self.results.diff_frame
        fig = go.Figure()

        for group in DesignGroup:
            sub_df = df[df[Col.GROUP] == group]
            values = (
                sub_df.drop(columns=[Col.GROUP], errors="ignore").to_numpy().ravel()
            )
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
            title="Distribution of β-values max–min ranges by design group",
            x_title="β-values max–min range",
            y_title="Number of observations",
            height=400,
        )

    def generate(self) -> None:
        self.figures.append(self._plot())


class DifferencesBoxplot(ReportSubSection):
    @property
    def title(self) -> str:
        return "β-values Range Distribution (Box Plots)"

    @property
    def description(self) -> str:
        return (
            "Box plots of the β-values max–min range pooled over all sites and samples "
        )

    def _plot(self) -> Figure:
        diff = self.results.diff_frame
        value_cols = [col for col in diff.columns if col != Col.GROUP]

        fig = go.Figure()

        for group in DesignGroup:
            group_frame = diff[diff[Col.GROUP] == group]
            values = group_frame[value_cols].to_numpy().ravel()
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

        return _layout_figure(
            fig,
            title="β-values max–min range by design group",
            x_title="Design group",
            y_title="β-values max–min range",
            height=400,
        )

    def generate(self) -> None:
        self.figures.append(self._plot())


class _GroupHeatmapSubSection(ReportSubSection):
    GROUP_KEY: DesignGroup = DesignGroup.PAIR_TYPE

    @property
    def title(self) -> str:
        return f"Heatmap - {self.GROUP_KEY.value}"

    @property
    def description(self) -> str:
        return (
            f"Per-site, per-sample β-values max-min ranges for the {self.GROUP_KEY.value} group. "
            "Colour intensity reflects the magnitude of the range "
            "(square-root-transformed for visual clarity). "
            "Rows represent CpG sites and the columns represent "
            "individual samples. "
        )

    def _plot_for_group(
        self,
        group: DesignGroup,
        value_cols: List[str],
        global_vmax: float,
    ) -> Optional[Figure]:
        df = self.results.diff_frame
        subset = df[df[Col.GROUP] == group].sort_index()
        if subset.empty:
            return None

        y_labels = [str(site) for site in subset.index]
        z_raw = subset[value_cols].to_numpy()
        z_values = np.sqrt(z_raw)  # spreads out the small values

        vmax = max(0.0, global_vmax)
        tick_orig = np.linspace(0, vmax, 6) if vmax > 0 else np.array([0.0])
        tickvals = np.sqrt(tick_orig)
        ticktext = [f"{v:.3f}" for v in tick_orig]

        fig = go.Figure(
            go.Heatmap(
                z=z_values,
                x=value_cols,
                y=y_labels,
                zmin=0,
                zmax=float(np.sqrt(vmax)) if vmax > 0 else None,
                colorscale=BS_HEATMAP_COLORSCALE,
                customdata=z_raw,
                colorbar=dict(
                    title="β-values max-min",
                    tickvals=tickvals,
                    ticktext=ticktext,
                ),
                hovertemplate=(
                    "Sample: %{x}<br>"
                    "CpG: %{y}<br>"
                    "Difference: %{customdata:.4f}<extra></extra>"
                ),
            )
        )

        fig.update_layout(
            template="plotly_white",
            title=dict(
                text=f"β-values max-min range heatmap - {group.value}",
                x=0.5,
                xanchor="center",
            ),
            xaxis_title="Sample",
            yaxis_title="CpG site",
            height=max(420, min(1200, 14 * len(y_labels))),
            font=dict(family=REPORT_FONT_FAMILY),
        )
        fig.update_yaxes(autorange="reversed")
        return fig

    def generate(self) -> None:
        df = self.results.diff_frame
        value_cols = [col for col in df.columns if col != Col.GROUP]
        if len(value_cols) == 0:
            return

        non_exact_df = df[df[Col.GROUP] != DesignGroup.EXACT_REPLICATES]
        if non_exact_df.empty:
            return

        global_vmax = float(np.nanmax(non_exact_df[value_cols].to_numpy()))
        fig = self._plot_for_group(self.GROUP_KEY, value_cols, global_vmax)
        if fig is not None:
            self.figures.append(fig)


class PairTypeHeatmapSubSection(_GroupHeatmapSubSection):
    GROUP_KEY = DesignGroup.PAIR_TYPE


class PairDesignHeatmapSubSection(_GroupHeatmapSubSection):
    GROUP_KEY = DesignGroup.PAIR_DESIGN


class TripletHeatmapSubSection(_GroupHeatmapSubSection):
    GROUP_KEY = DesignGroup.TRIPLET


class QuadrupletHeatmapSubSection(_GroupHeatmapSubSection):
    GROUP_KEY = DesignGroup.QUADRUPLET
