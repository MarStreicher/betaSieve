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
from betasieve.report.plots import (
    ABS_Y_TITLE,
    add_abs_pct_hist_toggle,
    _layout_figure,
)
from betasieve.report.report_section import ReportMainSection, ReportSubSection

_X_MAX = 0.3
_X_TITLE = "β-values max-min difference"


def _hex_fill(color: str, alpha: float) -> str:
    r, g, b = int(color[1:3], 16), int(color[3:5], 16), int(color[5:7], 16)
    return f"rgba({r},{g},{b},{alpha})"


def _add_overlay_histogram(
    fig: go.Figure, values, name: str, color: str, shared_bins: dict
) -> None:
    finite = np.asarray(values).ravel()
    finite = finite[np.isfinite(finite)]
    fig.add_trace(
        go.Histogram(
            x=finite,
            name=name,
            marker=dict(color=_hex_fill(color, 0.4), line=dict(color=color, width=2)),
            opacity=0.7,
            xbins=shared_bins,
        )
    )


class DifferencesSection(ReportMainSection):
    @property
    def title(self) -> str:
        return "Max-Min Differences"

    @property
    def description(self) -> str:
        return "For each CpG-site-sample pair, the max-min difference was computed. "

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
            ExactVsOtherDifferencesHistogram,
            GroupedDifferencesHistogram,
            GroupedDifferencesBoxplot,
            *available_heatmaps,
        ]


class ExactVsOtherDifferencesHistogram(ReportSubSection):
    @property
    def title(self) -> str:
        return "β-values max-min difference: exact replicates vs design replicates"

    @property
    def description(self) -> str:
        return (
            "Distribution of β-value max-min differences, comparing exact "
            "replicates with all design replicate subgroups pooled together. "
            "Use the Absolute / Percentage "
            "buttons to switch the y-axis between counts and percentages."
        )

    def _plot(self) -> Figure:
        df = self.results.diff_frame
        fig = go.Figure()

        value_cols = [column for column in df.columns if column != Col.GROUP]
        shared_bins = dict(start=0, end=_X_MAX, size=_X_MAX / 100)

        exact_mask = df[Col.GROUP] == DesignGroup.EXACT_REPLICATES
        distributions = (
            ("Exact replicates", df.loc[exact_mask, value_cols], BS_GREEN),
            ("Design replicates", df.loc[~exact_mask, value_cols], "#F97316"),
        )

        for name, subset, color in distributions:
            _add_overlay_histogram(fig, subset.to_numpy(), name, color, shared_bins)

        fig.update_layout(barmode="overlay", bargap=0.05)
        fig.update_xaxes(range=[0, _X_MAX])
        fig = _layout_figure(
            fig,
            title="",
            x_title=_X_TITLE,
            y_title=ABS_Y_TITLE,
            height=400,
        )
        return add_abs_pct_hist_toggle(fig, n_traces=len(fig.data))

    def _figures(self):
        return [self._plot()]


class GroupedDifferencesHistogram(ReportSubSection):
    @property
    def title(self) -> str:
        return "β-values max-min difference: exact replicates vs design replicates subgroups"

    @property
    def description(self) -> str:
        return (
            "Distribution of the β-values max-min difference across all CpG-site-sample pairs, "
            "stratified by exact replicates and design replicate subgroups. "
            "Use the Absolute / Percentage buttons to switch the y-axis between "
            "counts and percentages."
        )

    def _plot(self) -> Figure:
        df = self.results.diff_frame
        fig = go.Figure()
        shared_bins = dict(start=0, end=_X_MAX, size=_X_MAX / 100)

        for group in DesignGroup:
            sub_df = df[df[Col.GROUP] == group]
            if sub_df.empty:
                continue
            values = sub_df.drop(columns=[Col.GROUP], errors="ignore").to_numpy()
            color = DESIGN_GROUP_COLORS.get(group, BS_GREEN)
            _add_overlay_histogram(fig, values, group.value, color, shared_bins)

        fig.update_layout(barmode="overlay", bargap=0.05)
        fig.update_xaxes(range=[0, _X_MAX])
        fig = _layout_figure(
            fig,
            title="",
            x_title=_X_TITLE,
            y_title=ABS_Y_TITLE,
            height=400,
        )
        return add_abs_pct_hist_toggle(fig, n_traces=len(fig.data))

    def _figures(self):
        return [self._plot()]


class GroupedDifferencesBoxplot(ReportSubSection):
    @property
    def title(self) -> str:
        return "β-values max-min difference: exact replicates vs design replicates subgroups (box plots)"

    @property
    def description(self) -> str:
        return (
            "Box plots of the β-values max-min difference pooled over all sites and "
            "samples, stratified by exact replicates and design replicate subgroups. Click a group in the legend to include or exclude it."
        )

    def _plot(self) -> Figure:
        diff = self.results.diff_frame
        value_cols = [col for col in diff.columns if col != Col.GROUP]

        fig = go.Figure()

        for group in DesignGroup:
            group_frame = diff[diff[Col.GROUP] == group]
            if group_frame.empty:
                continue
            values = group_frame[value_cols].to_numpy().ravel()
            color = DESIGN_GROUP_COLORS.get(group, BS_GREEN)
            fig.add_trace(
                go.Box(
                    y=values,
                    name=group.value,
                    boxpoints="outliers",
                    marker=dict(color=color, size=3, opacity=0.5),
                    line=dict(color=color),
                    fillcolor=_hex_fill(color, 0.15),
                )
            )

        return _layout_figure(
            fig,
            title="",
            x_title="Replicate group",
            y_title=_X_TITLE,
            height=400,
        )

    def _figures(self):
        return [self._plot()]


class _GroupHeatmapSubSection(ReportSubSection):
    GROUP_KEY: DesignGroup = DesignGroup.PAIR_TYPE

    @property
    def title(self) -> str:
        return f"Heatmap - {self.GROUP_KEY.value}"

    @property
    def description(self) -> str:
        return (
            f"Per CpG-site-sample pair β-values max-min differences for the "
            f"{self.GROUP_KEY.value} design replicate subgroups. "
            "Colour intensity reflects the magnitude of the max-min difference "
            "(square-root-transformed for visual clarity). "
            "Rows represent CpG-sites and the columns represent "
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
                text=f"β-values max-min difference heatmap - {group.value}",
                x=0.5,
                xanchor="center",
            ),
            xaxis_title="Sample",
            yaxis_title="CpG-site",
            height=max(420, min(1200, 14 * len(y_labels))),
            font=dict(family=REPORT_FONT_FAMILY),
        )
        fig.update_yaxes(autorange="reversed")
        return fig

    def _figures(self):
        df = self.results.diff_frame
        value_cols = [col for col in df.columns if col != Col.GROUP]
        if len(value_cols) == 0:
            return []

        non_exact_df = df[df[Col.GROUP] != DesignGroup.EXACT_REPLICATES]
        if non_exact_df.empty:
            return []

        global_vmax = float(np.nanmax(non_exact_df[value_cols].to_numpy()))
        fig = self._plot_for_group(self.GROUP_KEY, value_cols, global_vmax)
        return [fig] if fig is not None else []


class PairTypeHeatmapSubSection(_GroupHeatmapSubSection):
    GROUP_KEY = DesignGroup.PAIR_TYPE


class PairDesignHeatmapSubSection(_GroupHeatmapSubSection):
    GROUP_KEY = DesignGroup.PAIR_DESIGN


class TripletHeatmapSubSection(_GroupHeatmapSubSection):
    GROUP_KEY = DesignGroup.TRIPLET


class QuadrupletHeatmapSubSection(_GroupHeatmapSubSection):
    GROUP_KEY = DesignGroup.QUADRUPLET
