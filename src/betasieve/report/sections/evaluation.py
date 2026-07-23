from __future__ import annotations

from typing import List

import numpy as np
import plotly.graph_objects as go
from plotly.graph_objects import Figure
from scipy import stats

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
from betasieve.report.tables import _summary_table_figure


def _bb_thresholds(flagged_frame):
    er_mask = flagged_frame[Col.GROUP] == DesignGroup.EXACT_REPLICATES
    values = flagged_frame.loc[er_mask, Col.P_HAT].to_numpy()
    n = int(flagged_frame.loc[er_mask, Col.N].iloc[0])
    p0 = float(flagged_frame.loc[er_mask, Col.P0].iloc[0])
    confidence = float(flagged_frame[Col.CONFIDENCE].iloc[0])
    alpha = 1.0 - confidence

    k_binom = int(stats.binom.ppf(1 - alpha, n, p0)) + 1
    x_binom_crit = k_binom / n

    emp_var = float(np.var(values, ddof=1))
    rho = max(0.0, (emp_var * n / (p0 * (1 - p0)) - 1) / (n - 1))
    concentration = (1.0 / rho - 1.0) if rho > 0 else None

    x_bb_crit = 0.0
    if concentration is not None:
        k_bb = (
            int(
                stats.betabinom.ppf(
                    1 - alpha, n, p0 * concentration, (1 - p0) * concentration
                )
            )
            + 1
        )
        x_bb_crit = k_bb / n

    return x_binom_crit, x_bb_crit, concentration


class EvaluationSection(ReportMainSection):
    @property
    def title(self) -> str:
        return "Null Model Evaluation"

    @property
    def description(self) -> str:
        return (
            "Evaluation of the null distribution used for site-level statistical testing. "
            "The observed exceedance rate p̂ at exact-replicate sites provides an empirical estimate of the null distribution. "
            "Under random binomial sampling with rate p₀, these p̂ values are expected to follow a Binomial(n, p₀) distribution. "
            "Depending on the chosen threshold, the resulting distributions may appear skewed or exhibit overdispersion in the subsequent figures."
        )

    @property
    def subsection_types(self) -> List[type[ReportSubSection]]:
        return [
            H0Histogram,
            AllGroupsThresholdHistogram,
            EvaluationSummaryTableSubSection,
        ]


class H0Histogram(ReportSubSection):
    @property
    def title(self) -> str:
        return "Exact Replicates Distribution of p̂ Exact Replicates"

    @property
    def description(self) -> str:
        return (
            "Observed distribution of the exceedance rate p̂ at exact-replicate "
            "sites (bars), with the fitted Binomial(n, p₀) PMF (dashed black) and "
            "Beta-Binomial PMF (dotted green) overlaid. "
            "Both theoretical curves are scaled to the total number of "
            "exact-replicate sites for direct comparison with observed counts. "
            "The Beta-Binomial parameters are estimated by the method of moments "
            "from the empirical variance of p̂. "
            "Vertical lines mark the minimum p̂ at which a one-sided test at "
            "significance level α would declare a site significant. "
        )

    def _plot(self) -> Figure:
        df = self.results.flagged_frame
        fig = go.Figure()

        er_mask = df[Col.GROUP] == DesignGroup.EXACT_REPLICATES
        values = df.loc[er_mask, Col.P_HAT].to_numpy().ravel()
        color = DESIGN_GROUP_COLORS.get(DesignGroup.EXACT_REPLICATES, BS_GREEN)
        r, g, b = int(color[1:3], 16), int(color[3:5], 16), int(color[5:7], 16)
        fill = f"rgba({r},{g},{b},0.4)"

        n = int(df.loc[er_mask, Col.N].iloc[0])
        p0 = float(df.loc[er_mask, Col.P0].iloc[0])
        N = len(values)
        k = np.arange(0, n + 1)
        x_pmf = k / n

        bin_size = 1.0 / n
        fig.add_trace(
            go.Histogram(
                x=values,
                name=DesignGroup.EXACT_REPLICATES.value,
                marker=dict(color=fill, line=dict(color=color, width=2)),
                opacity=0.7,
                xbins=dict(start=-bin_size / 2, end=1 + bin_size / 2, size=bin_size),
            )
        )

        binom_y = stats.binom.pmf(k, n, p0) * N
        fig.add_trace(
            go.Scatter(
                x=x_pmf,
                y=binom_y,
                mode="lines",
                name=f"Binomial null (p₀={p0:.3f})",
                line=dict(color=BS_DARK, width=2, dash="dash"),
            )
        )

        emp_var = float(np.var(values, ddof=1))
        rho = max(0.0, (emp_var * n / (p0 * (1 - p0)) - 1) / (n - 1))
        concentration = (1.0 / rho - 1.0) if rho > 0 else None
        if concentration is not None:
            bb_y = (
                stats.betabinom.pmf(k, n, p0 * concentration, (1 - p0) * concentration)
                * N
            )
            fig.add_trace(
                go.Scatter(
                    x=x_pmf,
                    y=bb_y,
                    mode="lines",
                    name=f"Beta-Binomial fit (ρ={rho:.3f})",
                    line=dict(color=BS_GREEN, width=2, dash="dot"),
                )
            )

        y_max = float(max(binom_y[1:]))  # skip k=0

        x_binom_crit, x_bb_crit, _ = _bb_thresholds(self.results.flagged_frame)
        alpha = 1.0 - float(self.results.flagged_frame[Col.CONFIDENCE].iloc[0])

        nonzero = values[values > 0]
        x_max = (
            float(
                max(
                    np.quantile(nonzero, 0.99) if len(nonzero) else alpha,
                    x_binom_crit,
                    x_bb_crit,
                )
            )
            + bin_size
        )

        fig.add_trace(
            go.Scatter(
                x=[x_binom_crit, x_binom_crit],
                y=[0, y_max],
                mode="lines",
                name=f"Binomial threshold (p̂={x_binom_crit:.3f})",
                line=dict(color="red", width=2, dash="dash"),
            )
        )

        if concentration is not None:
            # if x_bb_crit > x_binom_crit:
            #    fig.add_vrect(
            #        x0=x_binom_crit,
            #        x1=x_bb_crit,
            #        fillcolor="rgba(220, 50, 50, 0.12)",
            #        layer="below",
            #        line_width=0,
            #        annotation_text="Binomial only",
            #        annotation_position="top left",
            #        annotation_font_size=11,
            #    )

            fig.add_trace(
                go.Scatter(
                    x=[x_bb_crit, x_bb_crit],
                    y=[0, y_max],
                    mode="lines",
                    name=f"Beta-Binomial threshold (p̂={x_bb_crit:.3f})",
                    line=dict(color=BS_GREEN, width=2, dash="dash"),
                )
            )

        fig_out = _layout_figure(
            fig,
            title="Observed vs theoretical p̂ distribution — exact-replicate sites",
            x_title="p̂ (observed exceedance rate)",
            y_title="Number of exact-replicate sites",
            height=400,
        )
        fig_out.update_layout(
            xaxis=dict(range=[0, x_max]),
            yaxis=dict(range=[0, y_max * 1.15]),
        )
        return fig_out

    def generate(self) -> None:
        self.figures.append(self._plot())


class AllGroupsThresholdHistogram(ReportSubSection):
    @property
    def title(self) -> str:
        return "p̂ Distribution Across All Design Groups vs Beta-Binomial Threshold"

    @property
    def description(self) -> str:
        return (
            "Distribution of the observed exceedance rate p̂ for non-replicate "
            "design groups (pairs, triplets, quadruplets), with the Beta-Binomial "
            "significance threshold (green dashed) and Binomial "
            "threshold (red dashed) overlaid. "
            "Sites to the right of the Beta-Binomial threshold are declared "
            "significant. "
        )

    def _plot(self) -> Figure:
        flagged = self.results.flagged_frame
        x_binom_crit, x_bb_crit, _ = _bb_thresholds(flagged)

        fig = go.Figure()
        all_values = []

        for group in DesignGroup:
            if group == DesignGroup.EXACT_REPLICATES:
                continue  # shown separately in H0Histogram
            mask = flagged[Col.GROUP] == group
            values = flagged.loc[mask, Col.P_HAT].dropna().to_numpy()
            if len(values) == 0:
                continue
            all_values.append(values)

            color = DESIGN_GROUP_COLORS.get(group, BS_GREEN)
            r, g, b = int(color[1:3], 16), int(color[3:5], 16), int(color[5:7], 16)
            fill = f"rgba({r},{g},{b},0.4)"
            fig.add_trace(
                go.Histogram(
                    x=values,
                    name=group.value,
                    marker=dict(color=fill, line=dict(color=color, width=2)),
                    opacity=0.7,
                    nbinsx=60,
                )
            )

        fig.update_layout(barmode="overlay")

        all_concat = np.concatenate(all_values) if all_values else np.array([0.0])
        x_max = float(max(np.quantile(all_concat, 0.99), x_bb_crit)) + 0.02

        fig.add_vline(
            x=x_binom_crit,
            line_color="red",
            line_width=2,
            line_dash="dash",
        )
        fig.add_trace(
            go.Scatter(
                x=[None],
                y=[None],
                mode="lines",
                name=f"Binomial threshold (p̂={x_binom_crit:.3f})",
                line=dict(color="red", width=2, dash="dash"),
                showlegend=True,
            )
        )

        if x_bb_crit > 0:
            # if x_bb_crit > x_binom_crit:
            #    fig.add_vrect(
            #        x0=x_binom_crit,
            #        x1=x_bb_crit,
            #        fillcolor="rgba(220, 50, 50, 0.10)",
            #        layer="below",
            #        line_width=0,
            #        annotation_text="Binomial only",
            #        annotation_position="top left",
            #        annotation_font_size=11,
            #    )
            fig.add_vline(
                x=x_bb_crit,
                line_color=BS_GREEN,
                line_width=2,
                line_dash="dash",
            )
            fig.add_trace(
                go.Scatter(
                    x=[None],
                    y=[None],
                    mode="lines",
                    name=f"Beta-Binomial threshold (p̂={x_bb_crit:.3f})",
                    line=dict(color=BS_GREEN, width=2, dash="dash"),
                    showlegend=True,
                )
            )

        fig_out = _layout_figure(
            fig,
            title="p̂ distribution by design group with significance thresholds",
            x_title="p̂ (observed exceedance rate)",
            y_title="Number of CpG sites",
            height=420,
        )
        fig_out.update_layout(xaxis=dict(range=[0, x_max]))
        return fig_out

    def generate(self) -> None:
        self.figures.append(self._plot())


class EvaluationSummaryTableSubSection(ReportSubSection):
    @property
    def title(self) -> str:
        return "Flagging Summary (Classical Tests)"

    @property
    def description(self) -> str:
        return (
            "Summary of flagged sites under the classical z-test and "
            "Wilson confidence interval approaches.\n"
            "<strong>Note:</strong> this table is provided for reference only. "
            "The primary flagging criterion used by betaSieve is the "
            "FDR-adjusted permutation p-value, which does not depend "
            "on this approximation of a normal distribution and is reported in the Analysis Configuration "
            "section."
        )

    def _summary_table(self) -> Figure:
        flagged = self.results.flagged_frame
        rows = [
            ("Threshold", round(float(flagged[Col.THRESHOLD].iloc[0]), 4)),
            ("p₀", round(float(flagged[Col.P0].iloc[0]), 4)),
            ("Sites with CI flag", f"{int(flagged[Col.CI_FLAG].sum()):,}"),
            ("Sites with p-flag", f"{int(flagged[Col.P_FLAG].sum()):,}"),
            ("Sites with adjusted p-flag", f"{int(flagged[Col.P_ADJ_FLAG].sum()):,}"),
        ]
        return _summary_table_figure(rows)

    def generate(self) -> None:
        self.figures.append(self._summary_table())
