from __future__ import annotations

from typing import List, NamedTuple

import numpy as np
import plotly.graph_objects as go
from plotly.graph_objects import Figure

from betasieve.analysis import SieveResults
from betasieve.columns import Col
from betasieve.cg_probe_table import DesignGroup
from betasieve.null_models import NullModels
from betasieve.report.domain.mappings import (
    BS_DARK,
    BS_GREEN,
    DESIGN_GROUP_COLORS,
)
from betasieve.report.plots import _layout_figure
from betasieve.report.report_section import ReportMainSection, ReportSubSection
from betasieve.report.tables import _summary_table_figure


def _null_models(results: SieveResults) -> NullModels:
    return results.null_models


class _Cutoff(NamedTuple):
    label: str
    p_hat: float
    color: str
    dash: str


def _null_cutoffs(models: NullModels) -> List[_Cutoff]:
    alpha = models.alpha
    specs = (
        ("Binomial", models.binomial.critical_p_hat(alpha), "red", "dash"),
        ("Beta-Binomial", models.beta_binomial.critical_p_hat(alpha), BS_GREEN, "dash"),
        ("Empirical", models.empirical.critical_p_hat(alpha), "#2563EB", "dashdot"),
    )
    return [
        _Cutoff(f"{name} threshold (<i>p̂</i>={p_hat:.3f})", p_hat, color, dash)
        for name, p_hat, color, dash in specs
        if p_hat is not None
    ]


def _add_cutoff_vlines(fig: Figure, cutoffs: List[_Cutoff]) -> None:
    """Draw cutoffs as traces so legend clicks hide the lines."""
    fig.update_layout(
        yaxis2=dict(
            overlaying="y",
            side="right",
            range=[0, 1],
            showticklabels=False,
            showgrid=False,
            zeroline=False,
            showline=False,
            ticks="",
        )
    )
    for cutoff in cutoffs:
        fig.add_trace(
            go.Scatter(
                x=[cutoff.p_hat, cutoff.p_hat],
                y=[0, 1],
                mode="lines",
                yaxis="y2",
                name=cutoff.label,
                line=dict(color=cutoff.color, width=2, dash=cutoff.dash),
                hoverinfo="skip",
                showlegend=True,
            )
        )


def _group_fill(color: str) -> str:
    r, g, b = int(color[1:3], 16), int(color[3:5], 16), int(color[5:7], 16)
    return f"rgba({r},{g},{b},0.4)"


class EvaluationSection(ReportMainSection):
    @property
    def title(self) -> str:
        return "Null Model Evaluation"

    @property
    def description(self) -> str:
        return (
            "Comparison of three null representations for CpG-site-level exceedance rates: "
            "the observed exact replicate distribution, the Binomial(<em>n</em>, <em>p₀</em>) model, "
            "and a method-of-moments Beta-Binomial model. The adjusted empirical "
            "upper-tail test is the primary flagging criterion and the model-based "
            "distributions are shown as diagnostics."
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
        return "Empirical and Modelled exact replicate Null Distributions"

    @property
    def description(self) -> str:
        return (
            "Observed distribution of the exceedance rate <em>p̂</em> at exact replicate "
            "CpG-sites (bars), with the fitted Binomial(<em>n</em>, <em>p₀</em>) PMF (dashed black) and "
            "Beta-Binomial PMF (dotted green) overlaid. "
            "All three are expressed as a percentage of exact replicate CpG-sites, so "
            "bar heights and PMF values are directly comparable. "
            "The Beta-Binomial parameters are estimated by the method of moments "
            "from the empirical variance of <em>p̂</em>. "
            "Vertical lines mark the minimum <em>p̂</em> meeting each unadjusted one-sided "
            "significance criterion at level <em>α</em>. Final candidate selection additionally "
            "applies the configured multiple-testing correction."
        )

    def _plot(self) -> Figure:
        models = _null_models(self.results)
        df = self.results.flagged_frame
        fig = go.Figure()

        er_mask = df[Col.GROUP] == DesignGroup.EXACT_REPLICATES
        values = df.loc[er_mask, Col.P_HAT].to_numpy().ravel()
        color = DESIGN_GROUP_COLORS.get(DesignGroup.EXACT_REPLICATES, BS_GREEN)

        n = models.n
        k = np.arange(0, n + 1)
        x_pmf = k / n
        bin_size = 1.0 / n

        fig.add_trace(
            go.Histogram(
                x=values,
                histnorm="percent",
                name=DesignGroup.EXACT_REPLICATES.value,
                marker=dict(color=_group_fill(color), line=dict(color=color, width=2)),
                opacity=0.7,
                xbins=dict(start=-bin_size / 2, end=1 + bin_size / 2, size=bin_size),
            )
        )

        binom_y = models.binomial.pmf(k) * 100.0
        fig.add_trace(
            go.Scatter(
                x=x_pmf,
                y=binom_y,
                mode="lines",
                name=f"Binomial null (<i>p₀</i>={models.binomial.p0:.3f})",
                line=dict(color=BS_DARK, width=2, dash="dash"),
            )
        )

        bb_y = models.beta_binomial.pmf(k) * 100.0
        fig.add_trace(
            go.Scatter(
                x=x_pmf,
                y=bb_y,
                mode="lines",
                name=(
                    f"Beta-Binomial fit "
                    f"(<i>α</i>={models.beta_binomial.a:.3g}, <i>β</i>={models.beta_binomial.b:.3g})"
                ),
                line=dict(color=BS_GREEN, width=2, dash="dot"),
            )
        )

        observed_pct = (
            100.0
            * np.bincount(np.rint(values * n).astype(int), minlength=n + 1)
            / max(len(values), 1)
        )
        y_max = float(max(binom_y[1:].max(), bb_y[1:].max(), observed_pct[1:].max()))

        cutoffs = _null_cutoffs(models)
        for cutoff in cutoffs:
            fig.add_trace(
                go.Scatter(
                    x=[cutoff.p_hat, cutoff.p_hat],
                    y=[0, y_max],
                    mode="lines",
                    name=cutoff.label,
                    line=dict(color=cutoff.color, width=2, dash=cutoff.dash),
                )
            )

        nonzero = values[values > 0]
        x_max = (
            float(
                max(
                    [np.quantile(nonzero, 0.99) if len(nonzero) else models.alpha]
                    + [cutoff.p_hat for cutoff in cutoffs]
                )
            )
            + bin_size
        )

        fig_out = _layout_figure(
            fig,
            title="Empirical, Binomial, and Beta-Binomial null comparison",
            x_title="<i>p̂</i> (observed exceedance rate)",
            y_title="% of exact replicate CpG-sites",
            height=400,
        )
        fig_out.update_layout(
            xaxis=dict(range=[0, x_max]),
            yaxis=dict(range=[0, y_max * 1.15]),
        )
        return fig_out

    def _figures(self):
        return [self._plot()]


class AllGroupsThresholdHistogram(ReportSubSection):
    @property
    def title(self) -> str:
        return "<em>p̂</em> distribution across design replicates subgroups and null cutoffs"

    @property
    def description(self) -> str:
        return (
            "Distribution of the observed exceedance rate <em>p̂</em> for design replicate "
            "subgroups, with the "
            "unadjusted empirical, Binomial, and Beta-Binomial significance cutoffs "
            "overlaid. Each group is normalised to its own site count, since the groups "
            "differ by orders of magnitude in size. "
            "Click a cutoff in the legend to hide or show its vertical line."
        )

    def _plot(self) -> Figure:
        models = _null_models(self.results)
        flagged = self.results.flagged_frame
        cutoffs = _null_cutoffs(models)
        bin_size = 1.0 / models.n

        fig = go.Figure()
        all_values = []

        for group in DesignGroup:
            if group == DesignGroup.EXACT_REPLICATES:
                continue
            mask = flagged[Col.GROUP] == group
            values = flagged.loc[mask, Col.P_HAT].dropna().to_numpy()
            if len(values) == 0:
                continue
            all_values.append(values)

            color = DESIGN_GROUP_COLORS.get(group, BS_GREEN)
            fig.add_trace(
                go.Histogram(
                    x=values,
                    name=group.value,
                    histnorm="percent",
                    marker=dict(
                        color=_group_fill(color), line=dict(color=color, width=2)
                    ),
                    opacity=0.7,
                    # p̂ is discrete on multiples of 1/n; auto-binning would smear
                    # those atoms unevenly across bins.
                    xbins=dict(
                        start=-bin_size / 2, end=1 + bin_size / 2, size=bin_size
                    ),
                )
            )

        fig.update_layout(barmode="overlay")

        all_concat = np.concatenate(all_values) if all_values else np.array([0.0])
        x_max = (
            float(
                max(
                    [np.quantile(all_concat, 0.99)]
                    + [cutoff.p_hat for cutoff in cutoffs]
                )
            )
            + 0.02
        )

        _add_cutoff_vlines(fig, cutoffs)

        fig_out = _layout_figure(
            fig,
            title="<i>p̂</i> distribution by design replicate group with three null cutoffs",
            x_title="<i>p̂</i> (observed exceedance rate)",
            y_title="% of CpG-sites in replicate group",
            height=420,
        )
        fig_out.update_layout(xaxis=dict(range=[0, x_max]))
        return fig_out

    def _figures(self):
        return [self._plot()]


class EvaluationSummaryTableSubSection(ReportSubSection):
    @property
    def title(self) -> str:
        return "Flagging Summary by Statistical Criterion"

    @property
    def description(self) -> str:
        return (
            "Fitted null parameters and site counts under the empirical upper-tail, "
            "Beta-Binomial, classical z-test, and Wilson lower-bound criteria. "
            "The adjusted empirical upper-tail criterion is used to create the "
            "candidate list."
        )

    def _summary_table(self) -> Figure:
        models = _null_models(self.results)
        flagged = self.results.flagged_frame
        beta_binomial = models.beta_binomial
        rows = [
            ("Threshold", round(models.threshold, 4)),
            ("𝑝₀", round(models.binomial.p0, 4)),
            ("exact replicate CpG-sites (𝑚)", f"{models.empirical.m:,}"),
            ("Samples (𝑛)", f"{models.n:,}"),
            (
                "Beta-Binomial fit",
                f"𝛼={beta_binomial.a:.4g}, 𝛽={beta_binomial.b:.4g}",
            ),
            (
                "CpG-sites with empirical flag",
                f"{int(flagged[Col.P_EMPIR_FLAG].sum()):,}",
            ),
            (
                "CpG-sites with adjusted empirical flag",
                f"{int(flagged[Col.P_EMPIR_ADJ_FLAG].sum()):,}",
            ),
            (
                "CpG-sites with adjusted Beta-Binomial flag",
                f"{int(flagged[Col.P_BETA_ADJ_FLAG].sum()):,}",
            ),
            ("CpG-sites with CI flag", f"{int(flagged[Col.CI_FLAG].sum()):,}"),
            ("CpG-sites with p-flag", f"{int(flagged[Col.P_FLAG].sum()):,}"),
            (
                "CpG-sites with adjusted p-flag",
                f"{int(flagged[Col.P_ADJ_FLAG].sum()):,}",
            ),
        ]
        return _summary_table_figure(rows)

    def _figures(self):
        return [self._summary_table()]
