"""Plotly figure factory for betaSieve HTML reports."""

from __future__ import annotations

import pandas as pd
import plotly.graph_objects as go
from plotly.graph_objects import Figure

from betasieve.analysis import Col
from betasieve.report.domain.mappings import (
    BS_DARK,
    BS_GREEN,
    BS_MID,
    DESIGN_GROUP_COLORS,
    REPORT_FONT_FAMILY,
)


def _layout_figure(
    fig: Figure,
    *,
    title: str,
    x_title: str,
    y_title: str,
    height: int,
    show_legend: bool = True,
) -> Figure:
    fig.update_layout(
        template="plotly_white",
        title=dict(text=title, x=0.5, xanchor="center"),
        xaxis_title=x_title,
        yaxis_title=y_title,
        yaxis_rangemode="tozero",
        height=height,
        font=dict(family=REPORT_FONT_FAMILY),
        showlegend=show_legend,
        legend_title_text="Design group" if show_legend else None,
    )
    return fig
