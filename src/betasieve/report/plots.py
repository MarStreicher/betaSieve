from __future__ import annotations

from plotly.graph_objects import Figure

from betasieve.report.domain.mappings import BS_GREEN, REPORT_FONT_FAMILY

ABS_Y_TITLE = "Number of observations"
PCT_Y_TITLE = "Percentage of observations"


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
        legend_title_text="Replicate group" if show_legend else None,
    )
    return fig


def add_abs_pct_hist_toggle(
    fig: Figure,
    *,
    n_traces: int,
    accent_color: str = BS_GREEN,
    abs_y_title: str = ABS_Y_TITLE,
    pct_y_title: str = PCT_Y_TITLE,
) -> Figure:
    """Add Plotly buttons that switch histogram traces between counts and percent."""
    fig.update_layout(
        margin=dict(t=100),
        updatemenus=[
            dict(
                type="buttons",
                direction="left",
                showactive=True,
                x=1.0,
                xanchor="right",
                y=1.0,
                yanchor="bottom",
                pad=dict(t=0, r=0, b=8, l=0),
                bgcolor="white",
                bordercolor=accent_color,
                borderwidth=1,
                font=dict(size=12, color=accent_color),
                buttons=[
                    dict(
                        label="Absolute",
                        method="update",
                        args=[
                            {"histnorm": [""] * n_traces},
                            {
                                "yaxis.title.text": abs_y_title,
                                "yaxis.ticksuffix": "",
                            },
                        ],
                    ),
                    dict(
                        label="Percentage",
                        method="update",
                        args=[
                            {"histnorm": ["percent"] * n_traces},
                            {
                                "yaxis.title.text": pct_y_title,
                                "yaxis.ticksuffix": "%",
                            },
                        ],
                    ),
                ],
            )
        ],
    )
    return fig
