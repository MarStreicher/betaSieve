from typing import Any, Optional, Sequence, Tuple

import plotly.graph_objects as go
from plotly.graph_objects import Figure

from betasieve.report.domain.mappings import (
    BS_DARK,
    BS_LIGHTEST,
    REPORT_FONT_FAMILY,
    TABLE_HEADER_BG,
)

_ROW_WHITE = "#FFFFFF"
_CELL_LINE = "rgba(0,0,0,0)"


def _summary_table_figure(rows: Sequence[Tuple[str, Any]]) -> Figure:
    """Parameter/value table styled with betaSieve report colors (no grid lines)."""
    labels = [r[0] for r in rows]
    values = [str(r[1]) for r in rows]
    n = len(rows)
    row_colors = [_ROW_WHITE if i % 2 == 0 else BS_LIGHTEST for i in range(n)]

    fig = go.Figure(
        data=[
            go.Table(
                columnwidth=[0.42, 0.58],
                header=dict(
                    values=["Parameter", "Value"],
                    fill_color=TABLE_HEADER_BG,
                    font=dict(
                        color=_ROW_WHITE,
                        size=14,
                        family=REPORT_FONT_FAMILY,
                    ),
                    align=["left", "left"],
                    line=dict(color=TABLE_HEADER_BG, width=0),
                    height=32,
                ),
                cells=dict(
                    values=[labels, values],
                    fill_color=[row_colors, row_colors],
                    font=dict(
                        color=BS_DARK,
                        size=14,
                        family=REPORT_FONT_FAMILY,
                    ),
                    align=["left", "left"],
                    line=dict(color=_CELL_LINE, width=0),
                    height=28,
                ),
            )
        ]
    )
    fig.update_layout(
        margin=dict(l=12, r=12, t=8, b=8),
        height=48 + 28 * n,
        paper_bgcolor=_ROW_WHITE,
        plot_bgcolor=_ROW_WHITE,
    )
    return fig


def _data_dict_figure(
    headers: Sequence[str],
    rows: Sequence[Tuple],
    col_widths: Optional[Sequence[float]] = None,
    row_height: int = 30,
) -> Figure:
    """
    Multi-column data-dictionary table.

    Parameters
    ----------
    headers : column header labels
    rows    : sequence of tuples, one per data row — must match len(headers)
    col_widths : relative column widths (defaults to equal distribution)
    row_height : pixel height per data row
    """
    n = len(rows)
    n_cols = len(headers)
    row_colors = [_ROW_WHITE if i % 2 == 0 else BS_LIGHTEST for i in range(n)]
    col_values = [[str(row[i]) for row in rows] for i in range(n_cols)]

    fig = go.Figure(
        data=[
            go.Table(
                columnwidth=col_widths or [1] * n_cols,
                header=dict(
                    values=list(headers),
                    fill_color=TABLE_HEADER_BG,
                    font=dict(
                        color=_ROW_WHITE,
                        size=13,
                        family=REPORT_FONT_FAMILY,
                    ),
                    align=["left"] * n_cols,
                    line=dict(color=TABLE_HEADER_BG, width=0),
                    height=32,
                ),
                cells=dict(
                    values=col_values,
                    fill_color=[row_colors] * n_cols,
                    font=dict(
                        color=BS_DARK,
                        size=12,
                        family=REPORT_FONT_FAMILY,
                    ),
                    align=["left"] * n_cols,
                    line=dict(color=_CELL_LINE, width=0),
                    height=row_height,
                ),
            )
        ]
    )
    fig.update_layout(
        margin=dict(l=12, r=12, t=8, b=8),
        height=48 + row_height * n,
        paper_bgcolor=_ROW_WHITE,
        plot_bgcolor=_ROW_WHITE,
    )
    return fig
