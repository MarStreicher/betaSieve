from typing import Any, Optional, Sequence, Tuple

import plotly.graph_objects as go
from plotly.graph_objects import Figure

from betasieve.report.domain.mappings import (
    TABLE_CELL_TEXT,
    TABLE_FONT_FAMILY,
    TABLE_HEADER_BG,
    TABLE_HEADER_TEXT,
    TABLE_ROW_ALT,
    TABLE_ROW_BG,
    TABLE_RULE,
)


def _summary_table_figure(rows: Sequence[Tuple[str, Any]]) -> Figure:
    """Render a compact, publication-style parameter/value table."""
    labels = [r[0] for r in rows]
    values = [str(r[1]) for r in rows]
    n = len(rows)
    row_colors = [TABLE_ROW_BG if i % 2 == 0 else TABLE_ROW_ALT for i in range(n)]

    fig = go.Figure(
        data=[
            go.Table(
                columnwidth=[0.42, 0.58],
                header=dict(
                    values=["Parameter", "Value"],
                    fill_color=TABLE_HEADER_BG,
                    font=dict(
                        color=TABLE_HEADER_TEXT,
                        size=13,
                        family=TABLE_FONT_FAMILY,
                    ),
                    align=["left", "left"],
                    line=dict(color=TABLE_RULE, width=1),
                    height=34,
                ),
                cells=dict(
                    values=[labels, values],
                    fill_color=[row_colors, row_colors],
                    font=dict(
                        color=TABLE_CELL_TEXT,
                        size=13,
                        family=TABLE_FONT_FAMILY,
                    ),
                    align=["left", "left"],
                    line=dict(color=TABLE_RULE, width=0.5),
                    height=30,
                ),
            )
        ]
    )
    fig.update_layout(
        margin=dict(l=12, r=12, t=8, b=8),
        height=50 + 30 * n,
        paper_bgcolor=TABLE_ROW_BG,
        plot_bgcolor=TABLE_ROW_BG,
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
    row_colors = [TABLE_ROW_BG if i % 2 == 0 else TABLE_ROW_ALT for i in range(n)]
    col_values = [[str(row[i]) for row in rows] for i in range(n_cols)]

    fig = go.Figure(
        data=[
            go.Table(
                columnwidth=col_widths or [1] * n_cols,
                header=dict(
                    values=list(headers),
                    fill_color=TABLE_HEADER_BG,
                    font=dict(
                        color=TABLE_HEADER_TEXT,
                        size=13,
                        family=TABLE_FONT_FAMILY,
                    ),
                    align=["left"] * n_cols,
                    line=dict(color=TABLE_RULE, width=1),
                    height=34,
                ),
                cells=dict(
                    values=col_values,
                    fill_color=[row_colors] * n_cols,
                    font=dict(
                        color=TABLE_CELL_TEXT,
                        size=12,
                        family=TABLE_FONT_FAMILY,
                    ),
                    align=["left"] * n_cols,
                    line=dict(color=TABLE_RULE, width=0.5),
                    height=row_height,
                ),
            )
        ]
    )
    fig.update_layout(
        margin=dict(l=12, r=12, t=8, b=8),
        height=50 + row_height * n,
        paper_bgcolor=TABLE_ROW_BG,
        plot_bgcolor=TABLE_ROW_BG,
    )
    return fig
