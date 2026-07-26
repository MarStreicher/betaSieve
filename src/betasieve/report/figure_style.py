"""Shared Matplotlib rcParams for publication figures (LaTeX-like serif + CM math)."""

from __future__ import annotations

from typing import Any

import matplotlib.pyplot as plt

DEFAULT_DPI = 300
DEFAULT_FIG_WIDTH = 14.0
DEFAULT_ROW_HEIGHT = 8.0
DEFAULT_FONT_SIZE = 14
PANEL_LETTER_FONT_SIZE = 16

_LATEX_LIKE_RCPARAMS: dict[str, Any] = {
    "font.family": "serif",
    "font.serif": ["CMU Serif", "STIXGeneral", "DejaVu Serif"],
    "mathtext.fontset": "cm",
    "axes.unicode_minus": False,
    "font.size": DEFAULT_FONT_SIZE,
    "axes.titlesize": DEFAULT_FONT_SIZE,
    "axes.labelsize": DEFAULT_FONT_SIZE,
    "legend.fontsize": DEFAULT_FONT_SIZE,
    "xtick.labelsize": DEFAULT_FONT_SIZE,
    "ytick.labelsize": DEFAULT_FONT_SIZE,
}


def configure_matplotlib() -> None:
    """Apply shared typography and sizes; call before creating figures."""
    plt.rcParams.update(_LATEX_LIKE_RCPARAMS)
