from __future__ import annotations

from betasieve.cg_probe_table import DesignGroup

# Neutral base with logo-green report accents
BS_DARK = "#2F3437"
BS_MID = "#6B7280"
BS_GREEN = "#4A934A"
BS_LIGHT = "#9CA3AF"
BS_LIGHTEST = "#F4F5F6"

BS_HEATMAP_COLORSCALE: list[list[float | str]] = [
    [0.0, BS_LIGHTEST],
    [0.25, "#E2F0E2"],
    [0.5, "#8FBC8F"],
    [0.75, BS_GREEN],
    [1.0, "#2563EB"],
]


REPORT_FONT_FAMILY = "Arial, Helvetica, sans-serif"
TABLE_FONT_FAMILY = "Montserrat, Arial, sans-serif"
TABLE_HEADER_BG = "#E8F3E8"
TABLE_HEADER_TEXT = "#285C2D"
TABLE_ROW_BG = "#FFFFFF"
TABLE_ROW_ALT = "#F7FAF7"
TABLE_CELL_TEXT = BS_DARK
TABLE_RULE = "#C9DDC9"

DESIGN_GROUP_COLORS: dict[DesignGroup, str] = {
    DesignGroup.PAIR_TYPE: "#2563EB",
    DesignGroup.PAIR_DESIGN: "#F97316",
    DesignGroup.TRIPLET: "#DC2626",
    DesignGroup.QUADRUPLET: "#6B7280",
    DesignGroup.EXACT_REPLICATES: BS_GREEN,
}

BS_CHART_PALETTE = [
    BS_GREEN,
    "#5C6370",
    "#9CA3AF",
    BS_DARK,
    "#B8BFC6",
    "#7A8490",
    "#3D4449",
    "#6FA86F",
]
