from __future__ import annotations

from betasieve.cg_probe_table import DesignGroup

# Neutral base with logo-green report accents
BS_DARK = "#2F3437"      # charcoal — headings, primary text
BS_MID = "#6B7280"       # slate gray — secondary text, muted UI
BS_GREEN = "#4A934A"     # brand green — links, borders, key highlights
BS_LIGHT = "#9CA3AF"     # cool gray — tertiary series / subtle accents
BS_LIGHTEST = "#F4F5F6"  # off-white gray — panel / TOC background

# Sequential heatmap: panel background → soft green → brand green → report blue
BS_HEATMAP_COLORSCALE: list[list[float | str]] = [
    [0.0, BS_LIGHTEST],
    [0.25, "#E2F0E2"],
    [0.5, "#8FBC8F"],
    [0.75, BS_GREEN],
    [1.0, "#2563EB"],
]

# Match template.css body typography
REPORT_FONT_FAMILY = "Arial, Helvetica, sans-serif"
TABLE_HEADER_BG = "#7A8490"  # softer than BS_DARK for Plotly table headers

# One color per design group, keyed by DesignGroup enum member
DESIGN_GROUP_COLORS: dict[DesignGroup, str] = {
    DesignGroup.PAIR_TYPE: "#2563EB",
    DesignGroup.PAIR_DESIGN: "#F97316",
    DesignGroup.TRIPLET: "#DC2626",
    DesignGroup.QUADRUPLET: "#6B7280",
    DesignGroup.EXACT_REPLICATES: BS_GREEN,
}

# Ordered palette for charts that need many distinct colors
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
