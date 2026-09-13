from typing import List

from plotly.graph_objects import Figure

from betasieve.columns import Col
from betasieve.cg_probe_table import DesignGroup
from betasieve.report.report_section import ReportMainSection, ReportSubSection
from betasieve.report.tables import _summary_table_figure


class ConfigSection(ReportMainSection):
    @property
    def title(self) -> str:
        return "Configuration"

    @property
    def description(self) -> str:
        return (
            "This section summarises the parameters used for the analysis and key "
            "descriptive statistics derived from the input data."
        )

    @property
    def subsection_types(self) -> List[type[ReportSubSection]]:
        return [ConfigSubSection, ConfigResultsSubSection]


class ConfigSubSection(ReportSubSection):
    @property
    def title(self) -> str:
        return "Parameters"

    @property
    def description(self) -> str:
        return (
            "The table summarizes the analysis parameters used by <strong>betaSieve</strong>. "
            "Please note that you can specify the threshold <em>t</em> directly or determine it automatically through a threshold sweep. "
        )

    def _summary_table(self) -> Figure:
        args = self.args
        config = args.analysis

        rows = [
            ("Betas file", str(args.betas_path)),
            (
                ("Threshold 𝑡", str(config.threshold))
                if config.threshold is not None
                else ("Threshold 𝑡", "not used (threshold sweep)")
            ),
            ("Confidence", str(config.confidence)),
            ("FDR method", config.fdr),
            (
                "Target 𝑝₀ (threshold sweep)",
                (
                    str(config.target_p0)
                    if config.threshold is None
                    else "not used (fixed threshold)"
                ),
            ),
            (
                "Threshold sweep min / max / step",
                (
                    f"{config.threshold_min} / {config.threshold_max} / "
                    f"{config.threshold_step}"
                    if config.threshold_min is not None
                    else "not used (fixed threshold)"
                ),
            ),
        ]
        return _summary_table_figure(rows)

    def _figures(self):
        return [self._summary_table()]


class ConfigResultsSubSection(ReportSubSection):
    @property
    def title(self) -> str:
        return "Dataset Summary"

    @property
    def description(self) -> str:
        return (
            "Key statistics derived from the input data and analysis results. "
            "Please note that individual CpG-sites may appear in more than one "
            "replicate group (e.g. a CpG-site belonging to both pair type "
            "and to exact replicates). "
            "The selected threshold <em>t</em> and corresponding <em>p₀</em> are the values used for "
            "all downstream statistical tests and figures."
        )

    def _summary_table(self) -> Figure:
        flagged = self.results.flagged_frame
        group = flagged[Col.GROUP]

        rows = [
            ("Number of samples (𝑛)", int(flagged[Col.N].iloc[0])),
            ("Total CpG-sites analysed", f"{len(flagged):,}"),
            (
                f"{DesignGroup.PAIR_TYPE.value} sites",
                int((group == DesignGroup.PAIR_TYPE).sum()),
            ),
            (
                f"{DesignGroup.PAIR_DESIGN.value} sites",
                int((group == DesignGroup.PAIR_DESIGN).sum()),
            ),
            (
                f"{DesignGroup.TRIPLET.value} sites",
                int((group == DesignGroup.TRIPLET).sum()),
            ),
            (
                f"{DesignGroup.QUADRUPLET.value} sites",
                int((group == DesignGroup.QUADRUPLET).sum()),
            ),
            (
                f"{DesignGroup.EXACT_REPLICATES.value} sites",
                int((group == DesignGroup.EXACT_REPLICATES).sum()),
            ),
            (
                "Selected threshold 𝑡",
                round(float(flagged[Col.THRESHOLD].iloc[0]), 4),
            ),
            (
                "Background exceedance rate 𝑝₀",
                round(float(flagged[Col.P0].iloc[0]), 4),
            ),
            (
                "Sites flagged",
                f"{int(flagged[Col.P_EMPIR_FLAG].sum()):,}",
            ),
            (
                "Sites flagged (multiple-testing adjusted)",
                f"{int(flagged[Col.P_EMPIR_ADJ_FLAG].sum()):,}",
            ),
        ]
        return _summary_table_figure(rows)

    def _figures(self):
        return [self._summary_table()]
