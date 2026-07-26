from typing import List

from plotly.graph_objects import Figure

from betasieve.analysis import Col
from betasieve.cg_probe_table import DesignGroup
from betasieve.report.report_section import ReportMainSection, ReportSubSection
from betasieve.report.tables import _summary_table_figure


class ConfigSection(ReportMainSection):
    @property
    def title(self) -> str:
        return "Analysis Configuration"

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
            "The table summarizes the analysis parameters used by betaSieve. "
            "You can specify the threshold t directly or determine it automatically through a threshold sweep. "
            "If you provide t, the empirical background exceedance rate p₀ is calculated as the fraction of sites in the exact replicates that exceed this threshold. "
            "Alternatively, you can provide a threshold range and a target p₀. "
            "betaSieve then performs a threshold sweep and selects the threshold t that most closely matches your target background exceedance rate."
        )

    def _summary_table(self) -> Figure:
        args = self.args

        rows = [
            ("Betas file", str(args.betas_path)),
            (
                ("Threshold", str(args.threshold))
                if args.threshold is not None
                else ("Threshold", "not used (threshold sweep)")
            ),
            ("Confidence", str(args.confidence)),
            ("FDR method", args.fdr),
            (
                "Target p₀ (threshold sweep)",
                (
                    str(args.target_p0)
                    if args.threshold is None
                    else "not used (fixed threshold)"
                ),
            ),
            (
                "Threshold sweep min / max / step",
                (
                    f"{args.threshold_min} / {args.threshold_max} / {args.threshold_step}"
                    if args.threshold_min is not None
                    else "not used (fixed threshold)"
                ),
            ),
        ]
        return _summary_table_figure(rows)

    def generate(self) -> None:
        self.figures.append(self._summary_table())


class ConfigResultsSubSection(ReportSubSection):
    @property
    def title(self) -> str:
        return "Dataset Summary"

    @property
    def description(self) -> str:
        return (
            "Key statistics derived from the input data and analysis results. "
            "Site counts are reported per design group. Please note that individual CpG "
            "sites may appear in more than one group (e.g. a site belonging to both "
            "a pair and to exact replicates). "
            "The selected threshold t and corresponding p₀ are the values used for "
            "all downstream statistical tests and figures."
        )

    def _summary_table(self) -> Figure:
        args = self.args
        flagged = self.results.flagged_frame
        group = flagged[Col.GROUP]

        rows = [
            ("Number of samples (n)", int(flagged[Col.N].iloc[0])),
            ("Total CpG sites analysed", f"{len(flagged):,}"),
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
            ("FDR correction method", args.fdr),
            (
                "Selected threshold t",
                round(float(flagged[Col.THRESHOLD].iloc[0]), 4),
            ),
            (
                "Background exceedance rate p₀",
                round(float(flagged[Col.P0].iloc[0]), 4),
            ),
            (
                "Sites flagged",
                f"{int(flagged[Col.P_BETA_FLAG].sum()):,}",
            ),
            (
                "Sites flagged (FDR-corrected)",
                f"{int(flagged[Col.P_BETA_ADJ_FLAG].sum()):,}",
            ),
        ]
        return _summary_table_figure(rows)

    def generate(self) -> None:
        self.figures.append(self._summary_table())
