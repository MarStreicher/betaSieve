from __future__ import annotations

import base64
from datetime import datetime
from functools import lru_cache
from pathlib import Path
from typing import Callable, List, Optional, TYPE_CHECKING

from betasieve.columns import Col
from betasieve.config import PipelineConfig
from betasieve.report.report_section import ReportMainSection
from betasieve.report.sections.experiment_config import ConfigSection
from betasieve.report.sections.differences import DifferencesSection
from betasieve.report.sections.threshold_sweep import ThresholdSweepSection
from betasieve.report.sections.p_hat_distribution import PhatDistributionSection
from betasieve.report.sections.evaluation import EvaluationSection
from betasieve.report.sections.output_description import OutputDescriptionSection

if TYPE_CHECKING:
    from betasieve.analysis import SieveResults


_SECTION_REGISTRY: List[
    Callable[["SieveResults", PipelineConfig], ReportMainSection]
] = [
    ConfigSection,
    DifferencesSection,
    ThresholdSweepSection,
    PhatDistributionSection,
    EvaluationSection,
    OutputDescriptionSection,
]


_RESOURCES_DIR = Path(__file__).resolve().parent / "resources"

_FONT_FILES = {
    "resources/fonts/Montserrat-VariableFont_wght.ttf": "fonts/Montserrat-VariableFont_wght.ttf",
    "resources/fonts/Montserrat-Italic-VariableFont_wght.ttf": "fonts/Montserrat-Italic-VariableFont_wght.ttf",
}


@lru_cache(maxsize=1)
def _resources_dir() -> Path:
    if not _RESOURCES_DIR.is_dir():
        raise FileNotFoundError(
            f"Report resources directory not found: {_RESOURCES_DIR}"
        )
    return _RESOURCES_DIR


def _resource_path(name: str) -> Path:
    return _resources_dir() / name


def _read_resource(filename: str) -> str:
    path = _resources_dir() / filename
    return path.read_text(encoding="utf-8")


def _file_to_data_uri(path: Path, mime: str) -> str:
    b64 = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:{mime};base64,{b64}"


def _embed_image(filename: str, mime: str) -> str:
    path = _resources_dir() / "images" / filename
    if not path.is_file():
        return ""
    raw = path.read_bytes()
    b64 = base64.b64encode(raw).decode("ascii")
    return f"data:{mime};base64,{b64}"


def _embed_logo() -> str:
    for filename, mime in (
        ("betasieve_logo_transparent.png", "image/png"),
        ("betasieve_logo_transparent.svg", "image/svg+xml"),
    ):
        uri = _embed_image(filename, mime)
        if uri:
            return uri
    return ""


class SieveReportGenerator:
    def __init__(self, results: "SieveResults", args: PipelineConfig) -> None:
        self.results = results
        self.args = args

    def _generate_style(self) -> str:
        css = _read_resource("template.css")
        for src, resource_name in _FONT_FILES.items():
            css = css.replace(
                src, _file_to_data_uri(_resource_path(resource_name), "font/ttf")
            )
        return f"<style>{css}</style>"

    def _generate_preamble(self) -> str:
        r = self.results
        flagged = r.flagged_frame

        n_sites = len(flagged)
        p0_val = float(flagged[Col.P0].dropna().iloc[0])
        confidence = float(flagged[Col.CONFIDENCE].iloc[0])
        alpha = 1.0 - confidence
        threshold_val = float(flagged[Col.THRESHOLD].iloc[0])
        threshold_source = (
            "automatically selected by threshold sweep"
            if r.sweep_df is not None
            else "user-specified"
        )
        n_flagged = int(flagged[Col.P_EMPIR_ADJ_FLAG].sum())
        pct_flagged = 100.0 * n_flagged / n_sites if n_sites else 0.0

        n_flagged_probes = int(len(r.candidate_cpgs))

        return (
            f"<p>This report documents the <strong>betaSieve</strong> analysis "
            f"of cg-only probes "
            f"on data from the "
            f"<strong>Infinium&#x2122; MethylationEPIC v2.0 BeadChip (EPICv2)</strong>. "
            f"Please refer to [COMING SOON] for further information.</p>"
            f"<p><strong>Results:</strong> "
            f"In total, <strong>{n_sites:,}</strong> CpG-sites were analysed. "
            f"A difference threshold of <em>t</em>&nbsp;=&nbsp;{threshold_val:.2f} "
            f"({threshold_source}) was applied and the empirical background exceedance "
            f"rate at exact-replicate CpG-sites was "
            f"<em>p&#x2080;</em>&nbsp;=&nbsp;{p0_val:.2f}. "
            f"Statistical significance was assessed at "
            f"<em>&#x03B1;</em>&nbsp;=&nbsp;{alpha:.2f}, with FDR "
            f"correction applied across design replicate CpG-sites. "
            f"<strong>{n_flagged:,}</strong> CpG-sites ({pct_flagged:.1f}&#x25;), therefore, "
            f"<strong>{n_flagged_probes}</strong> probes, were "
            f"flagged as discordant based on the multiple-testing-adjusted empirical "
            f"upper-tail p-value.</p>"
            f"<p><strong>Definitions:</strong> "
            f"CpG-sites were analysed in two replicate classes: "
            f"<strong>exact replicates</strong> and <strong>design replicates</strong>. "
            f"Design replicates comprise four subgroups of sites with different probe designs: "
            f"<strong>pair type</strong> (two designs that differ only in probe type, "
            f"TC1 and TC2, or BC1 and BC2), "
            f"<strong>pair design</strong> (all remaining sites with exactly two designs), "
            f"<strong>triplets</strong>, and <strong>quadruplets</strong> "
            f"(three or four different designs, respectively). "
            f"<strong>Exact replicates</strong> are technical replicates of the same "
            f"probe design.</p>"
        )

    def _generate_toc(
        self, sections: List[ReportMainSection], prefix: str = "", level: int = 0
    ) -> str:
        toc = ""
        for i, section in enumerate(sections, start=1):
            number = f"{prefix}{i}"
            toc += (
                f'<details class="toc-group">'
                f'  <summary class="toc-group-title">'
                f'    <a class="toc-item level-0" href="#{section.id}">'
                f"      {number} {section.title}"
                f"    </a>"
                f"  </summary>"
                f'  <div class="toc-subitems">'
            )
            for j, sub in enumerate(section.subsections, start=1):
                sub_number = f"{number}.{j}"
                toc += (
                    f'<a class="toc-item level-1" href="#{sub.id}">'
                    f"  {sub_number} {sub.title}"
                    f"</a>\n"
                )
            toc += "</div></details>\n"
        return toc

    def _generate_content(self, sections: List[ReportMainSection]) -> str:
        html = ""
        for i, section in enumerate(sections, start=1):
            html += (
                f'<div id="{section.id}" class="section level-0">'
                f"  <h2>{i} {section.title}</h2>"
                f"  <p>{section.description}</p>"
                f"</div>\n"
            )
            for j, sub in enumerate(section.subsections, start=1):
                sub.generate()
                figures_html = ""
                for fig in sub.figures:
                    if isinstance(fig, str):
                        figures_html += fig
                    else:
                        figures_html += fig.to_html(
                            full_html=False, include_plotlyjs=False
                        )
                if figures_html:
                    html += (
                        f'<div id="{sub.id}" class="section level-1">'
                        f"  <h3>{i}.{j} {sub.title}</h3>"
                        f"  <p>{sub.description}</p>"
                        f"  {figures_html}"
                        f"</div>\n"
                    )
        return html

    def build_report(self, output_path: Optional[Path] = None) -> Path:
        output_path = Path(output_path or self.args.report_dir)

        sections = [
            section_cls(self.results, self.args) for section_cls in _SECTION_REGISTRY
        ]

        style = self._generate_style()
        toc = self._generate_toc(sections)
        preamble = self._generate_preamble()
        content = self._generate_content(sections)
        logo_uri = _embed_logo()

        template = _read_resource("template.html")
        html = template.format(
            tab_title="betaSieve Report",
            report_title="Replicate Analysis",
            style=style,
            logo_uri=logo_uri,
            toc_items=toc,
            preamble=preamble,
            sections=content,
            generated=datetime.now().strftime("%Y-%m-%d %H:%M"),
        )

        out_file = output_path.parent / (output_path.name + ".html")
        out_file.parent.mkdir(parents=True, exist_ok=True)
        out_file.write_text(html, encoding="utf-8")
        print(f"Report written to {out_file}")
        return out_file
