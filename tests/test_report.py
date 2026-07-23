from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import plotly.graph_objects as go
import pytest

from betasieve.analysis import Col, SieveResults
from betasieve.config import SieveArgs
from betasieve.report.figure_style import configure_matplotlib
from betasieve.report.plots import _layout_figure
from betasieve.report.report_generator import (
    SieveReportGenerator,
    _embed_image,
    _embed_logo,
    _read_resource,
    _resources_dir,
)
from betasieve.report.report_section import ReportMainSection, ReportSubSection
from betasieve.report.sections.output_description import OutputDescriptionSection
from betasieve.report.sections.threshold_sweep import ThresholdSweepSection
from betasieve.report.tables import _data_dict_figure, _summary_table_figure


class ExampleSubSection(ReportSubSection):
    @property
    def title(self) -> str:
        return "Example Child"

    def generate(self) -> None:
        self.figures.append("<p>generated</p>")


class ExampleMainSection(ReportMainSection):
    @property
    def title(self) -> str:
        return "Example Main"

    @property
    def subsection_types(self):
        return [ExampleSubSection]


def test_report_section_ids_and_tree_contracts(
    sieve_args: SieveArgs, sieve_results: SieveResults
) -> None:
    section = ExampleMainSection(sieve_results, sieve_args)
    child = section.subsections[0]

    assert section.id == "example-main"
    assert child.id == "example-child"
    assert child.description == ""
    with pytest.raises(NotImplementedError, match="main section"):
        section.generate()
    with pytest.raises(AttributeError, match="must not define subsections"):
        _ = child.subsections


def test_table_factories_encode_values_and_size() -> None:
    summary = _summary_table_figure([("Alpha", 1), ("Beta", "two")])
    data_dict = _data_dict_figure(
        ["Name", "Meaning"],
        [("x", "first"), ("y", "second")],
        col_widths=[0.3, 0.7],
        row_height=40,
    )

    assert summary.data[0].cells.values == (
        ["Alpha", "Beta"],
        ["1", "two"],
    )
    assert summary.layout.height == 104
    assert data_dict.data[0].header.values == ("Name", "Meaning")
    assert data_dict.data[0].columnwidth == (0.3, 0.7)
    assert data_dict.layout.height == 128


def test_layout_figure_sets_shared_chart_style() -> None:
    figure = _layout_figure(
        go.Figure(),
        title="Title",
        x_title="X",
        y_title="Y",
        height=321,
        show_legend=False,
    )

    assert figure.layout.title.text == "Title"
    assert figure.layout.xaxis.title.text == "X"
    assert figure.layout.yaxis.title.text == "Y"
    assert figure.layout.height == 321
    assert figure.layout.showlegend is False


def test_configure_matplotlib_applies_shared_defaults() -> None:
    configure_matplotlib()

    assert plt.rcParams["font.family"] == ["serif"]
    assert plt.rcParams["mathtext.fontset"] == "cm"
    assert plt.rcParams["axes.unicode_minus"] is False


def test_resource_helpers_read_text_and_embed_images(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    resources = tmp_path / "resources"
    images = resources / "images"
    images.mkdir(parents=True)
    (resources / "template.css").write_text("body {}", encoding="utf-8")
    (images / "logo_transparent.svg").write_bytes(b"<svg/>")

    import betasieve.report.report_generator as report_generator

    monkeypatch.setattr(report_generator, "_RESOURCES_DIR", resources)
    _resources_dir.cache_clear()
    try:
        assert _read_resource("template.css") == "body {}"
        assert _embed_image("missing.png", "image/png") == ""
        assert _embed_logo() == "data:image/svg+xml;base64,PHN2Zy8+"
    finally:
        _resources_dir.cache_clear()


def test_resources_dir_raises_for_missing_directory(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import betasieve.report.report_generator as report_generator

    monkeypatch.setattr(report_generator, "_RESOURCES_DIR", tmp_path / "missing")
    _resources_dir.cache_clear()
    try:
        with pytest.raises(FileNotFoundError, match="resources directory not found"):
            _resources_dir()
    finally:
        _resources_dir.cache_clear()


def test_dynamic_sections_reflect_optional_results(
    sieve_args: SieveArgs, sieve_results: SieveResults
) -> None:
    fixed = ThresholdSweepSection(sieve_results, sieve_args)
    outputs = OutputDescriptionSection(sieve_results, sieve_args)

    assert [section.__name__ for section in fixed.subsection_types] == [
        "FixedThresholdPlaceholderSubSection",
        "ThresholdSummaryTableSubSection",
    ]
    assert [section.__name__ for section in outputs.subsection_types] == [
        "MinMaxDiffCsvSubSection",
        "CandidateCpgsCsvSubSection",
    ]

    sieve_results.sweep_df = pd.DataFrame(
        {
            Col.GROUP: ["Exact replicates"],
            Col.THRESHOLD: [0.1],
            Col.P0: [0.05],
            Col.PCT_BETA_ADJ_FLAGGED: [0.0],
        }
    )
    sweep = ThresholdSweepSection(sieve_results, sieve_args)
    outputs = OutputDescriptionSection(sieve_results, sieve_args)
    assert [section.__name__ for section in sweep.subsection_types] == [
        "P0SubSection",
        "BetaSection",
        "ThresholdSummaryTableSubSection",
    ]
    assert "SweepSummaryCsvSubSection" in [
        section.__name__ for section in outputs.subsection_types
    ]


@pytest.mark.integration
def test_build_report_writes_standalone_html(
    tmp_path: Path,
    sieve_args: SieveArgs,
    sieve_results: SieveResults,
) -> None:
    output = tmp_path / "report" / "betasieve"

    result = SieveReportGenerator(sieve_results, sieve_args).build_report(output)

    assert result == output.with_suffix(".html")
    html = result.read_text(encoding="utf-8")
    assert "<title>betaSieve Report</title>" in html
    assert "Analysis Configuration" in html
    assert "Probe Range Analysis" in html
    assert "Output Files" in html
    assert "<strong>2</strong> probes" in html
    assert "data:image/" in html
