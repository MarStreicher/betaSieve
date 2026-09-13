from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, List, Sequence, final

if TYPE_CHECKING:
    from betasieve.analysis import SieveResults
    from betasieve.config import PipelineConfig


class ReportSection(ABC):
    """Base class for all report sections."""

    @property
    @abstractmethod
    def title(self) -> str: ...

    @property
    def description(self) -> str:
        return ""

    @property
    def id(self) -> str:
        return self.title.lower().replace(" ", "-")

    def __init__(self, results: "SieveResults", args: "PipelineConfig") -> None:
        self.results = results
        self.args = args


class ReportMainSection(ReportSection):
    """
    A section that groups subsections.  It has no figures of its own
    and must not implement ``generate()``.
    """

    @property
    @abstractmethod
    def subsection_types(self) -> List[type["ReportSubSection"]]: ...

    def __init__(self, results: "SieveResults", args: "PipelineConfig") -> None:
        super().__init__(results, args)
        self.subsections: List[ReportSubSection] = [
            cls(results, args) for cls in self.subsection_types
        ]

    @final
    def generate(self) -> None:
        raise NotImplementedError(
            f"{type(self).__name__} is a main section — call generate() on its subsections."
        )


class ReportSubSection(ReportSection):
    """
    A leaf section that owns figures.

    Implement ``_figures()`` to return the Plotly figures (or HTML snippets)
    for this subsection. ``generate()`` appends them to ``self.figures``.
    """

    @property
    @final
    def subsections(self):
        raise AttributeError(
            f"{type(self).__name__} is a subsection and must not define subsections."
        )

    def __init__(self, results: "SieveResults", args: "PipelineConfig") -> None:
        super().__init__(results, args)
        self.figures: list = []

    @abstractmethod
    def _figures(self) -> Sequence[object]: ...

    @final
    def generate(self) -> None:
        self.figures.extend(self._figures())
