from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, List, final

if TYPE_CHECKING:
    from betasieve.analysis import SieveResults
    from betasieve.config import SieveArgs


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

    def __init__(self, results: "SieveResults", args: "SieveArgs") -> None:
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

    def __init__(self, results: "SieveResults", args: "SieveArgs") -> None:
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
    A leaf section that owns figures.  Implement ``generate()`` to
    populate ``self.figures`` with Plotly Figure objects.
    """

    @property
    @final
    def subsections(self):
        raise AttributeError(
            f"{type(self).__name__} is a subsection and must not define subsections."
        )

    def __init__(self, results: "SieveResults", args: "SieveArgs") -> None:
        super().__init__(results, args)
        self.figures: list = []

    @abstractmethod
    def generate(self) -> None:
        """Populate ``self.figures`` from ``self.results`` and ``self.args``."""
        ...
