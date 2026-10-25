from typing import List


def raise_validation_errors(title: str, errors: List[str]) -> None:
    """Raise a ``ValueError`` listing every collected error as a bullet."""
    if errors:
        message = f"{title}:\n" + "\n".join(f"  • {err}" for err in errors)
        raise ValueError(message)
