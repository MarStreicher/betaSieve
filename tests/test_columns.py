import pandas as pd

from betasieve.columns import Col, _diff_value_columns


def test_diff_value_columns_returns_only_numeric_sample_columns(
    diff_frame: pd.DataFrame,
) -> None:
    frame = diff_frame.copy()
    frame[Col.N] = 4
    frame["note"] = "text"

    assert _diff_value_columns(frame) == [
        "Sample_A",
        "Sample_B",
        "Sample_C",
        "Sample_D",
    ]
