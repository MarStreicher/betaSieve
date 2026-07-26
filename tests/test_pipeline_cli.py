import pickle
from pathlib import Path

import pandas as pd
import pytest

from betasieve import cli, pipeline
from betasieve.analysis import Col, SieveResults
from betasieve.config import SieveArgs


def test_write_csv_outputs_writes_all_available_frames(
    sieve_args: SieveArgs,
    sieve_results: SieveResults,
) -> None:
    sieve_results.sweep_df = pd.DataFrame({Col.THRESHOLD: [0.1], Col.P0: [0.05]})

    pipeline._write_csv_outputs(sieve_args, sieve_results)

    csv_dir = sieve_args.csv_dir
    assert (csv_dir / "min_max_difference_0.1.csv").is_file()
    assert (csv_dir / "threshold_sweep_summary.csv").is_file()
    assert (csv_dir / "candidate_cpgs.csv").is_file()
    candidate = pd.read_csv(csv_dir / "candidate_cpgs.csv")
    assert candidate["IlmnID"].tolist() == [
        "cg_pair1_TC11",
        "cg_pair1_TC21",
    ]


def test_write_csv_outputs_omits_optional_frames(
    sieve_args: SieveArgs,
    sieve_results: SieveResults,
) -> None:
    sieve_results.sweep_df = None
    sieve_results.candidate_cpgs = None

    pipeline._write_csv_outputs(sieve_args, sieve_results)

    assert (sieve_args.csv_dir / "min_max_difference_0.1.csv").is_file()
    assert not (sieve_args.csv_dir / "threshold_sweep_summary.csv").exists()
    assert not (sieve_args.csv_dir / "candidate_cpgs.csv").exists()


def test_pickle_intermediate_results_round_trips_payloads(
    sieve_args: SieveArgs,
    sieve_results: SieveResults,
) -> None:
    pipeline._pickle_intermediate_results(sieve_args, sieve_results)

    with (sieve_args.pkl_dir / "args.pkl").open("rb") as file:
        restored_args = pickle.load(file)
    with (sieve_args.pkl_dir / "results.pkl").open("rb") as file:
        restored_results = pickle.load(file)

    assert restored_args == sieve_args
    pd.testing.assert_frame_equal(
        restored_results.flagged_frame, sieve_results.flagged_frame
    )


def test_write_report_builds_generator(
    sieve_args: SieveArgs,
    sieve_results: SieveResults,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    observed = {}

    class FakeGenerator:
        def __init__(self, results, args):
            observed["init"] = (results, args)

        def build_report(self):
            observed["built"] = True

    import betasieve.report

    monkeypatch.setattr(betasieve.report, "SieveReportGenerator", FakeGenerator)

    pipeline._write_report(sieve_args, sieve_results)

    assert observed == {
        "init": (sieve_results, sieve_args),
        "built": True,
    }


@pytest.mark.parametrize(
    ("write_pickle", "write_report"),
    [(False, False), (True, True)],
)
def test_run_beta_sieve_orchestrates_optional_outputs(
    sieve_args: SieveArgs,
    sieve_results: SieveResults,
    monkeypatch: pytest.MonkeyPatch,
    write_pickle: bool,
    write_report: bool,
) -> None:
    sieve_args.pkl = write_pickle
    sieve_args.report = write_report
    calls = []

    monkeypatch.setattr(
        pipeline, "validate_sieve_args", lambda args: calls.append("validate")
    )
    monkeypatch.setattr(
        pipeline,
        "run_duplicate_analysis",
        lambda args: calls.append("analyze") or sieve_results,
    )
    monkeypatch.setattr(
        pipeline,
        "_pickle_intermediate_results",
        lambda args, results: calls.append("pickle"),
    )
    monkeypatch.setattr(
        pipeline,
        "_write_csv_outputs",
        lambda args, results: calls.append("csv"),
    )
    monkeypatch.setattr(
        pipeline,
        "_write_report",
        lambda args, results: calls.append("report"),
    )

    result = pipeline.run_beta_sieve(sieve_args)

    expected = ["validate", "analyze"]
    if write_pickle:
        expected.append("pickle")
    expected.append("csv")
    if write_report:
        expected.append("report")
    assert calls == expected
    assert result is sieve_results
    assert sieve_args.out_dir.is_dir()


def test_parser_defaults_to_automatic_threshold_search(tmp_path: Path) -> None:
    betas = tmp_path / "betas.tsv"
    namespace = cli.build_parser().parse_args(["--betas", str(betas)])

    assert namespace.betas == betas
    assert namespace.threshold is None
    assert namespace.threshold_min == 0.01
    assert namespace.threshold_max == 0.1
    assert namespace.threshold_step == 0.01
    assert namespace.report is True
    assert namespace.pkl is False


def test_parser_accepts_boolean_and_analysis_options(tmp_path: Path) -> None:
    betas = tmp_path / "betas.tsv"
    namespace = cli.build_parser().parse_args(
        [
            "--betas",
            str(betas),
            "--threshold",
            "0.2",
            "--confidence",
            "0.9",
            "--fdr",
            "holm",
            "--no-report",
            "--pkl",
        ]
    )

    assert namespace.threshold == 0.2
    assert namespace.confidence == 0.9
    assert namespace.fdr == "holm"
    assert namespace.report is False
    assert namespace.pkl is True


def test_main_converts_namespace_and_runs_pipeline(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    betas = tmp_path / "betas.tsv"
    observed = []
    monkeypatch.setattr(cli, "run_beta_sieve", observed.append)

    cli.main(["--betas", str(betas), "--threshold", "0.2", "--no-report"])

    assert len(observed) == 1
    assert observed[0].betas_path == betas
    assert observed[0].threshold == 0.2
    assert observed[0].report is False
