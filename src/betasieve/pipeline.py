import pickle

from pathlib import Path

from .analysis import SieveResults, run_duplicate_analysis
from .config import SieveArgs, validate_sieve_args


def _pickle_intermediate_results(args: SieveArgs, results: SieveResults):
    pkl_dir = args.pkl_dir
    pkl_dir.mkdir(parents=True, exist_ok=True)

    for payload, filename in [(args, "args"), (results, "results")]:
        with open(pkl_dir / (filename + ".pkl"), "wb") as file:
            pickle.dump(payload, file)

    print(f"SieveArgs and SieveResults written to {pkl_dir}.")


def _write_csv_outputs(args: SieveArgs, results: SieveResults) -> None:
    csv_dir = args.csv_dir
    csv_dir.mkdir(parents=True, exist_ok=True)

    threshold_label = round(results.threshold, 4)
    results.flagged_frame.to_csv(csv_dir / f"min_max_difference_{threshold_label}.csv")

    if results.sweep_df is not None:
        results.sweep_df.to_csv(csv_dir / "threshold_sweep_summary.csv", index=False)

    if results.candidate_cpgs is not None:
        results.candidate_cpgs.to_csv(csv_dir / "candidate_cpgs.csv", index=False)

    print(f"CSV outputs written to {csv_dir}")


def _write_report(args: SieveArgs, results: SieveResults) -> None:
    from betasieve.report import SieveReportGenerator

    gen = SieveReportGenerator(results, args)
    gen.build_report()


def run_beta_sieve(args: SieveArgs) -> SieveResults:

    validate_sieve_args(args)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    results = run_duplicate_analysis(args)

    if args.pkl:
        _pickle_intermediate_results(args, results)

    _write_csv_outputs(args, results)

    if args.report:
        _write_report(args, results)

    return results


__all__ = ["run_beta_sieve"]
