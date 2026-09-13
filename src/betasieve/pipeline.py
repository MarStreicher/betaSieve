import pickle

from epicv2io import BetasLoader
from .analysis import SieveResults, sieve_betas
from .config import PipelineConfig, validate_pipeline_config


def _pickle_intermediate_results(config: PipelineConfig, results: SieveResults):
    pkl_dir = config.pkl_dir
    pkl_dir.mkdir(parents=True, exist_ok=True)

    for payload, filename in [(config, "args"), (results, "results")]:
        with open(pkl_dir / (filename + ".pkl"), "wb") as file:
            pickle.dump(payload, file)

    print(f"PipelineConfig and SieveResults written to {pkl_dir}.")


def _write_csv_outputs(config: PipelineConfig, results: SieveResults) -> None:
    csv_dir = config.csv_dir
    csv_dir.mkdir(parents=True, exist_ok=True)

    threshold_label = round(results.threshold, 4)
    results.flagged_frame.to_csv(csv_dir / f"min_max_difference_{threshold_label}.csv")

    if results.sweep_df is not None:
        results.sweep_df.to_csv(csv_dir / "threshold_sweep_summary.csv", index=False)

    if results.candidate_cpgs is not None:
        results.candidate_cpgs.to_csv(csv_dir / "candidate_cpgs.csv", index=False)

    print(f"CSV outputs written to {csv_dir}")


def _write_report(config: PipelineConfig, results: SieveResults) -> None:
    from betasieve.report import SieveReportGenerator

    gen = SieveReportGenerator(results, config)
    gen.build_report()


def run_sieve_pipeline(config: PipelineConfig) -> SieveResults:
    """Load betas from disk, run sieve_betas, and write any requested outputs."""

    cg_by_sample = BetasLoader(config.betas_path).load_data()
    results = sieve_betas(cg_by_sample, config.analysis)

    if config.pkl:
        _pickle_intermediate_results(config, results)

    if config.csv_files:
        _write_csv_outputs(config, results)

    if config.report:
        _write_report(config, results)

    return results


__all__ = ["run_sieve_pipeline"]
