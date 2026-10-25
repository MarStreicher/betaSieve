import os
from pathlib import Path
from dotenv import load_dotenv

from betasieve import PipelineConfig, SieveConfig, run_sieve_pipeline

ROOT = Path(__file__).resolve().parent.parent
if not load_dotenv(ROOT / ".env"):
    raise FileNotFoundError(
        "Missing .env. Copy .env.example to .env and set the local paths."
    )

betas_path = Path(os.environ["BETAS_PATH"])
out_dir = Path(os.environ.get("OUT_DIR", ROOT / "results"))

pipeline_config = PipelineConfig(
    betas_path=betas_path,
    out_dir=out_dir,
    analysis=SieveConfig(threshold=0.01),
    report=True,
    pkl=False,
)

results = run_sieve_pipeline(pipeline_config)
n_candidates = 0 if results.candidate_cpgs is None else len(results.candidate_cpgs)
print(f"Threshold: {results.threshold}, candidates: {n_candidates}")
