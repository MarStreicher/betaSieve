<p align="right">
  <img src="https://raw.githubusercontent.com/MarStreicher/betaSieve/main/assets/betasieve_logo.png" alt="betaSieve logo" width="250">
</p>

<h1 align="left">betaSieve</h1>

<p align="left">
  EPICv2 probe designs analysis
</p>

betaSieve is a Python package that identifies HumanMethylationEPIC v2.0 BeadChip (EPICv2) probes that exhibit high variability between measurements of different probe designs. To do so, the package evaluates probe design agreement by analyzing beta-value differences for each (CpG site, sample) pair.

## Why should I care?

EPICv2 can measure the **same CpG site with several probe designs**. Those designs appear as separate `IlmnID` rows, and we observed that their β-values can disagree for the same sample. 

<p align="center">
  <img src="https://raw.githubusercontent.com/MarStreicher/betaSieve/main/assets/epicv2_ilmnid_problem.svg" alt="EPICv2 IlmnID naming: one CpG site, multiple probe designs, and suffix encoding" width="780">
</p>

<p align="center"><em>
IlmnID naming and replicate definitions adapted from Peters et al. (2024),
BMC Genomics 25:251, and Supplementary File 4
(<a href="https://doi.org/10.1186/s12864-024-10027-5">doi:10.1186/s12864-024-10027-5</a>).
</em></p>

betaSieve finds sites where design disagreement is larger than expected from exact technical replicates. It can thus help to refine downstream analyses. 

Please have a look at the [wiki](https://github.com/MarStreicher/betaSieve/wiki) for a more detailed explanation of the background and theory behind this package.

---

## Features

1. Analysis of EPICv2 IlmnIDs with different probe designs.
2. Optional threshold sweep across user-defined beta-value differences.
3. Identification of highly variable candidate IlmnIDs that should be flagged.

---

## Installation

### From GitHub

```bash
pip install betasieve
```

Dependencies are installed automatically.

---

## Quick Start

### Python API

```python
from pathlib import Path
from betasieve import SieveArgs, run_beta_sieve

args = SieveArgs(
    betas_path=Path("betas.csv"),
    threshold_min=0.03,
    threshold_max=0.07,
    threshold_step=0.01,
)

results = run_beta_sieve(args)

n_candidates = (
    0 if results.candidate_cpgs is None
    else len(results.candidate_cpgs)
)

print(
    f"Threshold: {results.threshold}, "
    f"candidates: {n_candidates}"
)
```

### Command Line

```bash
betasieve \
  --betas betas.csv \
  --threshold_min 0.03 \
  --threshold_max 0.07 \
  --threshold_step 0.01
```

---

## Citation

If you use betaSieve in published work, please cite:

> Streicher M. betaSieve: Filter-first analysis of EPICv2 duplicate CpG probes. GitHub repository. 2026.

For EPICv2 probe naming and replicate definitions (as used in the overview figure above), please also cite:

> Peters TJ, Meyer B, Ryan L, Achinger-Kawecka J, Song J, Campbell EM, Qu W, Nair S, Loi-Luu P, Stricker P, Lim E, Stirzaker C, Clark SJ, Pidsley R. Characterisation and reproducibility of the HumanMethylationEPIC v2.0 BeadChip for DNA methylation profiling. *BMC Genomics*. 2024;25:251. https://doi.org/10.1186/s12864-024-10027-5

---

## License

BSD 3-Clause License (see [`LICENSE`](LICENSE)).

---

## Parameters

betaSieve can be configured either through the Python API (`SieveArgs`) or through the command-line interface.

### Required Parameters

#### `betas_path`

Path to the input beta-value matrix.

Expected format:

- Rows correspond to CpG probes.
- Columns correspond to samples.
- Values must be beta values between 0 and 1.

Example:

```python
betas_path=Path("betas.csv")
```

---

### Threshold Selection

betaSieve requires either:

1. A fixed threshold (`threshold`)
2. An automatic threshold search (`threshold_min`, `threshold_max`, `threshold_step`)

#### `threshold`

Fixed beta-value difference threshold used to determine whether duplicate probes are considered concordant.

Range:

```text
(0, 1]
```

Example:

```python
threshold=0.05
```

---

#### `threshold_min`

Minimum threshold evaluated during automatic threshold search.

Example:

```python
threshold_min=0.03
```

---

#### `threshold_max`

Maximum threshold evaluated during automatic threshold search.

Example:

```python
threshold_max=0.07
```

---

#### `threshold_step`

Increment between tested thresholds.

Example:

```python
threshold_step=0.01
```

The above configuration evaluates:

```text
0.03, 0.04, 0.05, 0.06, 0.07
```

---

### Statistical Parameters

#### `fdr`

Method used for multiple-testing correction. betaSieve does not implement these procedures itself; it passes the chosen method to [`statsmodels.stats.multitest.multipletests`](https://www.statsmodels.org/stable/generated/statsmodels.stats.multitest.multipletests.html) from the [statsmodels](https://www.statsmodels.org/) package ([source](https://github.com/statsmodels/statsmodels/blob/main/statsmodels/stats/multitest.py)).

Default:

```python
fdr="fdr_bh"
```

Available methods (as supported by `statsmodels.stats.multitest.multipletests`):

| Method         | Description                           |
| -------------- | ------------------------------------- |
| bonferroni     | Bonferroni correction                 |
| sidak          | Sidak correction                      |
| holm-sidak     | Holm-Sidak procedure                  |
| holm           | Holm procedure                        |
| simes-hochberg | Simes-Hochberg procedure              |
| hommel         | Hommel procedure                      |
| fdr_bh         | Benjamini-Hochberg FDR                |
| fdr_by         | Benjamini-Yekutieli FDR               |
| fdr_tsbh       | Two-stage Benjamini-Hochberg          |
| fdr_tsbky      | Two-stage Benjamini-Krieger-Yekutieli |

---

#### `confidence`

Confidence level used for confidence interval calculations.

Default:

```python
confidence=0.95
```

Typical values:

| Value | Interpretation          |
| ----- | ----------------------- |
| 0.90  | 90% confidence interval |
| 0.95  | 95% confidence interval |
| 0.99  | 99% confidence interval |

---

#### `target_p0`

Target false-positive rate used internally during threshold optimization.

Default:

```python
target_p0=0.05
```

---

### Output Parameters

#### `out_dir`

Directory where all results are written.

Default:

```python
out_dir=Path("results")
```

Generated structure:

```text
results/
├── csv/
├── figures/
├── report/
└── pkl/
```

---

#### `report`

Generate an HTML summary report.

Default:

```python
report=True
```

Disable report generation:

```python
report=False
```

---

#### `pkl`

Save intermediate Python objects as pickle files.

Default:

```python
pkl=False
```

---

## Examples

### Fixed Threshold

```python
args = SieveArgs(
    betas_path=Path("betas.csv"),
    threshold=0.05,
)
```

### Automatic Threshold Search

```python
args = SieveArgs(
    betas_path=Path("betas.csv"),
    threshold_min=0.03,
    threshold_max=0.07,
    threshold_step=0.01,
)
```
