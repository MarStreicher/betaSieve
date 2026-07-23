<p align="right">
  <img src="assets/betasieve_logo.png" alt="betaSieve logo" width="250">
</p>

<h1 align="left">betaSieve</h1>

<p align="left">
  Filter-first analysis of EPICv2 duplicate CpG probes
</p>

betaSieve is a Python package for identifying probes from Illumina EPICv2 methylation arrays that exhibit high variability between duplicate measurements. To do so, the package evaluates duplicate probe agreement by analyzing beta-value differences for each (CpG site, sample) pair.

---

## Features

- Analysis of EPICv2 duplicate probes
- Optional threshold sweep across user-defined beta-value differences
- Identification of highly variable candidate probes

---

## Installation

### From GitHub

```bash
pip install betasieve
```

### Requirements

- Python ≥ 3.9
- NumPy
- Pandas
- SciPy
- Statsmodels
- Matplotlib
- Plotly

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

## Data and Annotation Resources

betaSieve depends on [EpicV2IO](https://github.com/MarStreicher/EpicV2IO) for the
shipped Peters et al. EPICv2 probe-annotation subset (Parquet). That file is a
**derived column subset** of third-party open materials (not original betaSieve
content). See EpicV2IO’s [`NOTICE`](https://github.com/MarStreicher/EpicV2IO/blob/main/NOTICE)
for full provenance.

When you use betaSieve analyses that rely on that annotation, cite Peters et al.
(2024) as in the [Citation](#citation) section below.

## Citation

If you use betaSieve in published work, please cite:

> Streicher M. betaSieve: Filter-first analysis of EPICv2 duplicate CpG probes. GitHub repository. 2026.

Please also cite the EPICv2 annotation resource used via EpicV2IO (CC BY 4.0):

> Peters TJ, Meyer B, Ryan L, Achinger-Kawecka J, Song J, Campbell EM, Qu W, Nair S, Loi-Luu P, Stricker P, Lim E, Stirzaker C, Clark SJ, Pidsley R. Characterisation and reproducibility of the HumanMethylationEPIC v2.0 BeadChip for DNA methylation profiling. *BMC Genomics*. 2024;25:251. https://doi.org/10.1186/s12864-024-10027-5

---

## License

BSD 3-Clause License (see [`LICENSE`](LICENSE)) covers **this software only**.

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

Method used for multiple-testing correction.

Default:

```python
fdr="fdr_bh"
```

Available methods:

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
