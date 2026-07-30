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
3. Cohort-specific calibration against exact technical replicates.
4. Empirical upper-tail testing with multiple-testing correction.
5. An HTML report comparing empirical, binomial, and beta-binomial null distributions.
6. Identification of highly variable candidate IlmnIDs that can be reviewed for masking.

---

## Method overview

### Probe-design groups

betaSieve parses EPICv2 IlmnIDs and classifies duplicate measurements into:

- **Pair type:** two designs from the same chemistry family (`TC1`/`TC2` or `BC1`/`BC2`);
- **Pair design:** any other site with exactly two design types;
- **Triplet:** three design types;
- **Quadruplet:** four design types;
- **Exact replicates:** repeated measurements using the same design.

These categories are not mutually exclusive. A site can contribute both a design-comparison row and an exact-replicate row.

### Max–min ranges and threshold calibration

For every CpG-site–sample pair, betaSieve calculates the max–min β-value range across the probes in the relevant group. Given a threshold \(t\), an exceedance is recorded when this range is strictly greater than \(t\).

The empirical background exceedance rate is

\[
p_0 =
\frac{\text{exact-replicate site–sample pairs with range}>t}
{\text{all exact-replicate site–sample pairs}}.
\]

When a threshold sweep is requested, betaSieve selects the smallest evaluated threshold for which \(p_0\) is less than or equal to `target_p0`. If none of the evaluated thresholds reaches the target, it selects the threshold with the smallest observed \(p_0\).

### Site-level exceedance rate

For site \(i\), betaSieve calculates

\[
\hat p_i =
\frac{\text{number of samples with range}>t}
{\text{number of samples}}.
\]

The primary empirical upper-tail p-value compares \(\hat p_i\) with the observed exact-replicate reference distribution:

\[
p_i^{\mathrm{emp}} =
\frac{1 + \#\{j:\hat p_j^{ER}\geq\hat p_i\}}
{m+1},
\]

where \(m\) is the number of exact-replicate reference sites. The plus-one correction prevents zero p-values. Multiple-testing correction is applied across non-exact-replicate comparisons, and `p_empir_adj_flagged` is the current primary flag.

The output also contains a one-sided normal-approximation z-test and a Wilson lower-bound criterion for reference. The HTML report overlays the observed exact-replicate distribution with binomial and method-of-moments beta-binomial models. These modelled distributions are diagnostic and do not currently determine the candidate list.

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
  --threshold-min 0.03 \
  --threshold-max 0.07 \
  --threshold-step 0.01
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

Fixed β-value max–min threshold. A site–sample observation is an exceedance when its range is strictly greater than this value.

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

Method used for multiple-testing correction. Despite the parameter name, the supported choices include both false-discovery-rate and family-wise-error-rate procedures. betaSieve passes the selected method to [`statsmodels.stats.multitest.multipletests`](https://www.statsmodels.org/stable/generated/statsmodels.stats.multitest.multipletests.html).

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

Defines the testing level as `alpha = 1 - confidence` and supplies the normal quantile used by the Wilson calculation. With the default value of `0.95`, the implementation uses a one-sided 95% lower-bound criterion.

Default:

```python
confidence=0.95
```

#### `target_p0`

Target empirical background exceedance rate used during automatic threshold selection. It is not itself a site-level false-positive rate.

Default:

```python
target_p0=0.05
```

`target_p0` can currently be set through the Python API. Command-line runs use the default value.

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
    target_p0=0.05,
)
```

---

## Outputs

By default, betaSieve writes:

- `results/csv/min_max_difference_{threshold}.csv`: site-level statistics, per-sample ranges, raw and adjusted p-values, and flag columns;
- `results/csv/candidate_cpgs.csv`: all IlmnIDs belonging to sites with `p_empir_adj_flagged=True`;
- `results/csv/threshold_sweep_summary.csv`: threshold-specific background and flagging rates when a sweep was performed;
- `results/report.html`: interactive analysis and null-model diagnostics.

The principal statistical columns are:

- `above_threshold`: number of samples whose max–min range exceeds \(t\);
- `p_hat`: site-level sample exceedance rate;
- `p0`: pooled exact-replicate background exceedance rate;
- `p_empir`: empirical upper-tail p-value;
- `p_empir_adj`: multiple-testing-adjusted empirical p-value;
- `p_empir_adj_flagged`: primary discordance flag.

`candidate_cpgs.csv` is a recommendation list; betaSieve does not modify the input β-value matrix automatically.

## Statistical considerations

- Empirical p-values are discrete, with minimum possible value \(1/(m+1)\).
- The empirical test assumes that exact-replicate sites provide an appropriate reference distribution for non-replicate design comparisons.
- Max–min ranges can increase with the number of probes in a group, so comparisons involving pairs, triplets, and quadruplets should be interpreted with that difference in mind.
- The current implementation assumes a complete sample matrix; missing β-values require careful preprocessing.
