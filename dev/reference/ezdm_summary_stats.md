# Compute Robust Summary Statistics for EZ-Diffusion Model

Computes robust summary statistics for the EZ-Diffusion Model by fitting
mixture models to raw trial-level RT data, separating contaminant
responses from true responses.

## Usage

``` r
ezdm_summary_stats(
  data,
  rt,
  response,
  .by = NULL,
  version = c("3par", "4par"),
  distribution = c("exgaussian", "lognormal", "invgaussian"),
  method = c("mixture", "simple", "robust"),
  robust_scale = c("iqr", "mad"),
  contaminant_bound = c(0.1, 3),
  min_trials = 10,
  init_contaminant = 0.05,
  max_contaminant = 0.5,
  maxit = 100,
  tol = 1e-06,
  adjust_accuracy = FALSE,
  guess_rate = 0.5
)
```

## Arguments

- data:

  A `data.frame` containing trial-level data with RT and accuracy
  columns

- rt:

  Character. The name of the column containing reaction times (in
  seconds)

- response:

  Character. The name of the column containing response indicators.
  Accepts multiple formats:

  - Numeric: 1 = upper/correct, 0 = lower/error

  - Logical: TRUE = upper/correct, FALSE = lower/error

  - Character/Factor: "upper"/"lower", "correct"/"error", "acc"/"err",
    "hit"/"miss", "yes"/"no" (case-insensitive)

- .by:

  A character vector of column names to group by before computing
  summary statistics (e.g., `.by = c("subject", "condition")`). If NULL
  (default), computes statistics across all data without grouping.

- version:

  Character. Either "3par" (default) for pooled RTs or "4par" for
  separate upper/lower boundary RTs

- distribution:

  Character. The parametric distribution for the RT component. One of
  "exgaussian" (default), "lognormal", or "invgaussian"

- method:

  Character. One of "mixture" (default) for robust estimation via
  mixture modeling, "robust" for non-parametric robust estimation using
  median and IQR/MAD-based variance, or "simple" for standard moment
  calculation. The "robust" method is faster and requires no
  distributional assumptions, but note that the EZ equations were
  derived for mean and variance, so using median may introduce some bias
  for skewed distributions.

- robust_scale:

  Character. Scale estimator for robust method. Either "iqr" (default)
  for IQR-based variance estimation (variance = (IQR/1.349)^2) or "mad"
  for MAD-based estimation (variance = MAD^2, where MAD is scaled to be
  consistent with SD for normal data). Only used when method = "robust".

- contaminant_bound:

  Vector of length 2 specifying the bounds (in seconds) for the uniform
  contaminant distribution. Can be numeric values or the special strings
  "min" and "max" to use data-driven bounds:

  - Numeric: Fixed bounds, e.g., c(0.1, 3.0) (default)

  - "min": Use the minimum RT in each group, minus a 50\\

  - "max": Use the maximum RT in each group, plus a 50\\

  The buffer extends data-driven bounds to ensure conservative
  estimates. Examples: c(0.1, 3.0), c("min", "max"), c(0.1, "max"),
  c("min", 3.0)

- min_trials:

  Integer. Minimum number of trials required for fitting. Groups with
  fewer trials will return NA. Default is 10

- init_contaminant:

  Numeric. Initial proportion of contaminants for EM algorithm. Default
  is 0.05

- max_contaminant:

  Numeric. Maximum allowed contaminant proportion (0 \< max \<= 1).
  Estimates are clipped to this value to prevent inflated contaminant
  proportions. Default is 0.5

- maxit:

  Integer. Maximum number of EM iterations. Default is 100

- tol:

  Numeric. Convergence tolerance for EM algorithm. Default is 1e-6

- adjust_accuracy:

  Logical. If TRUE and method = "mixture", adjust accuracy counts by
  removing estimated contaminant guesses using binomial sampling.
  Default is FALSE

- guess_rate:

  Numeric. Assumed accuracy rate for contaminant trials (random
  guessing). Default is 0.5 (appropriate for 2AFC tasks)

## Value

A `data.frame` with summary statistics. For version = "3par": grouping
variables, mean_rt, var_rt, n_upper, n_trials, contaminant_prop. When
adjust_accuracy = TRUE, also includes n_upper_adj and n_trials_adj. For
version = "4par": grouping variables, mean_rt_upper, mean_rt_lower,
var_rt_upper, var_rt_lower, n_upper, n_trials, contaminant_prop_upper,
contaminant_prop_lower.

## Details

RT outliers and contaminant responses (fast guesses, lapses of
attention) can distort the mean and variance estimates used as input to
the EZ-Diffusion equations. This function addresses this by fitting a
mixture model with two components: a uniform distribution for
contaminants and a parametric RT distribution for true responses. Robust
moments are then extracted from the fitted parametric component.

## Examples

``` r
# Generate example data
set.seed(123)
test_data <- data.frame(
  subject = rep(1:3, each = 100),
  condition = rep(c("A", "B"), 150),
  rt = rgamma(300, shape = 5, rate = 10) + 0.3,
  correct = rbinom(300, 1, 0.8)
)

# Compute summary statistics grouped by subject
result <- ezdm_summary_stats(test_data,
  rt = "rt", response = "correct",
  .by = "subject"
)
print(result)
#>   subject   mean_rt     var_rt n_upper n_trials contaminant_prop
#> 1       1 0.7752000 0.03617616      85      100     5.232355e-09
#> 2       2 0.7722113 0.06167891      78      100     1.008148e-08
#> 3       3 0.7906609 0.04786826      85      100     1.694455e-08

# Group by multiple variables using simple method
result_multi <- ezdm_summary_stats(test_data,
  rt = "rt",
  response = "correct",
  .by = c("subject", "condition"),
  method = "simple"
)
```
