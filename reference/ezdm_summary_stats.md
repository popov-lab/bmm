# Compute Robust Summary Statistics for EZ-Diffusion Model

Computes robust summary statistics for the EZ-Diffusion Model by fitting
mixture models to raw trial-level RT data, separating contaminant
responses from true responses.

## Usage

``` r
ezdm_summary_stats(
  rt,
  response,
  version = c("3par", "4par"),
  distribution = c("exgaussian", "lognormal", "invgaussian"),
  method = c("mixture", "simple", "robust"),
  robust_scale = c("iqr", "mad"),
  contaminant_bound = c("min", "max"),
  min_trials = 10,
  init_contaminant = 0.05,
  max_contaminant = 0.5,
  maxit = 100,
  tol = 1e-06
)
```

## Arguments

- rt:

  Numeric vector of reaction times in seconds.

- response:

  Vector of response indicators. Accepts multiple formats:

  - Numeric: 1 = upper/correct, 0 = lower/error

  - Logical: TRUE = upper/correct, FALSE = lower/error

  - Character/Factor: "upper"/"lower", "correct"/"error", "acc"/"err",
    "hit"/"miss", "yes"/"no" (case-insensitive)

- version:

  Character. Either "3par" (default) for pooled RTs or "4par" for
  separate upper/lower boundary RTs. Controls the output columns.

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
  "min" and "max" to use data-driven bounds (default):

  - "min": Use the minimum RT in each group, minus a 50\\

  - "max": Use the maximum RT in each group, plus a 50\\

  - Numeric: Fixed bounds, e.g., c(0.1, 3.0)

  The buffer extends data-driven bounds to ensure conservative
  estimates. Examples: c(0.1, 3.0), c("min", "max"), c(0.1, "max"),
  c("min", 3.0)

- min_trials:

  Integer. Minimum number of trials required for fitting. Returns NA if
  fewer trials are available. Default is 10

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

## Value

A 1-row `data.frame`. For version = "3par": `mean_rt`, `var_rt`,
`n_upper`, `n_trials`, `contaminant_prop`. For version = "4par":
`mean_rt_upper`, `mean_rt_lower`, `var_rt_upper`, `var_rt_lower`,
`n_upper`, `n_trials`, `contaminant_prop_upper`,
`contaminant_prop_lower`.

## Details

RT outliers and contaminant responses (fast guesses, lapses of
attention) can distort the mean and variance estimates used as input to
the EZ-Diffusion equations. This function addresses this by fitting a
mixture model with two components: a uniform distribution for
contaminants and a parametric RT distribution for true responses. Robust
moments are then extracted from the fitted parametric component.

This function is designed to work with
[`dplyr::group_by()`](https://dplyr.tidyverse.org/reference/group_by.html)
and
[`dplyr::reframe()`](https://dplyr.tidyverse.org/reference/reframe.html)
for grouped operations. Use
[`adjust_ezdm_accuracy()`](https://venpopov.com/bmm/reference/adjust_ezdm_accuracy.md)
as a separate step if you need to adjust accuracy counts for
contamination.

## See also

[`adjust_ezdm_accuracy()`](https://venpopov.com/bmm/reference/adjust_ezdm_accuracy.md)
for adjusting accuracy counts,
[`flag_contaminant_rts()`](https://venpopov.com/bmm/reference/flag_contaminant_rts.md)
for trial-level contamination probabilities,
[`ezdm()`](https://venpopov.com/bmm/reference/ezdm.md) for fitting the
EZ-Diffusion Model

## Examples

``` r
# Generate example data
set.seed(123)
rt <- rgamma(100, shape = 5, rate = 10) + 0.3
response <- rbinom(100, 1, 0.8)

# 3par summary stats
ezdm_summary_stats(rt, response)
#>      mean_rt     var_rt n_upper n_trials contaminant_prop
#> mu 0.7751987 0.03617632      82      100      4.68105e-08

# With dplyr for grouped operations
# library(dplyr)
# mydata |>
#   group_by(subject) |>
#   reframe(ezdm_summary_stats(rt, response))

# 4par version with separate upper/lower moments
ezdm_summary_stats(rt, response, version = "4par")
#>    mean_rt_upper mean_rt_lower var_rt_upper var_rt_lower n_upper n_trials
#> mu     0.7751185     0.7756404   0.04176804   0.01614607      82      100
#>    contaminant_prop_upper contaminant_prop_lower
#> mu           5.349368e-08           5.258552e-08
```
