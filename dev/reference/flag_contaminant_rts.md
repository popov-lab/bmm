# Flag contaminant reaction times using mixture modeling

Identifies contaminant RTs (fast guesses, attention lapses) at the trial
level using mixture modeling. For each trial, it computes the posterior
probability of being a contaminant given a mixture of a uniform
distribution (contaminants) and an RT distribution.

The function takes a numeric vector of RTs and returns a numeric vector
of contamination probabilities, making it compatible with
[`dplyr::mutate()`](https://dplyr.tidyverse.org/reference/mutate.html)
and
[`dplyr::group_by()`](https://dplyr.tidyverse.org/reference/group_by.html)
workflows.

## Usage

``` r
flag_contaminant_rts(
  rt,
  distribution = c("exgaussian", "lognormal", "invgaussian"),
  contaminant_bound = c("min", "max"),
  init_contaminant = 0.05,
  max_contaminant = 0.5,
  maxit = 100,
  tol = 1e-06
)
```

## Arguments

- rt:

  Numeric vector. Reaction times in seconds. Must be positive.

- distribution:

  Character. RT distribution for the mixture model: "exgaussian"
  (default), "lognormal", or "invgaussian".

- contaminant_bound:

  Vector of length 2. Bounds `[lower, upper]` for the uniform
  contaminant distribution. Can be numeric values or "min"/"max" for
  data-driven bounds. Default `c("min", "max")`.

- init_contaminant:

  Numeric. Initial contaminant proportion for EM algorithm. Must be in
  (0, 1). Default 0.05.

- max_contaminant:

  Numeric. Maximum allowed contaminant proportion. Values exceeding this
  are clipped with a warning. Must be in (0, 1\]. Default 0.5.

- maxit:

  Integer. Maximum EM iterations. Default 100.

- tol:

  Numeric. Convergence tolerance for log-likelihood. Default 1e-6.

## Value

Numeric vector of posterior contamination probabilities P(contaminant \|
RT), with a `"diagnostics"` attribute containing a one-row data.frame
with columns: `mixture_params` (list), `contaminant_prop`, `converged`,
`iterations`, `loglik`, `n_trials`, `distribution`, `method`.

## Details

### Mixture Model

The function fits:
`f(RT) = pi_c * Uniform(a,b) + (1-pi_c) * f_RT(RT|theta)`

where pi_c is the contaminant proportion, Uniform(a,b) is the
contaminant distribution over `contaminant_bound`, and f_RT is the
specified RT distribution with parameters theta.

### Grouping

To fit separate mixtures by condition or response boundary, use
[`dplyr::group_by()`](https://dplyr.tidyverse.org/reference/group_by.html)
before calling this function inside
[`dplyr::mutate()`](https://dplyr.tidyverse.org/reference/mutate.html).

### Diagnostics

Mixture fit diagnostics (parameters, convergence, log-likelihood) are
attached as the `"diagnostics"` attribute of the returned vector. Access
them with `attr(result, "diagnostics")`.

## See also

[`ezdm_summary_stats()`](https://venpopov.com/bmm/dev/reference/ezdm_summary_stats.md)
for aggregated RT statistics with contamination handling,
[`validate_fast_guesses()`](https://venpopov.com/bmm/dev/reference/validate_fast_guesses.md)
for testing whether flagged contaminants show random guessing behavior

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate data with contaminants
library(bmm)
set.seed(123)
rt_clean <- rgamma(150, shape = 5, rate = 10)
rt_contam <- runif(50, 0.1, 0.2)

data <- data.frame(
  rt = c(rt_clean, rt_contam),
  subject = 1,
  response = sample(c("upper", "lower"), 200, replace = TRUE)
)

# Basic usage with mutate
library(dplyr)
data <- data |>
  mutate(contam_prob = flag_contaminant_rts(rt))

# Hard threshold: remove trials with P(contaminant) > 0.5
data_clean <- data |> filter(contam_prob <= 0.5)

# Separate fits by response boundary
data <- data |>
  group_by(subject, response) |>
  mutate(contam_prob = flag_contaminant_rts(rt))

# Access diagnostics
probs <- flag_contaminant_rts(data$rt)
attr(probs, "diagnostics")
} # }
```
