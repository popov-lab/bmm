# Adjust Accuracy Counts for Contamination (deprecated)

**Deprecated.**
[`ezdm_summary_stats()`](https://venpopov.com/bmm/dev/reference/ezdm_summary_stats.md)
now returns contaminant-free `n_upper` and `n_trials`, so this second
step is no longer needed; applying it corrects the same cell twice. Set
the guess rate with the `guess_rate` argument of
[`ezdm_summary_stats()`](https://venpopov.com/bmm/dev/reference/ezdm_summary_stats.md)
instead.

Adjusts accuracy counts (`n_upper`, `n_trials`) by removing estimated
contaminant trials using binomial sampling. Contaminant trials are
assumed to produce correct responses at a fixed guess rate (e.g., 0.5
for 2AFC tasks).

## Usage

``` r
adjust_ezdm_accuracy(n_upper, n_trials, contaminant_prop, guess_rate = 0.5)
```

## Arguments

- n_upper:

  Numeric. Count of upper boundary (correct) responses.

- n_trials:

  Numeric. Total number of trials.

- contaminant_prop:

  Numeric. Estimated proportion of contaminant trials (e.g., from the
  `contaminant_prop` column of
  [`ezdm_summary_stats()`](https://venpopov.com/bmm/dev/reference/ezdm_summary_stats.md)).

- guess_rate:

  Numeric. Assumed accuracy rate for contaminant trials (random
  guessing). Default is 0.5 (appropriate for 2AFC tasks).

## Value

A 1-row `data.frame` with columns `n_upper_adj` and `n_trials_adj`
(integers). When `contaminant_prop` is `NA` or \<= 0, returns the
original counts unchanged.

## Details

Uses binomial sampling to estimate the number of contaminant trials and
contaminant correct responses, then subtracts these from the raw counts.
Because of the stochastic sampling, results will vary across calls
unless a seed is set by the user.

## See also

[`ezdm_summary_stats()`](https://venpopov.com/bmm/dev/reference/ezdm_summary_stats.md)
for computing the summary statistics and contamination proportions

## Examples

``` r
# Deprecated. ezdm_summary_stats() returns contaminant-free counts already:
set.seed(42)
rt <- c(rnorm(80, 0.55, 0.05), runif(20, 0.1, 4))
response <- c(rbinom(80, 1, 0.85), rbinom(20, 1, 0.5))
ezdm_summary_stats(rt, response, contaminant_bound = c(0.1, 4))
#>      mean_rt      var_rt n_upper n_trials contaminant_prop
#> mu 0.5532669 0.002680182      67       80        0.1995584
```
