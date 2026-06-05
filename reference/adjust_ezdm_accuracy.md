# Adjust Accuracy Counts for Contamination

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
  [`ezdm_summary_stats()`](https://venpopov.com/bmm/reference/ezdm_summary_stats.md)).

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

[`ezdm_summary_stats()`](https://venpopov.com/bmm/reference/ezdm_summary_stats.md)
for computing the summary statistics and contamination proportions

## Examples

``` r
# Adjust accuracy for estimated 10% contamination
set.seed(42)
adjust_ezdm_accuracy(n_upper = 80, n_trials = 100, contaminant_prop = 0.1)
#>   n_upper_adj n_trials_adj
#> 1          70           86

# In a pipeline with ezdm_summary_stats
# library(dplyr)
# mydata |>
#   group_by(subject) |>
#   reframe(ezdm_summary_stats(rt, response)) |>
#   mutate(adjust_ezdm_accuracy(n_upper, n_trials, contaminant_prop))
```
