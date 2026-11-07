# Thresholds and slopes of a fitted psychometric function

Computes, per posterior draw and condition, the stimulus intensity at
which a
[`psychometric_function()`](https://venpopov.com/bmm/dev/reference/psychometric_function.md)
fit reaches a target response probability, the slope of the psychometric
function there, and the guess and lapse rates, including the lapse rate
as Wichmann and Hill (2001) define it.

## Usage

``` r
psychometric_threshold(
  fit,
  p = 0.75,
  conditions = NULL,
  probs = c(0.025, 0.975),
  ...
)
```

## Arguments

- fit:

  A `bmmfit` object returned by
  [`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md) with a
  [`psychometric_function()`](https://venpopov.com/bmm/dev/reference/psychometric_function.md)
  model.

- p:

  Numeric vector of target response probabilities, each strictly between
  0 and 1. The default 0.75 is the conventional 2AFC threshold.

- conditions:

  A data frame of predictor values at which to evaluate the parameters.
  The default uses every unique combination of the population-level
  predictors in the data, and of the guess or lapse column when a rate
  is fixed per row. A user-supplied data frame must have a column for
  each of these and no other: the thresholds are population-level, so a
  grouping column such as a participant id is refused.

- probs:

  Numeric vector of length 2. Lower and upper quantiles for the credible
  interval in the summary. Default `c(0.025, 0.975)`.

- ...:

  Further arguments passed to
  [`brms::posterior_linpred()`](https://mc-stan.org/rstantools/reference/posterior_linpred.html),
  such as `draw_ids`.

## Value

A data frame of class `"bmm_psychometric_threshold"` with one row per
draw, measure and condition, and columns `measure` (`"threshold"`,
`"slope"`, `"guess"`, `"lapse"`, `"lapse_wh"`), `p` (`NA` for the three
rates, which do not depend on it), `value`, `.draw` and the condition
columns. `lapse_wh` is `lapse * (1 - guess)`, the gap between the upper
asymptote and 1, which psignifit and Palamedes report as the lapse rate.
Attribute `summary` holds the posterior mean, the credible interval and
the share of `undefined` draws per measure, `p` and condition.

## Details

The threshold at target probability \\p\\ solves \\\gamma + (1 -
\gamma)(1 - \lambda) F(x) = p\\. It exists only between the lower
asymptote, `guess`, and the upper asymptote, `1 - lapse * (1 - guess)`.
A draw whose asymptotes do not enclose `p` has no threshold: its
threshold and slope are `NA`, the summary reports the share of such
draws as `undefined`, and the function warns. Asking for 0.75 in a 2AFC
task with `guess = 0.5` is the usual case; asking for 0.5 is not,
because no 2AFC observer falls below chance.

The slope is \\d\psi / dx\\ at the threshold, in probability per unit of
intensity.

Thresholds, slopes and rates come from the population-level effects only
(`re_formula = NA` in
[`brms::posterior_linpred()`](https://mc-stan.org/rstantools/reference/posterior_linpred.html)):
group-level effects are left out, and `re_formula` cannot be passed
through `...`. The closed form needs `guess` and `lapse` constant along
the intensity axis, so a fit whose formula predicts either by the
intensity is refused.

## References

Wichmann, F. A., & Hill, N. J. (2001). The psychometric function: I.
Fitting, sampling, and goodness of fit. *Perception & Psychophysics*,
*63*(8), 1293–1313.
[doi:10.3758/BF03194544](https://doi.org/10.3758/BF03194544)

## See also

[`psychometric_function()`](https://venpopov.com/bmm/dev/reference/psychometric_function.md)

## Examples

``` r
if (FALSE) { # \dontrun{
dat <- expand.grid(contrast = c(0.01, 0.02, 0.04, 0.08, 0.16),
                   cond = c("a", "b"))
dat$n_trials <- 60L
dat$n_correct <- rpsychometric_function(nrow(dat), dat$contrast,
                                        midpoint = ifelse(dat$cond == "a", 0.03, 0.05),
                                        width = 2, guess = 0.5, lapse = 0.02,
                                        n_trials = 60, sigmoid = "weibull")
fit <- bmm(
  bmf(midpoint ~ 0 + cond, width ~ 1, lapse ~ 1), dat,
  psychometric_function("n_correct", "contrast", "n_trials", sigmoid = "weibull",
                        guess = 0.5),
  cores = 4, backend = "cmdstanr"
)
thr <- psychometric_threshold(fit, p = c(0.75, 0.9))
thr

# the threshold difference between the conditions, per draw
wide <- subset(thr, measure == "threshold" & p == 0.75)
diff_ab <- wide$value[wide$cond == "b"] - wide$value[wide$cond == "a"]
quantile(diff_ab, c(0.025, 0.5, 0.975))
} # }
```
