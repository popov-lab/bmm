# Distribution functions for the psychometric function model

Density and random generation for the
[`psychometric_function()`](https://venpopov.com/bmm/dev/reference/psychometric_function.md)
model: the number of positive responses ("yes", "correct", "longer") out
of `n_trials` at a stimulus intensity, a binomial with probability
\$\$\psi(x) = \gamma + (1 - \gamma)(1 - \lambda) F(x),\$\$ where \\F\\
is the sigmoid with its midpoint and width.

## Usage

``` r
dpsychometric_function(
  x,
  intensity,
  midpoint,
  width,
  guess = 0,
  lapse = 0,
  n_trials = 1,
  sigmoid = c("normal", "logistic", "gumbel_min", "gumbel_max", "weibull", "lognormal",
    "loglogistic"),
  log = FALSE
)

rpsychometric_function(
  n,
  intensity,
  midpoint,
  width,
  guess = 0,
  lapse = 0,
  n_trials = 1,
  sigmoid = c("normal", "logistic", "gumbel_min", "gumbel_max", "weibull", "lognormal",
    "loglogistic")
)
```

## Arguments

- x:

  Integer vector. Number of positive responses, between 0 and
  `n_trials`.

- intensity:

  Numeric vector. Stimulus intensity. Must be positive for the sigmoids
  that read it on the log scale (`"weibull"`, `"lognormal"`,
  `"loglogistic"`).

- midpoint:

  Numeric. Intensity at which the sigmoid \\F\\ is 0.5, in the units of
  `intensity`. Positive for the log-scale sigmoids.

- width:

  Numeric, positive. Distance on the intensity axis between the points
  where \\F\\ is 0.05 and 0.95; for the log-scale sigmoids, in log units
  of intensity.

- guess:

  Numeric in \[0, 1). Guess rate, the lower asymptote: 0.5 in 2AFC,
  `1 / m` in m-AFC.

- lapse:

  Numeric in \[0, 1). Proportion of trials on which the observer does
  not attend and responds at the guess rate. The upper asymptote is
  `1 - lapse * (1 - guess)`.

- n_trials:

  Integer vector. Number of trials per observation (default 1, a single
  binary response).

- sigmoid:

  The sigmoid \\F\\; see
  [`psychometric_function()`](https://venpopov.com/bmm/dev/reference/psychometric_function.md).

- log:

  Logical. If `TRUE`, returns the log-density (default `FALSE`).

- n:

  Integer. Number of observations to generate. `intensity`, `n_trials`
  and the parameters are recycled to this length. As in
  [`stats::rnorm()`](https://rdrr.io/r/stats/Normal.html), a fraction is
  truncated and `n = 0` gives an empty vector; unlike
  [`stats::rnorm()`](https://rdrr.io/r/stats/Normal.html), a vector `n`
  stops.

## Value

`dpsychometric_function()` returns the (log-)density, a binomial
probability. `rpsychometric_function()` returns an integer vector with
the number of positive responses per observation.

## Parameter scales

These functions take every argument on the **natural** scale, while
[`psychometric_function()`](https://venpopov.com/bmm/dev/reference/psychometric_function.md)
*estimates* `width` on the log scale and `guess` and `lapse` on the
logit scale (and `midpoint` on the log scale for the log-scale
sigmoids). Transform a fitted value back before passing it here.

## Vector arguments

- `r*(n, ...)` returns exactly `n` draws: a vector, or a data frame or
  matrix with `n` rows. Every per-observation argument (a parameter, a
  stimulus, a number of trials) is recycled to length `n` with
  [`rep_len()`](https://rdrr.io/r/base/rep.html), as in
  [`stats::rnorm()`](https://rdrr.io/r/stats/Normal.html): draw `i` uses
  element `i`, and a matrix is read element by element; only a set of
  values per draw, such as `thresholds`, takes row `i` of a matrix (see
  below). A shorter vector is repeated, also when `n` is not a multiple
  of its length. Until bmm 1.6.0 at the earliest,
  [`rddm()`](https://venpopov.com/bmm/dev/reference/ddm_dist.md) with
  `n = 1` and parameter vectors still returns one draw per element, with
  a deprecation warning.

- `n` is a single number of at least 0. A fraction is truncated, as in
  [`stats::rnorm()`](https://rdrr.io/r/stats/Normal.html), and `n = 0`
  gives an empty vector, or a matrix or data frame with zero rows and
  the usual columns. Unlike
  [`stats::rnorm()`](https://rdrr.io/r/stats/Normal.html), which takes
  the length of a vector `n`, these functions stop on one;
  [`rddm()`](https://venpopov.com/bmm/dev/reference/ddm_dist.md) still
  takes its length, with a deprecation warning.

- `d*(x, ...)` returns one value per observation: as many as the longest
  of `x` (its rows, when the observations are rows of counts) and the
  per-observation arguments, the others recycled to that length. Value
  `i` is what the call on the `i`-th elements alone returns. A matrix in
  an argument that takes one number per observation, such as `x` or `d`,
  holds one observation per element, as in
  [`stats::dnorm()`](https://rdrr.io/r/stats/Normal.html). When `x` or
  any per-observation argument is empty, the result is empty, as in
  [`stats::dnorm()`](https://rdrr.io/r/stats/Normal.html). A `NULL`
  argument, such as a misspelt column of a data frame, counts as empty,
  except where the default is `NULL`, as for `criterion` or `kcrit`:
  there it means that the argument is not given. Where an observation is
  a vector of counts
  ([`dsdt_rating()`](https://venpopov.com/bmm/dev/reference/sdt_rating_dist.md),
  [`dsdt_dpsdt()`](https://venpopov.com/bmm/dev/reference/sdt_dpsdt_dist.md),
  [`dsdt_metad()`](https://venpopov.com/bmm/dev/reference/sdt_metad_dist.md),
  [`dsdt_ranking()`](https://venpopov.com/bmm/dev/reference/sdt_ranking_dist.md)
  and
  [`dsdt_cdp()`](https://venpopov.com/bmm/dev/reference/sdt_cdp_dist.md)),
  or the pair of `mean_rt` or `var_rt` values of
  `dezdm(version = "4par")`, a vector is always one observation: pass a
  matrix with zero rows for none, as an empty vector stops. An empty
  `kcrit` of
  [`dsdt_cdp()`](https://venpopov.com/bmm/dev/reference/sdt_cdp_dist.md)
  and
  [`rsdt_cdp()`](https://venpopov.com/bmm/dev/reference/sdt_cdp_dist.md)
  stops as well; `kcrit = NULL` gives the R/K model. Counts in a data
  frame, a tibble included, are read as the matrix of its rows, and so
  are `mean_rt` and `var_rt` of `dezdm(version = "4par")`.
  [`dsdm()`](https://venpopov.com/bmm/dev/reference/SDMdist.md) with a
  data frame `x` is the exception: it returns a data frame of the same
  shape, one density per cell, and cuts longer parameters to the number
  of cells.

- A set of values that describes one observation, such as the
  `thresholds` of the rating models, is a vector when every observation
  shares it and a matrix with one row per observation otherwise: row `i`
  belongs to observation `i`, or to draw `i`.

- Configuration arguments such as `version`, `dist`, `log`,
  `parametrization` and `n_new` are single values and are never
  recycled. An argument that a version does not use, such as `zr` of
  `dezdm(version = "3par")`, is ignored.

## References

Schütt, H. H., Harmeling, S., Macke, J. H., & Wichmann, F. A. (2016).
Painfree and accurate Bayesian estimation of psychometric functions for
(potentially) overdispersed data. *Vision Research*, *122*, 105–123.
[doi:10.1016/j.visres.2016.02.002](https://doi.org/10.1016/j.visres.2016.02.002)

Wichmann, F. A., & Hill, N. J. (2001). The psychometric function: I.
Fitting, sampling, and goodness of fit. *Perception & Psychophysics*,
*63*(8), 1293–1313.
[doi:10.3758/BF03194544](https://doi.org/10.3758/BF03194544)

## Examples

``` r
# Probability of 7 correct out of 10 2AFC trials at the midpoint
dpsychometric_function(7, intensity = 1, midpoint = 1, width = 2, guess = 0.5,
                       n_trials = 10)
#> [1] 0.2502823

# Simulate a method-of-constant-stimuli experiment with a Weibull observer
dat <- data.frame(contrast = rep(c(0.01, 0.02, 0.04, 0.08, 0.16), 4))
dat$n_correct <- rpsychometric_function(nrow(dat), dat$contrast, midpoint = 0.04,
                                        width = 2, guess = 0.5, lapse = 0.02,
                                        n_trials = 40, sigmoid = "weibull")
head(dat)
#>   contrast n_correct
#> 1     0.01        19
#> 2     0.02        21
#> 3     0.04        24
#> 4     0.08        40
#> 5     0.16        40
#> 6     0.01        23
```
