# Distribution functions for m-AFC SDT

Density and random generation for m-alternative forced choice signal
detection theory (DeCarlo, 2012). Models accuracy in tasks where one of
`m` alternatives contains the signal. Only the `d` parameter is
estimated (no criterion). All arguments are recycled to the length of
the longest one, so passing vectors of `d`, `m`, or `n_trials` generates
(or evaluates) one observation per element.

## Usage

``` r
dsdt_mafc(
  n_correct,
  n_trials,
  m,
  d,
  dist = c("normal", "gumbel_min", "gumbel_max", "logistic"),
  log = FALSE
)

rsdt_mafc(
  n,
  n_trials,
  m,
  d,
  dist = c("normal", "gumbel_min", "gumbel_max", "logistic")
)
```

## Arguments

- n_correct:

  Integer vector. Number of correct responses.

- n_trials:

  Integer vector. Total number of trials per observation.

- m:

  Integer vector. Number of alternatives per observation. Must be at
  least 2.

- d:

  Numeric vector. Sensitivity \\d'\\: the distance between the signal
  and distractor distributions in SD units. m-AFC assumes a common scale
  for the two distributions, so this is also the balanced index \\d_a\\
  that [`sdt_yn()`](https://venpopov.com/bmm/dev/reference/sdt_yn.md)
  reports.

- dist:

  Character. The distribution assumed for the latent evidence, given
  here by its cumulative distribution function:

  - "normal" (default): Gaussian, \\\Phi(x)\\

  - "gumbel_min": smallest extreme value, \\1 - \exp(-\exp(x))\\ (the
    complementary log-log distribution)

  - "gumbel_max": largest extreme value, \\\exp(-\exp(-x))\\ (the
    log-log distribution, as in
    [`evd::pgumbel`](https://rdrr.io/pkg/evd/man/gumbel.html))

  - "logistic": \\1 / (1 + \exp(-x))\\

- log:

  Logical. If `TRUE`, returns log-density (default `FALSE`).

- n:

  Integer. Number of observations to generate. `n_trials`, `m`, and `d`
  are recycled to this length.

## Value

`dsdt_mafc` returns the (log-)density (binomial probability).
`rsdt_mafc` returns an integer vector with the number of correct
responses per observation.

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

DeCarlo, L. T. (2012). On a signal detection approach to m-alternative
forced choice with bias, with maximum likelihood and Bayesian approaches
to estimation. *Journal of Mathematical Psychology*, *56*(3), 196–207.
[doi:10.1016/j.jmp.2012.02.004](https://doi.org/10.1016/j.jmp.2012.02.004)

## Examples

``` r
# 4-AFC density
dsdt_mafc(n_correct = 80, n_trials = 100, m = 4, d = 1.5)
#> [1] 0.008272764
# Generate 4-AFC data for 20 subjects with varying sensitivity
dat <- data.frame(id = 1:20, n_trials = 200L)
dat$n_correct <- rsdt_mafc(nrow(dat), dat$n_trials, m = 4,
                           d = rnorm(20, 1.5, 0.4))
head(dat)
#>   id n_trials n_correct
#> 1  1      200       158
#> 2  2      200       119
#> 3  3      200       145
#> 4  4      200       139
#> 5  5      200       131
#> 6  6      200       125
```
