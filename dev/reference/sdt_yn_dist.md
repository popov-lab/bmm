# Distribution functions for Yes/No SDT

Density and random generation for the yes/no signal detection theory
model, where the response is the number of "old"/"signal" responses out
of a fixed number of trials (a binomial likelihood).

## Usage

``` r
dsdt_yn(
  n_old,
  n_trials,
  stimulus,
  d,
  criterion = NULL,
  bias = NULL,
  sdratio = 1,
  dist = c("normal", "gumbel_min", "gumbel_max", "logistic"),
  log = FALSE
)

rsdt_yn(
  n,
  n_trials,
  stimulus,
  d,
  criterion = NULL,
  bias = NULL,
  sdratio = 1,
  dist = c("normal", "gumbel_min", "gumbel_max", "logistic")
)
```

## Arguments

- n_old:

  Integer vector. Number of "old"/"signal" responses.

- n_trials:

  Integer vector. Total number of trials per cell.

- stimulus:

  Numeric or logical vector (0/1). Stimulus type: 0 = noise, 1 = signal.

- d:

  Numeric. Sensitivity: \\d'\\ when `sdratio` is 1, and otherwise the
  balanced index \\d_a\\, the separation between the two distributions
  divided by the root-mean-square of their SDs (see
  [`sdt_yn()`](https://venpopov.com/bmm/dev/reference/sdt_yn.md)). The
  separation in noise units is `d * sqrt((1 + sdratio^2) / 2)`.

- criterion:

  Numeric. Location of the decision criterion, measured from the mean of
  the noise distribution in units of its SD. Give either `criterion` or
  `bias`, not both.

- bias:

  Numeric. Location of the decision criterion measured from the midpoint
  between the noise and signal means instead, in noise-SD units:
  `criterion - d * sqrt((1 + sdratio^2) / 2) / 2`. 0 is the neutral
  point, positive values a bias toward "new"/"noise".

- sdratio:

  Numeric. Ratio of signal to noise standard deviations (default 1,
  i.e., equal variance). Must be positive, and is on the **natural**
  scale — see the section below before reusing a fitted value.

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

  Integer. Number of observations to generate. `n_trials`, `stimulus`,
  and the model parameters are recycled to this length.

## Value

`dsdt_yn` returns the (log-)density (binomial probability). `rsdt_yn`
returns an integer vector with the number of "old"/"signal" responses
per observation.

## Parameter scales

As everywhere in `bmm`, these functions take their arguments on the
**natural** scale, while the parameters
[`sdt_yn()`](https://venpopov.com/bmm/dev/reference/sdt_yn.md)
*estimates* are on their link scale. `d` and `criterion` have identity
links and carry across unchanged, but `sdratio` has a log link: a fitted
`sdratio` of 0.375 is a ratio of `exp(0.375) = 1.455`, and passing
`0.375` here instead asks for a signal distribution 2.7 times narrower
than the noise. That is a legal value and raises no error, so
exponentiate first.

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

Green, D. M., & Swets, J. A. (1966). *Signal detection theory and
psychophysics*. Wiley.

## Examples

``` r
# Density of yes/no SDT data
dsdt_yn(n_old = 80, n_trials = 100, stimulus = 1,
        d = 1.5, criterion = 0.9)
#> [1] 0.02264629

# The same observer, with the criterion given relative to the midpoint
dsdt_yn(n_old = 80, n_trials = 100, stimulus = 1,
        d = 1.5, bias = 0.15)
#> [1] 0.02264629

# Vectorized over observations
dsdt_yn(n_old = c(30, 80), n_trials = c(100, 100),
        stimulus = c(0, 1), d = 1.5, criterion = 0.9,
        log = TRUE)
#> [1] -6.371711 -3.787759

# Unequal variance from a fitted model: sdt_yn() reports sdratio on its log
# link, so exponentiate before passing it here
dsdt_yn(n_old = 80, n_trials = 100, stimulus = 1,
        d = 1.5, criterion = 0.9, sdratio = exp(0.375))
#> [1] 0.04679876
# Generate yes/no SDT data for a design
dat <- expand.grid(id = 1:20, stimulus = c(0L, 1L))
dat$n_trials <- 100L
dat$n_old <- rsdt_yn(nrow(dat), dat$n_trials, dat$stimulus,
                     d = 1.5, criterion = 0.9)
head(dat)
#>   id stimulus n_trials n_old
#> 1  1        0      100    22
#> 2  2        0      100    16
#> 3  3        0      100    19
#> 4  4        0      100    16
#> 5  5        0      100    21
#> 6  6        0      100    18
```
