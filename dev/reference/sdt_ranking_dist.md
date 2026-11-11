# Distribution functions for Ranking SDT

Density and random generation for ranking signal detection theory
(Meyer-Grant et al., 2026). Models rank ordering of m items by perceived
strength. Only `d` is estimated (no criterion or stimulus column).
Supports Gumbel-min (closed form) and Gaussian UV-SDT (numerical
integration).

## Usage

``` r
dsdt_ranking(
  counts,
  m,
  d,
  sdratio = 1,
  dist = c("gumbel_min", "normal"),
  log = FALSE
)

rsdt_ranking(n, n_trials, m, d, sdratio = 1, dist = c("gumbel_min", "normal"))
```

## Arguments

- counts:

  Integer matrix with one row per observation and one rank-count column
  per rank position (1 = most likely target), or a vector for a single
  observation. A data frame is read as the matrix of its rows. Columns
  beyond a row's set size `m` must be 0.

- m:

  Integer vector. Number of ranked items per observation. Must be at
  least 2 and no larger than the number of count columns.

- d:

  Numeric vector. Sensitivity: the distance between the target and lure
  distributions. It is \\d'\\ when `sdratio` is 1 and, for
  `dist = "gumbel_min"`, the \\g'\\ of Meyer-Grant et al. (2026). With
  another `sdratio` it is the balanced index \\d_a\\ that
  [`sdt_yn()`](https://venpopov.com/bmm/dev/reference/sdt_yn.md) reports
  (in root-mean-square SD units).

- sdratio:

  Numeric vector. Ratio of signal to noise standard deviations (default
  1, i.e., equal variance). Must be positive. Only used when
  `dist = "normal"`.

- dist:

  Character. The distribution assumed for the latent evidence:
  "gumbel_min" (default), the smallest extreme value distribution with
  cumulative distribution function \\1 - \exp(-\exp(x))\\, evaluated in
  closed form; or "normal", Gaussian UV-SDT by numerical integration.

- log:

  Logical. If `TRUE`, returns log-density (default `FALSE`).

- n:

  Integer. Number of observations to generate. `n_trials`, `m`, `d`, and
  `sdratio` are recycled to this length.

- n_trials:

  Integer vector. Number of ranking trials per observation.

## Value

`dsdt_ranking` returns the (log-)density (multinomial probability).
`rsdt_ranking` returns an integer matrix with one row per observation
and one rank-count column per rank position, `rank1` to `rank<M>`, where
`M` is the largest `m` of the draws; with `n = 0` it is the largest of
all `m` given. Rows with a smaller set size have structural zeros in the
surplus columns, matching the wide format
[`sdt_ranking()`](https://venpopov.com/bmm/dev/reference/sdt_ranking.md)
expects.

## Parameter scales

These functions take `sdratio` as the ratio itself, while
[`sdt_ranking()`](https://venpopov.com/bmm/dev/reference/sdt_ranking.md)
estimates its logarithm (0 = equal variance): a fitted `sdratio` of
0.375 is a ratio of `exp(0.375) = 1.455`, and passing `0.375` here
instead asks for a signal distribution 2.7 times narrower than the
noise. That is a legal value and raises no error, so exponentiate first.
`d` carries across unchanged.

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
  `dsdt_ranking()` and
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

Meyer-Grant, C. G., Kellen, D., Harding, S. M., & Singmann, H. (2026).
Extreme-value signal detection theory for recognition memory: The
parametric road not taken. *Psychological Review*. Advance online
publication.
[doi:10.1037/rev0000615](https://doi.org/10.1037/rev0000615)

## Examples

``` r
# Gumbel-min ranking density
dsdt_ranking(counts = c(40, 30, 20, 10), m = 4, d = 1.0)
#> [1] 4.92936e-06
# Generate ranking data (m=4, Gumbel-min) for 10 subjects
dat <- data.frame(id = 1:10, set_size = 4L)
dat <- cbind(dat, rsdt_ranking(10, 100, m = 4, d = 1.0))
head(dat)
#>   id set_size rank1 rank2 rank3 rank4
#> 1  1        4    53    25    11    11
#> 2  2        4    63    14    10    13
#> 3  3        4    60    17    11    12
#> 4  4        4    56    16    21     7
#> 5  5        4    54    23    18     5
#> 6  6        4    53    27    10    10
```
