# Distribution functions for Continuous Dual-Process SDT (CDP)

Density and random generation for the continuous dual-process signal
detection theory model (Wixted & Mickes, 2010). Two correlated
continuous dimensions, Familiarity (F) and Recollection (R), generate
the aggregate strength S = F + R that drives old/new confidence;
Remember/Know judgments split "old" responses on R, and an optional
Know/Guess split uses F. The response is a vector of counts across the
response categories (multinomial likelihood), in the same format
[`sdt_cdp()`](https://venpopov.com/bmm/dev/reference/sdt_cdp.md) is fit
to. See [`sdt_cdp()`](https://venpopov.com/bmm/dev/reference/sdt_cdp.md)
for the model and
[`sdt_rating()`](https://venpopov.com/bmm/dev/reference/sdt_rating.md)
when no R/K split is available.

## Usage

``` r
dsdt_cdp(
  counts,
  stimulus,
  dfam,
  drec,
  thresholds,
  rcrit,
  kcrit = NULL,
  sigmar = 0,
  rho = 0,
  n_new = NULL,
  dist = "normal",
  log = FALSE
)

rsdt_cdp(
  n,
  n_trials,
  stimulus,
  dfam,
  drec,
  thresholds,
  rcrit,
  kcrit = NULL,
  sigmar = 0,
  rho = 0,
  n_new = NULL,
  dist = "normal"
)
```

## Arguments

- counts:

  Integer matrix with one row per observation and one column per
  response category in the canonical order `new`, `[guess]`, `know`,
  `remember` (each block ordered by confidence; see
  [`sdt_cdp()`](https://venpopov.com/bmm/dev/reference/sdt_cdp.md)), or
  a vector for a single observation. A data frame is read as the matrix
  of its rows. `n_new + 2 * n_old` columns (R/K) or `n_new + 3 * n_old`
  (R/K/G).

- stimulus:

  Integer vector (0/1). Stimulus type: 0 = new/lure, 1 = old/target.

- dfam, drec:

  Numeric vectors. Familiarity and recollection sensitivities: the
  target means on the F and R axes, each in units of the corresponding
  lure SD (both lure SDs are 1). They are component means of a bivariate
  latent space, not the balanced \\d_a\\ that the other SDT models
  report as `d` – the old/new decision here is read off the aggregate
  strength S = F + R, so the model's discriminability is a derived
  quantity rather than either of these.

- thresholds:

  Numeric vector of `n_new + n_old - 1` ordered confidence thresholds on
  the aggregate F+R axis, or a matrix with one row per observation.

- rcrit:

  Numeric vector. Remember criterion on the recollection axis.

- kcrit:

  Know criterion on the familiarity axis. `NULL` (default), or a
  non-finite value for every observation, gives the R/K model; finite
  values for every observation give the R/K/G model. A `kcrit` with any
  non-finite element gives the R/K model for the whole call, also for
  the observations whose `kcrit` is finite, and warns when there are
  such observations. Every element counts, also one past the last
  observation: with one draw, `kcrit = c(0.2, Inf)` gives R/K and a
  warning. The length of `kcrit` sets the number of observations only
  when every element is finite (the R/K/G model): with one row of
  `counts`, `kcrit = c(0.2, 0.3)` gives two values and
  `kcrit = c(0.2, Inf)` one.

- sigmar:

  Numeric vector. Log SD of the recollection target distribution (0 = SD
  1).

- rho:

  Numeric vector. F-R correlation on the unconstrained scale;
  `tanh(rho)` is the correlation. 0 (default) = independent processes
  (classic CDP).

- n_new:

  Integer number of "new" confidence levels. Defaults to half the number
  of confidence levels, rounded down (so 3 of 7); the number of "old"
  levels follows as `length(thresholds) + 1 - n_new`.

- dist:

  Noise distribution. Only `"normal"` is currently supported.

- log:

  Logical; if `TRUE` return the log-density (default `FALSE`).

- n:

  Integer. Number of observations to generate. `n_trials`, `stimulus`,
  `thresholds`, and the model parameters are recycled to this length.

- n_trials:

  Integer vector. Number of trials per observation.

## Value

`dsdt_cdp` returns the (log-)density (multinomial probability).
`rsdt_cdp` returns an integer matrix with one row per observation and
the canonical response count columns (`new1`, ..., `remember<K>`) that
[`sdt_cdp()`](https://venpopov.com/bmm/dev/reference/sdt_cdp.md) expects
– [`cbind()`](https://rdrr.io/r/base/cbind.html) it to a design data
frame for a ready-to-fit data set.

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
  and `dsdt_cdp()`), or the pair of `mean_rt` or `var_rt` values of
  `dezdm(version = "4par")`, a vector is always one observation: pass a
  matrix with zero rows for none, as an empty vector stops. An empty
  `kcrit` of `dsdt_cdp()` and `rsdt_cdp()` stops as well; `kcrit = NULL`
  gives the R/K model. Counts in a data frame, a tibble included, are
  read as the matrix of its rows, and so are `mean_rt` and `var_rt` of
  `dezdm(version = "4par")`.
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

Wixted, J. T., & Mickes, L. (2010). A continuous dual-process model of
remember/know judgments. *Psychological Review*, *117*(4), 1025–1054.
[doi:10.1037/a0020874](https://doi.org/10.1037/a0020874)

## Examples

``` r
# CDP density (R/K, 1 new + 3 old levels: 3 thresholds, 7 count columns)
dsdt_cdp(
  counts = c(40, 5, 12, 30, 10, 18, 60), stimulus = 1,
  dfam = 0.8, drec = 1.0,
  thresholds = c(-0.5, 0.3, 1.0),
  rcrit = 0.7, n_new = 1
)
#> [1] 1.870555e-27
# Generate CDP count data (R/K, 3 new + 3 old levels) for 10 subjects
dat <- expand.grid(id = 1:10, stimulus = c(0L, 1L))
dat <- cbind(dat, rsdt_cdp(nrow(dat), 100, dat$stimulus,
                           dfam = 0.8, drec = 1.0,
                           thresholds = c(-1.1, -0.5, 0, 0.6, 1.3),
                           rcrit = 0.5, n_new = 3))
head(dat)
#>   id stimulus new1 new2 new3 know4 know5 know6 remember4 remember5 remember6
#> 1  1        0   13   16   14    16     5     5         8         9        14
#> 2  2        0   16   13   10    16     5     8         7        14        11
#> 3  3        0   26   17   14     8    11     1         1        11        11
#> 4  4        0   18   18   11    11     9     3         3         8        19
#> 5  5        0   19   11   12    16     9     5         6        11        11
#> 6  6        0   18    9   22    13    10     6         2         6        14
```
