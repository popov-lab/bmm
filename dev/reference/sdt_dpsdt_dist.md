# Distribution functions for dual-process SDT (DPSDT)

Density and random generation for the dual-process signal detection
model (Yonelinas, 1994), with recall-to-reject of new items (`Rn`; as in
Yonelinas, 2024). Extends rating SDT with recollection probabilities
`Ro` (old items recollected as old) and `Rn` (new items recall-rejected)
that add mass to the most-confident rating category. These are the
simulation counterparts of the `dpsdt` version of
[`sdt_rating()`](https://venpopov.com/bmm/dev/reference/sdt_rating.md);
here `Ro`/`Rn` are supplied directly as probabilities in `[0, 1]`.

## Usage

``` r
dsdt_dpsdt(
  counts,
  stimulus,
  d,
  thresholds,
  Ro,
  Rn = 0,
  sdratio = 1,
  dist = c("normal", "gumbel_min", "gumbel_max", "logistic"),
  log = FALSE
)

rsdt_dpsdt(
  n,
  n_trials,
  stimulus,
  d,
  thresholds,
  Ro,
  Rn = 0,
  sdratio = 1,
  dist = c("normal", "gumbel_min", "gumbel_max", "logistic")
)
```

## Arguments

- counts:

  Integer matrix with one row per observation and one column per rating
  category, ordered from "definitely noise" (1) to "definitely signal"
  (K), or a vector for a single observation. A data frame is read as the
  matrix of its rows.

- stimulus:

  Integer vector (0/1). Stimulus type: 0 = noise, 1 = signal.

- d:

  Numeric vector. Sensitivity: \\d'\\ when `sdratio` is 1, and otherwise
  the balanced index \\d_a\\ (see
  [`sdt_yn()`](https://venpopov.com/bmm/dev/reference/sdt_yn.md)). The
  separation between the distributions in noise units is
  `d * sqrt((1 + sdratio^2) / 2)`.

- thresholds:

  Numeric vector of length K-1 with the ordered decision thresholds, or
  an n-by-(K-1) matrix with one row per observation. The thresholds are
  positions on the noise-standardized axis, measured from the mean of
  the noise distribution (the signal mean sits at
  `d * sqrt((1 + sdratio^2) / 2)`), and are not rescaled by `sdratio`.

- Ro:

  Numeric vector in `[0, 1]`. Recollection probability for old (signal)
  items.

- Rn:

  Numeric vector in `[0, 1]`. Recall-to-reject probability for new
  (noise) items. Defaults to 0, the classic one-sided model, as the
  `dpsdt` version fixes it off unless `Rn` is in the formula.

- sdratio:

  Numeric vector. Ratio of signal to noise standard deviations (default
  1, i.e., equal variance). Must be positive. Note that this is the
  natural scale: the `sdratio` parameter of
  [`sdt_rating()`](https://venpopov.com/bmm/dev/reference/sdt_rating.md)
  is sampled on the log scale, so it corresponds to `log(sdratio)` here.

- dist:

  Character. Noise distribution: "normal" (default), "logistic",
  "gumbel_min", or "gumbel_max".

- log:

  Logical. If `TRUE`, returns log-density (default `FALSE`).

- n:

  Integer. Number of observations to generate. `n_trials`, `stimulus`,
  `thresholds`, and the model parameters are recycled to this length.

- n_trials:

  Integer vector. Number of trials per observation.

## Value

`dsdt_dpsdt` returns the (log-)density (multinomial probability).
`rsdt_dpsdt` returns an integer matrix with one row per observation and
one rating-count column per category (`r1` ... `rK`).

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
  `dsdt_dpsdt()`,
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

Yonelinas, A. P. (1994). Receiver-operating characteristics in
recognition memory: Evidence for a dual-process model. *Journal of
Experimental Psychology: Learning, Memory, and Cognition*, *20*(6),
1341–1354.
[doi:10.1037/0278-7393.20.6.1341](https://doi.org/10.1037/0278-7393.20.6.1341)

Yonelinas, A. P. (2024). The role of recollection and familiarity in
visual working memory: A mixture of threshold and signal detection
processes. *Psychological Review*, *131*(2), 321–348.
[doi:10.1037/rev0000432](https://doi.org/10.1037/rev0000432)

## Examples

``` r
# Density for a single observation (K=4) with recollection of old items;
# Rn defaults to 0 (no recall-to-reject), the one-sided model
dsdt_dpsdt(counts = c(2, 8, 20, 70), stimulus = 1,
           d = 1.5, thresholds = c(0.25, 0.75, 1.25), Ro = 0.3)
#> [1] 3.023007e-05
# Generate DPSDT rating data (K=4) for 10 subjects and both stimulus types
dat <- expand.grid(id = 1:10, stimulus = c(0L, 1L))
dat <- cbind(dat, rsdt_dpsdt(nrow(dat), 100, dat$stimulus, d = 1.5,
                             thresholds = c(0.25, 0.75, 1.25),
                             Ro = 0.3, Rn = 0.1))
head(dat)
#>   id stimulus r1 r2 r3 r4
#> 1  1        0 64 19  8  9
#> 2  2        0 67 18  8  7
#> 3  3        0 73  7 10 10
#> 4  4        0 63 16 12  9
#> 5  5        0 68 15 11  6
#> 6  6        0 59 17 11 13
```
