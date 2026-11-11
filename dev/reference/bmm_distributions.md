# Vector arguments of the distribution functions

The density (`d*()`) and random generation (`r*()`) functions of the
sdm, ddm and ezdm models and of the signal detection models (`sdt_yn`,
`sdt_mafc`, `sdt_rating`, `sdt_dpsdt`, `sdt_metad`, `sdt_ranking` and
`sdt_cdp`) take their per-observation arguments as vectors, or as
matrices with one row per observation, and recycle them as
[`stats::dnorm()`](https://rdrr.io/r/stats/Normal.html) and
[`stats::rnorm()`](https://rdrr.io/r/stats/Normal.html) do. This page
states the rules; the help page of each of these distributions repeats
them under "Vector arguments".

The other distribution functions do not follow these rules yet:
[`dcswald()`](https://venpopov.com/bmm/dev/reference/cswald_dist.md) and
[`rcswald()`](https://venpopov.com/bmm/dev/reference/cswald_dist.md),
the circular mixtures
([`dmixture2p()`](https://venpopov.com/bmm/dev/reference/mixture2p_dist.md),
[`dmixture3p()`](https://venpopov.com/bmm/dev/reference/mixture3p_dist.md),
[`dimm()`](https://venpopov.com/bmm/dev/reference/IMMdist.md) and their
`r*()` functions),
[`dm3()`](https://venpopov.com/bmm/dev/reference/m3dist.md) and
[`rm3()`](https://venpopov.com/bmm/dev/reference/m3dist.md),
[`dmpt()`](https://venpopov.com/bmm/dev/reference/mptdist.md) and
[`rmpt()`](https://venpopov.com/bmm/dev/reference/mptdist.md), and
[`dpsychometric_function()`](https://venpopov.com/bmm/dev/reference/psychometric_function_dist.md).
[`rpsychometric_function()`](https://venpopov.com/bmm/dev/reference/psychometric_function_dist.md)
reads `n` by the rule below but does not yet recycle its other arguments
by these rules.

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
