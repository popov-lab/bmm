# Distribution functions for the two-parameter mixture model (mixture2p)

Density, distribution, and random generation functions for the
two-parameter mixture model with the location of `mu`, precision of
memory representations `kappa` and probability of recalling items from
memory `p_mem`.

## Usage

``` r
dmixture2p(
  x,
  mu = 0,
  kappa = 5,
  p_mem = 0.6,
  tau = 0,
  K = 3,
  set_size = NULL,
  version = "simple",
  vp_nodes = 81L,
  log = FALSE
)

rmixture2p(
  n,
  mu = 0,
  kappa = 5,
  p_mem = 0.6,
  tau = 0,
  K = 3,
  set_size = NULL,
  version = "simple"
)
```

## Arguments

- x:

  Vector of observed responses

- mu:

  Vector of target locations, one per observation. Unlike `mu` of
  [`dmixture3p()`](https://venpopov.com/bmm/dev/reference/mixture3p_dist.md)
  and [`dimm()`](https://venpopov.com/bmm/dev/reference/IMMdist.md), it
  holds no non-targets.

- kappa:

  Vector of precision values. With `tau > 0` this is the concentration
  of the average trial, whose Fisher information `J(kappa)` is the mean
  of `J`.

- p_mem:

  Vector of probabilities for memory recall, i.e. `thetat` of
  [`mixture2p()`](https://venpopov.com/bmm/dev/reference/mixture2p.md)
  on its native scale. Used by the `"simple"` version; the other
  versions derive it from `K` and `set_size`, and warn if `p_mem` is
  given.

- tau:

  Vector of scales for trial-to-trial variability in precision. The
  Fisher information `J` of a memory representation is drawn from
  `gamma(mean = J(kappa), scale = tau)`. `tau = 0`, the default, is
  constant precision.

- K:

  Vector of memory capacities. Used by the `"slot"` and
  `"slot_averaging"` versions; `"simple"` ignores it and warns if it is
  given.

- set_size:

  Vector of set sizes. Required by the `"slot"` and `"slot_averaging"`
  versions, which stop without it; `"simple"` ignores it and warns if it
  is given.

- version:

  A character string selecting the storage model: `"simple"`, `"slot"`
  or `"slot_averaging"`. See
  [`mixture2p()`](https://venpopov.com/bmm/dev/reference/mixture2p.md).

- vp_nodes:

  Number of quadrature nodes used when `tau > 0`; must be an odd number
  of at least 41. See
  [`mixture2p()`](https://venpopov.com/bmm/dev/reference/mixture2p.md).
  The random generation function needs no nodes: it draws `J` from its
  gamma distribution and then the response, which is exact.

- log:

  Logical; if `TRUE`, values are returned on the log scale. Pass it by
  name: in bmm 1.4 it was the fifth argument, which is now `tau`. Until
  bmm 1.6.0 a logical value in that position is still read as `log`,
  with a warning.

- n:

  Number of observations to generate data for

## Value

`dmixture2p` gives the density of the two-parameter mixture model, and
`rmixture2p` gives the random generation function for the two-parameter
mixture model.

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
  `thresholds` of the rating models or the item locations `mu` (and
  `dist`) of
  [`dmixture3p()`](https://venpopov.com/bmm/dev/reference/mixture3p_dist.md)
  and [`dimm()`](https://venpopov.com/bmm/dev/reference/IMMdist.md), is
  a vector when every observation shares it and a matrix with one row
  per observation otherwise: row `i` belongs to observation `i`, or to
  draw `i`. Item locations are padded with `NA` past each row's set
  size, and a data frame of them is read as the matrix of its rows. As
  for counts, an empty vector of them stops: pass a matrix with zero
  rows for no observations.

- Configuration arguments such as `version`, `dist` of the ranking
  models, `log`, `parametrization` and `n_new` are single values and are
  never recycled. An argument that a version does not use, such as `zr`
  of `dezdm(version = "3par")`, is ignored.

## References

Zhang, W., & Luck, S. J. (2008). Discrete fixed-resolution
representations in visual working memory. Nature, 453(7192), 233-235.

van den Berg, R., Shin, H., Chou, W.-C., George, R., & Ma, W. J. (2012).
Variability in encoding precision accounts for visual short-term memory
limitations. PNAS, 109(22), 8780-8785.

## Examples

``` r
# generate random samples from the mixture2p model and overlay the density
r <- rmixture2p(10000, mu = 0, kappa = 4, p_mem = 0.8)
x <- seq(-pi, pi, length.out = 10000)
d <- dmixture2p(x, mu = 0, kappa = 4, p_mem = 0.8)
hist(r, breaks = 60, freq = FALSE)
lines(x, d, type = "l", col = "red")

```
