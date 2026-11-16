# Distribution functions for the three-parameter mixture model (mixture3p)

Density, distribution, and random generation functions for the
three-parameter mixture model with the location of `mu`, precision of
memory representations `kappa`, probability of a target response
`p_mem`, and probability of a non-target response `p_nt`. The capacity
versions (`version = "slot"` and `"slot_averaging"`) replace `p_mem` and
`p_nt` by the capacity `K` and `ptarget`, as
[`mixture3p()`](https://venpopov.com/bmm/dev/reference/mixture3p.md)
does.

## Usage

``` r
dmixture3p(
  x,
  mu = c(0, 2, -1.5),
  kappa = 5,
  p_mem = 0.6,
  p_nt = 0.2,
  tau = 0,
  K = 3,
  ptarget = 0.75,
  version = "simple",
  vp_nodes = 81L,
  log = FALSE
)

rmixture3p(
  n,
  mu = c(0, 2, -1.5),
  kappa = 5,
  p_mem = 0.6,
  p_nt = 0.2,
  tau = 0,
  K = 3,
  ptarget = 0.75,
  version = "simple"
)
```

## Arguments

- x:

  Vector of observed responses

- mu:

  Item locations in radians. A vector is one item set shared by every
  observation: the target first, then the non-targets. A matrix or data
  frame holds one item set per observation in each row: the target in
  the first column, then the row's non-targets, and `NA` past its set
  size, so rows can differ in set size. The rows are recycled to the
  number of observations. The random generation function draws the rows
  of one set size at a time, so with a fixed seed, moving a row can
  change which draw it gets. Before bmm 1.5.0 a matrix was read as one
  long item set. The number of locations of a row, or of a vector,
  before any `NA` padding is the set size the `"slot"` and
  `"slot_averaging"` versions use, so pass the target and every
  non-target of the trial: for a trial of a fitted model, the target
  followed by the first `set_size - 1` columns of `nt_features`.

- kappa:

  Vector of precision values. With `tau > 0` this is the concentration
  of the average trial, whose Fisher information `J(kappa)` is the mean
  of `J`. Under `"slot_averaging"` it is the precision of an item held
  in a single slot.

- p_mem:

  Vector of probabilities of a target response, i.e. `thetat` of a
  [`mixture3p()`](https://venpopov.com/bmm/dev/reference/mixture3p.md)
  fit with `version = "simple"` on its native scale (see
  [`native_parameters()`](https://venpopov.com/bmm/dev/reference/native_parameters.md)).
  Used by the `"simple"` version; the other versions derive the
  probability of a response from memory from `K` and the set size, and
  warn if `p_mem` is given.

- p_nt:

  Vector of probabilities of a non-target response, summed over the
  non-targets, which share it equally: `thetant` of a `"simple"` fit on
  its native scale. Unlike `ptarget`, `p_mem` and `p_nt` are not
  conditional on a response from memory. Used by the `"simple"` version,
  like `p_mem`.

- tau:

  Vector of scales for trial-to-trial variability in precision. The
  Fisher information `J` of a memory representation is drawn from
  `gamma(mean = J(kappa), scale = tau)`. `tau = 0`, the default, is
  constant precision. See
  [`mixture3p()`](https://venpopov.com/bmm/dev/reference/mixture3p.md).

- K:

  Vector of memory capacities. Used by the `"slot"` and
  `"slot_averaging"` versions, which read each observation's set size as
  the number of locations in its item set in `mu` (the row's locations
  before the `NA` padding); `"simple"` ignores `K` and warns if it is
  given.

- ptarget:

  Vector of probabilities that a response from memory reports the target
  rather than a non-target. Used by the `"slot"` and `"slot_averaging"`
  versions; `"simple"` ignores it and warns if it is given.

- version:

  A character string selecting the storage model: `"simple"`, `"slot"`
  or `"slot_averaging"`. See
  [`mixture3p()`](https://venpopov.com/bmm/dev/reference/mixture3p.md).

- vp_nodes:

  Number of quadrature nodes used when `tau > 0`; must be an odd number
  of at least 41. See
  [`mixture2p()`](https://venpopov.com/bmm/dev/reference/mixture2p.md).
  The random generation function needs no nodes: it draws `J` from its
  gamma distribution and then the response, which is exact.

- log:

  Logical; if `TRUE`, values are returned on the log scale. Pass it by
  name: in bmm 1.4 it was the sixth argument, which is now `tau`. Until
  bmm 1.6.0 a logical value in that position is still read as `log`,
  with a warning.

- n:

  Number of observations to generate data for

## Value

`dmixture3p` gives the density of the three-parameter mixture model, and
`rmixture3p` gives the random generation function for the
three-parameter mixture model.

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
  `dist`) of `dmixture3p()` and
  [`dimm()`](https://venpopov.com/bmm/dev/reference/IMMdist.md), is a
  vector when every observation shares it and a matrix with one row per
  observation otherwise: row `i` belongs to observation `i`, or to draw
  `i`. Item locations are padded with `NA` past each row's set size, and
  a data frame of them is read as the matrix of its rows. As for counts,
  an empty vector of them stops: pass a matrix with zero rows for no
  observations.

- Configuration arguments such as `version`, `dist` of the ranking
  models, `log`, `parametrization` and `n_new` are single values and are
  never recycled. An argument that a version does not use, such as `zr`
  of `dezdm(version = "3par")`, is ignored.

## References

Bays, P. M., Catalao, R. F. G., & Husain, M. (2009). The precision of
visual working memory is set by allocation of a shared resource. Journal
of Vision, 9(10), Article 7.

van den Berg, R., Shin, H., Chou, W.-C., George, R., & Ma, W. J. (2012).
Variability in encoding precision accounts for visual short-term memory
limitations. PNAS, 109(22), 8780-8785.

Zhang, W., & Luck, S. J. (2008). Discrete fixed-resolution
representations in visual working memory. Nature, 453(7192), 233-235.

## Examples

``` r
# generate random samples from the mixture3p model and overlay the density
r <- rmixture3p(10000, mu = c(0, 2, -1.5), kappa = 4, p_mem = 0.6, p_nt = 0.2)
x <- seq(-pi, pi, length.out = 10000)
d <- dmixture3p(x, mu = c(0, 2, -1.5), kappa = 4, p_mem = 0.6, p_nt = 0.2)
hist(r, breaks = 60, freq = FALSE, ylim = c(0, 0.9))
lines(x, d, type = "l", col = "red")

# slot averaging: with 6 slots for 3 items each item holds two slots
d6 <- dmixture3p(x, mu = c(0, 2, -1.5), kappa = 4, K = 6, ptarget = 0.8,
                 version = "slot_averaging")
lines(x, d6, col = "blue")


# item locations per trial: one row per trial, NA past its set size
locations <- rbind(c(0, 2, -1.5), c(1, -2, NA), c(-0.5, NA, NA))
y <- rmixture3p(3, mu = locations, kappa = 4, p_mem = 0.6, p_nt = 0.2)
dmixture3p(y, mu = locations, kappa = 4, p_mem = 0.6, p_nt = 0.2)
#> [1] 0.1379080 0.4523002 0.4898005
```
