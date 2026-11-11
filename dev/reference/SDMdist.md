# Distribution functions for the Signal Discrimination Model (SDM)

Density, distribution function, and random generation for the Signal
Discrimination Model (SDM) Distribution with location `mu`, memory
strength `c`, and precision `kappa`. Currently only a single activation
source is supported.

## Usage

``` r
dsdm(x, mu = 0, c = 3, kappa = 3.5, log = FALSE, parametrization = "sqrtexp")

psdm(
  q,
  mu = 0,
  c = 3,
  kappa = 3.5,
  lower.tail = TRUE,
  log.p = FALSE,
  lower.bound = -pi,
  parametrization = "sqrtexp"
)

qsdm(
  p,
  mu = 0,
  c = 3,
  kappa = 3.5,
  lower.tail = TRUE,
  log.p = FALSE,
  lower.bound = -pi,
  parametrization = "sqrtexp"
)

rsdm(n, mu = 0, c = 3, kappa = 3.5, parametrization = "sqrtexp")
```

## Arguments

- x:

  Vector of quantiles, or a data frame of them. For a data frame,
  `dsdm()` returns a data frame of the same shape with one density per
  cell, its cells taken column by column; parameters longer than the
  number of cells are cut to it, and an empty parameter gives the data
  frame with no rows.

- mu:

  Vector of location values in radians

- c:

  Vector of memory strength values

- kappa:

  Vector of precision values

- log:

  Logical; if `TRUE`, values are returned on the log scale.

- parametrization:

  Character; either `"bessel"` or `"sqrtexp"` (default). See [the online
  article](https://venpopov.com/bmm/articles/bmm_sdm_simple.html) for
  details on the parameterization.

- q:

  Vector of quantiles

- lower.tail:

  Logical; If `TRUE` (default), return P(X \<= x). Else, return P(X \>
  x)

- log.p:

  Logical; if `TRUE`, probabilities are returned on the log scale.

- lower.bound:

  Numeric; Lower bound of integration for the cumulative distribution

- p:

  Vector of probabilities

- n:

  Number of observations to sample

## Value

`dsdm` gives the density, `psdm` gives the distribution function, `qsdm`
gives the quantile function, and `rsdm` generates random deviates.

## Details

**Parametrization**

See [the online
article](https://venpopov.com/bmm/articles/bmm_sdm_simple.html) for
details on the parameterization. Oberauer (2023) introduced the SDM with
the bessel parametrization. The sqrtexp parametrization is the default
in the `bmm` package for numerical stability and efficiency. The two
parametrizations are related by the functions
[`c_bessel2sqrtexp()`](https://venpopov.com/bmm/dev/reference/c_parametrizations.md)
and
[`c_sqrtexp2bessel()`](https://venpopov.com/bmm/dev/reference/c_parametrizations.md).

**The cumulative distribution function**

Since responses are on the circle, the cumulative distribution function
requires you to choose a lower bound of integration. The default is
\\-\pi\\, as for the brms::pvon_mises() function but you can choose any
value in the argument `lower.bound` of `psdm` and `qsdm`. Another useful
choice is the mean of the response distribution minus \\\pi\\, e.g.
`lower.bound = mu-pi`. This is the default in `circular::pvonmises()`,
and it ensures that 50% of the cumulative probability mass is below the
mean of the response distribution.

`psdm` is 0 at `lower.bound` and 1 at `lower.bound + 2 * pi`; a `q`
outside this range is wrapped around the circle. `qsdm` inverts `psdm`
by root finding and returns values between `lower.bound` and
`lower.bound + 2 * pi`. Values of `p` outside \[0, 1\] give `NaN` with a
warning. Missing values give `NA`; `NaN`, infinite values and positions
(`q`, `mu`, `lower.bound`) beyond 1e8 in absolute value give `NaN`, as
does the bessel parametrization for `kappa` above 1e5.

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
  are `mean_rt` and `var_rt` of `dezdm(version = "4par")`. `dsdm()` with
  a data frame `x` is the exception: it returns a data frame of the same
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

Oberauer, K. (2023). Measurement models for visual working memory - A
factorial model comparison. Psychological Review, 130(3), 841–852

## Examples

``` r
# plot the density of the SDM distribution
x <- seq(-pi, pi, length.out = 10000)
plot(x, dsdm(x, 0, 2, 3),
  type = "l", xlim = c(-pi, pi), ylim = c(0, 1),
  xlab = "Angle error (radians)",
  ylab = "density",
  main = "SDM density"
)
lines(x, dsdm(x, 0, 9, 1), col = "red")
lines(x, dsdm(x, 0, 2, 8), col = "green")
legend("topright", c(
  "c=2, kappa=3.0, mu=0",
  "c=9, kappa=1.0, mu=0",
  "c=2, kappa=8, mu=1"
),
col = c("black", "red", "green"), lty = 1, cex = 0.8
)


# plot the cumulative distribution function of the SDM distribution
p <- psdm(x, mu = 0, c = 3.1, kappa = 5)
plot(x, p, type = "l")


# generate random deviates from the SDM distribution and overlay the density
r <- rsdm(10000, mu = 0, c = 3.1, kappa = 5)
d <- dsdm(x, mu = 0, c = 3.1, kappa = 5)
hist(r, breaks = 60, freq = FALSE)
lines(x, d, type = "l", col = "red")
```
