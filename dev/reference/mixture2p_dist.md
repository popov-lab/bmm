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

pmixture2p(q, mu = 0, kappa = 7, p_mem = 0.8)

qmixture2p(p, mu = 0, kappa = 5, p_mem = 0.6)

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

  Vector of locations

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

- q:

  Vector of quantiles

- p:

  Vector of probability

- n:

  Number of observations to generate data for

## Value

`dmixture2p` gives the density of the two-parameter mixture model,
`pmixture2p` gives the cumulative distribution function of the
two-parameter mixture model, `qmixture2p` gives the quantile function of
the two-parameter mixture model, and `rmixture2p` gives the random
generation function for the two-parameter mixture model.

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
