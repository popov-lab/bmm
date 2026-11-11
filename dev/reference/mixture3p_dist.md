# Distribution functions for the three-parameter mixture model (mixture3p)

Density, distribution, and random generation functions for the
three-parameter mixture model with the location of `mu`, precision of
memory representations `kappa`, probability of a target response
`p_mem`, and probability of a non-target response `p_nt`.

## Usage

``` r
dmixture3p(
  x,
  mu = c(0, 2, -1.5),
  kappa = 5,
  p_mem = 0.6,
  p_nt = 0.2,
  tau = 0,
  vp_nodes = 81L,
  log = FALSE
)

rmixture3p(n, mu = c(0, 2, -1.5), kappa = 5, p_mem = 0.6, p_nt = 0.2, tau = 0)
```

## Arguments

- x:

  Vector of observed responses

- mu:

  Vector of locations. First value represents the location of the target
  item and any additional values indicate the location of non-target
  items.

- kappa:

  Vector of precision values. With `tau > 0` this is the concentration
  of the average trial, whose Fisher information `J(kappa)` is the mean
  of `J`.

- p_mem:

  Vector of probabilities of a target response. For a
  [`mixture3p()`](https://venpopov.com/bmm/dev/reference/mixture3p.md)
  fit these are `thetat` on its native scale (version `"simple"`, see
  [`native_parameters()`](https://venpopov.com/bmm/dev/reference/native_parameters.md))
  or `min(1, K / set_size) * ptarget` (versions `"slot"` and
  `"slot_averaging"`).

- p_nt:

  Vector of probabilities of a non-target response, summed over the
  non-targets, which share it equally. For a
  [`mixture3p()`](https://venpopov.com/bmm/dev/reference/mixture3p.md)
  fit these are `thetant` on its native scale (version `"simple"`) or
  `min(1, K / set_size) * (1 - ptarget)` (capacity versions). Unlike
  `ptarget`, `p_mem` and `p_nt` are not conditional on a response from
  memory.

- tau:

  Vector of scales for trial-to-trial variability in precision. The
  Fisher information `J` of a memory representation is drawn from
  `gamma(mean = J(kappa), scale = tau)`. `tau = 0`, the default, is
  constant precision. See
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

## References

Bays, P. M., Catalao, R. F. G., & Husain, M. (2009). The precision of
visual working memory is set by allocation of a shared resource. Journal
of Vision, 9(10), Article 7.

van den Berg, R., Shin, H., Chou, W.-C., George, R., & Ma, W. J. (2012).
Variability in encoding precision accounts for visual short-term memory
limitations. PNAS, 109(22), 8780-8785.

## Examples

``` r
# generate random samples from the mixture3p model and overlay the density
r <- rmixture3p(10000, mu = c(0, 2, -1.5), kappa = 4, p_mem = 0.6, p_nt = 0.2)
x <- seq(-pi, pi, length.out = 10000)
d <- dmixture3p(x, mu = c(0, 2, -1.5), kappa = 4, p_mem = 0.6, p_nt = 0.2)
hist(r, breaks = 60, freq = FALSE)
lines(x, d, type = "l", col = "red")

```
