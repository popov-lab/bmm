# Distribution functions for the three-parameter mixture model (mixture3p)

Density, distribution, and random generation functions for the
three-parameter mixture model with the location of `mu`, precision of
memory representations `kappa`, probability of recalling items from
memory `p_mem`, and probability of recalling non-targets `p_nt`.

## Usage

``` r
dmixture3p(
  x,
  mu = c(0, 2, -1.5),
  kappa = 5,
  p_mem = 0.6,
  p_nt = 0.2,
  log = FALSE
)

pmixture3p(q, mu = c(0, 2, -1.5), kappa = 5, p_mem = 0.6, p_nt = 0.2)

qmixture3p(p, mu = c(0, 2, -1.5), kappa = 5, p_mem = 0.6, p_nt = 0.2)

rmixture3p(n, mu = c(0, 2, -1.5), kappa = 5, p_mem = 0.6, p_nt = 0.2)
```

## Arguments

- x:

  Vector of observed responses

- mu:

  Vector of locations. First value represents the location of the target
  item and any additional values indicate the location of non-target
  items.

- kappa:

  Vector of precision values

- p_mem:

  Vector of probabilities for memory recall

- p_nt:

  Vector of probabilities for swap errors

- log:

  Logical; if `TRUE`, values are returned on the log scale.

- q:

  Vector of quantiles

- p:

  Vector of probability

- n:

  Number of observations to generate data for

## Value

`dmixture3p` gives the density of the three-parameter mixture model,
`pmixture3p` gives the cumulative distribution function of the
two-parameter mixture model, `qmixture3p` gives the quantile function of
the two-parameter mixture model, and `rmixture3p` gives the random
generation function for the two-parameter mixture model.

## References

Bays, P. M., Catalao, R. F. G., & Husain, M. (2009). The precision of
visual working memory is set by allocation of a shared resource. Journal
of Vision, 9(10), 7.

## Examples

``` r
# generate random samples from the mixture3p model and overlay the density
r <- rmixture3p(10000, mu = c(0, 2, -1.5), kappa = 4, p_mem = 0.6, p_nt = 0.2)
x <- seq(-pi, pi, length.out = 10000)
d <- dmixture3p(x, mu = c(0, 2, -1.5), kappa = 4, p_mem = 0.6, p_nt = 0.2)
hist(r, breaks = 60, freq = FALSE)
lines(x, d, type = "l", col = "red")

```
