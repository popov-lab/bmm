# Distribution functions for the two-parameter mixture model (mixture2p)

Density, distribution, and random generation functions for the
two-parameter mixture model with the location of `mu`, precision of
memory representations `kappa` and probability of recalling items from
memory `p_mem`.

## Usage

``` r
dmixture2p(x, mu = 0, kappa = 5, p_mem = 0.6, log = FALSE)

pmixture2p(q, mu = 0, kappa = 7, p_mem = 0.8)

qmixture2p(p, mu = 0, kappa = 5, p_mem = 0.6)

rmixture2p(n, mu = 0, kappa = 5, p_mem = 0.6)
```

## Arguments

- x:

  Vector of observed responses

- mu:

  Vector of locations

- kappa:

  Vector of precision values

- p_mem:

  Vector of probabilities for memory recall

- log:

  Logical; if `TRUE`, values are returned on the log scale.

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
representations in visual working memory. Nature, 453.

## Examples

``` r
# generate random samples from the mixture2p model and overlay the density
r <- rmixture2p(10000, mu = 0, kappa = 4, p_mem = 0.8)
x <- seq(-pi, pi, length.out = 10000)
d <- dmixture2p(x, mu = 0, kappa = 4, p_mem = 0.8)
hist(r, breaks = 60, freq = FALSE)
lines(x, d, type = "l", col = "red")

```
