# Rejection Sampling

Performs rejection sampling to generate samples from a target
distribution. Each draw can come from its own target: draw `i` is
sampled from `f` evaluated at the `i`-th element of every per-draw
argument in `...`, under the envelope `max_f[i]`.

## Usage

``` r
rejection_sampling(n, f, max_f, proposal_fun, ...)
```

## Arguments

- n:

  Integer. The number of samples to generate.

- f:

  Function. The target density divided by the proposal density, up to a
  constant; with a uniform proposal, the target density itself. Its
  first argument takes a vector of proposals, generally not of length
  `n`; `f` must be vectorized over it and over the per-draw arguments in
  `...`.

- max_f:

  Numeric. A finite upper bound of `f`, either a single value or one
  value per draw (length `n`). A bound below the maximum of `f` biases
  the draws without a warning.

- proposal_fun:

  Function. A function that generates samples from the proposal
  distribution.

- ...:

  Additional arguments to be passed to the target density function `f`.
  With `n > 1`, arguments of length `n` are taken per draw, so draw `i`
  uses their `i`-th elements. Arguments of any other length are passed
  whole to every call of `f`; recycle them with `rep_len(x, n)` to use
  them per draw. Pass constants whose length may equal `n`, such as a
  lookup table, through the closure of `f` instead of `...`.

## Value

A numeric vector of length `n` containing samples from the target
distribution.

## Examples

``` r
target_density <- function(x) brms::dvon_mises(x, mu = 0, kappa = 10)
proposal <- function(n) runif(n, min = -pi, max = pi)
samples <- rejection_sampling(10000, target_density, max_f = target_density(0), proposal)
hist(samples, freq = FALSE)
curve(target_density, col = "red", add = TRUE)


# one location per draw
mu <- rep(c(0, 2), 5000)
samples <- rejection_sampling(
  10000, brms::dvon_mises, max_f = brms::dvon_mises(0, 0, 10), proposal,
  mu = mu, kappa = 10
)
tapply(samples, mu, mean)
#>             0             2 
#> -0.0009627108  1.9951207670 
```
