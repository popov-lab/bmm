# Deprecated: rejection sampling

`rejection_sampling()` is deprecated since bmm 1.5.0 and will be removed
in bmm 1.6.0 at the earliest. It was the sampler behind
[`rsdm()`](https://venpopov.com/bmm/dev/reference/SDMdist.md), which no
longer needs the export; nothing replaces it. It still performs
rejection sampling from a target distribution, with a warning. Each draw
can come from its own target: draw `i` is sampled from `f` evaluated at
the `i`-th element of every per-draw argument in `...`, under the
envelope `max_f[i]`.

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
samples <- rejection_sampling(1000, target_density, max_f = target_density(0), proposal)
#> Warning: `rejection_sampling()` is deprecated since bmm 1.5.0 and will be removed in bmm 1.6.0 at the earliest.
```
