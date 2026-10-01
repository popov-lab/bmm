# Distribution functions for m-AFC SDT

Density and random generation for m-alternative forced choice signal
detection theory (DeCarlo, 2012). Models accuracy in tasks where one of
`m` alternatives contains the signal. Only the `d` parameter is
estimated (no criterion). All arguments are recycled to the length of
the longest one, so passing vectors of `d`, `m`, or `n_trials` generates
(or evaluates) one observation per element.

## Usage

``` r
dsdt_mafc(
  n_correct,
  n_trials,
  m,
  d,
  dist = c("normal", "gumbel_min", "gumbel_max", "logistic"),
  log = FALSE
)

rsdt_mafc(
  n,
  n_trials,
  m,
  d,
  dist = c("normal", "gumbel_min", "gumbel_max", "logistic")
)
```

## Arguments

- n_correct:

  Integer vector. Number of correct responses.

- n_trials:

  Integer vector. Total number of trials per observation.

- m:

  Integer vector. Number of alternatives per observation. Must be at
  least 2.

- d:

  Numeric vector. Sensitivity \\d'\\: the distance between the signal
  and distractor distributions in SD units. m-AFC assumes a common scale
  for the two distributions, so this is also the balanced index \\d_a\\
  that
  [`sdt_yn()`](https://popov-lab.github.io/bmm/dev/reference/sdt_yn.md)
  reports.

- dist:

  Character. The distribution assumed for the latent evidence, given
  here by its cumulative distribution function:

  - "normal" (default): Gaussian, \\\Phi(x)\\

  - "gumbel_min": smallest extreme value, \\1 - \exp(-\exp(x))\\ (the
    complementary log-log distribution)

  - "gumbel_max": largest extreme value, \\\exp(-\exp(-x))\\ (the
    log-log distribution, as in
    [`evd::pgumbel`](https://rdrr.io/pkg/evd/man/gumbel.html))

  - "logistic": \\1 / (1 + \exp(-x))\\

- log:

  Logical. If `TRUE`, returns log-density (default `FALSE`).

- n:

  Integer. Number of observations to generate. `n_trials`, `m`, and `d`
  are recycled to this length.

## Value

`dsdt_mafc` returns the (log-)density (binomial probability).
`rsdt_mafc` returns an integer vector with the number of correct
responses per observation.

## References

DeCarlo, L. T. (2012). On a signal detection approach to m-alternative
forced choice with bias, with maximum likelihood and Bayesian approaches
to estimation. *Journal of Mathematical Psychology*, *56*(3), 196–207.
[doi:10.1016/j.jmp.2012.02.004](https://doi.org/10.1016/j.jmp.2012.02.004)

## Examples

``` r
# 4-AFC density
dsdt_mafc(n_correct = 80, n_trials = 100, m = 4, d = 1.5)
#> [1] 0.008272764
# Generate 4-AFC data for 20 subjects with varying sensitivity
dat <- data.frame(id = 1:20, n_trials = 200L)
dat$n_correct <- rsdt_mafc(nrow(dat), dat$n_trials, m = 4,
                           d = rnorm(20, 1.5, 0.4))
head(dat)
#>   id n_trials n_correct
#> 1  1      200       146
#> 2  2      200       155
#> 3  3      200       159
#> 4  4      200       140
#> 5  5      200       137
#> 6  6      200       176
```
