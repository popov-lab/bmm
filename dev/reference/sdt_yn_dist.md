# Distribution functions for Yes/No SDT

Density and random generation for the yes/no signal detection theory
model, where the response is the number of "old"/"signal" responses out
of a fixed number of trials (a binomial likelihood).

## Usage

``` r
dsdt_yn(
  n_old,
  n_trials,
  stimulus,
  d,
  criterion,
  sdratio = 1,
  dist = c("normal", "gumbel_min", "gumbel_max", "logistic"),
  log = FALSE
)

rsdt_yn(
  n,
  n_trials,
  stimulus,
  d,
  criterion,
  sdratio = 1,
  dist = c("normal", "gumbel_min", "gumbel_max", "logistic")
)
```

## Arguments

- n_old:

  Integer vector. Number of "old"/"signal" responses.

- n_trials:

  Integer vector. Total number of trials per cell.

- stimulus:

  Numeric or logical vector (0/1). Stimulus type: 0 = noise, 1 = signal.

- d:

  Numeric. Sensitivity: \\d'\\ when `sdratio` is 1, and otherwise the
  balanced index \\d_a\\, the separation between the two distributions
  divided by the root-mean-square of their SDs (see
  [`sdt_yn()`](https://venpopov.com/bmm/dev/reference/sdt_yn.md)). The
  separation in noise units is `d * sqrt((1 + sdratio^2) / 2)`.

- criterion:

  Numeric. Response bias (decision boundary location), on the
  noise-standardized axis, i.e. in units of the noise distribution's SD.

- sdratio:

  Numeric. Ratio of signal to noise standard deviations (default 1,
  i.e., equal variance). Must be positive, and is on the **natural**
  scale — see the section below before reusing a fitted value.

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

  Integer. Number of observations to generate. `n_trials`, `stimulus`,
  and the model parameters are recycled to this length.

## Value

`dsdt_yn` returns the (log-)density (binomial probability). `rsdt_yn`
returns an integer vector with the number of "old"/"signal" responses
per observation.

## Parameter scales

As everywhere in `bmm`, these functions take their arguments on the
**natural** scale, while the parameters
[`sdt_yn()`](https://venpopov.com/bmm/dev/reference/sdt_yn.md)
*estimates* are on their link scale. `d` and `criterion` have identity
links and carry across unchanged, but `sdratio` has a log link: a fitted
`sdratio` of 0.375 is a ratio of `exp(0.375) = 1.455`, and passing
`0.375` here instead asks for a signal distribution 2.7 times narrower
than the noise. That is a legal value and raises no error, so
exponentiate first.

## References

Green, D. M., & Swets, J. A. (1966). *Signal detection theory and
psychophysics*. Wiley.

## Examples

``` r
# Density of yes/no SDT data
dsdt_yn(n_old = 80, n_trials = 100, stimulus = 1,
        d = 1.5, criterion = 0.2)
#> [1] 0.01136708

# Vectorized over observations
dsdt_yn(n_old = c(30, 80), n_trials = c(100, 100),
        stimulus = c(0, 1), d = 1.5, criterion = 0.2,
        log = TRUE)
#> [1] -7.463009 -4.477034

# Unequal variance from a fitted model: sdt_yn() reports sdratio on its log
# link, so exponentiate before passing it here
dsdt_yn(n_old = 80, n_trials = 100, stimulus = 1,
        d = 1.5, criterion = 0.2, sdratio = exp(0.375))
#> [1] 0.005538898
# Generate yes/no SDT data for a design
dat <- expand.grid(id = 1:20, stimulus = c(0L, 1L))
dat$n_trials <- 100L
dat$n_old <- rsdt_yn(nrow(dat), dat$n_trials, dat$stimulus,
                     d = 1.5, criterion = 0.2)
head(dat)
#>   id stimulus n_trials n_old
#> 1  1        0      100    15
#> 2  2        0      100    14
#> 3  3        0      100    25
#> 4  4        0      100    19
#> 5  5        0      100    14
#> 6  6        0      100    21
```
