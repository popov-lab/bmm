# m-Alternative Forced Choice Signal Detection Theory Model

m-Alternative Forced Choice Signal Detection Theory Model

## Usage

``` r
sdt_mafc(
  response,
  n_trials,
  m,
  dist = c("normal", "gumbel_min", "gumbel_max", "logistic"),
  links = NULL,
  ...
)
```

## Arguments

- response:

  A single string naming the column with counts of correct responses.

- n_trials:

  The name of the variable containing the total number of trials per
  cell.

- m:

  Either a single integer \>= 2 giving the number of alternatives
  (constant across all rows), or a single string naming a data column
  that gives the number of alternatives per row. A column lets trials
  with different set sizes be fit jointly.

- dist:

  The distribution assumed for the latent evidence, given here by its
  cumulative distribution function. One of:

  - "normal" (default): Gaussian m-AFC, \\\Phi(x)\\

  - "gumbel_min": smallest-extreme-value m-AFC, \\1 - \exp(-\exp(x))\\
    (complementary log-log)

  - "gumbel_max": largest-extreme-value m-AFC, \\\exp(-\exp(-x))\\
    (log-log, as in
    [`evd::pgumbel`](https://rdrr.io/pkg/evd/man/gumbel.html))

  - "logistic": logistic m-AFC, \\1 / (1 + \exp(-x))\\

- links:

  A named list of link functions for the parameters.

- ...:

  used internally for testing, ignore it

## Value

An object of class `bmmodel`

## Details

- **Domain:** Perception & Recognition Memory

- **Task:** m-Alternative Forced Choice

- **Name:** Signal Detection Theory (m-AFC)

- **Citation:**

  - Green, D. M., & Swets, J. A. (1966). Signal detection theory and
    psychophysics. Wiley.

  - DeCarlo, L. T. (2012). On a signal detection approach to
    m-alternative forced choice with bias, with maximum likelihood and
    Bayesian approaches to estimation. Journal of Mathematical
    Psychology, 56(3), 196-207.
    https://doi.org/10.1016/j.jmp.2012.02.004

- **Requirements:**

  Provide pre-aggregated accuracy data with the following columns:

&nbsp;

- Response counts (response): number of correct responses

- Number of trials (n_trials): total trials per cell No stimulus column
  needed (each trial has exactly one signal alternative)

&nbsp;

- **Parameters:**

  - `d`: Sensitivity: d', the distance between the signal and distractor
    distributions in SD units (m-AFC assumes a common SD, so this is
    also the d_a that sdt_yn reports under unequal variance)

- **Fixed parameters:**

- **Default parameter links:**

  - d = identity

- **Default priors:**

  - `d`:

    - `main`: normal(1, 1)

    - `effects`: normal(0, 0.5)

    - `sd`: exponential(1)

Models accuracy in m-AFC tasks where each trial presents one signal
among `m` alternatives and the observer chooses the strongest one. Only
`d` is estimated (m-AFC has no response bias). The probability correct,
\\P_c = \int f(x - d)\\ F(x)^{m-1}\\ dx\\, is computed per noise
distribution: a closed-form softmax for `gumbel_max`, a closed-form
Gamma ratio for `gumbel_min`, Gauss-Hermite quadrature for `normal`
(\\\Phi(d/\sqrt{2})\\ at m = 2), and Gauss-Legendre quadrature for
`logistic`.

## Sensitivity is on the same scale as [`sdt_yn()`](https://popov-lab.github.io/bmm/dev/reference/sdt_yn.md)

`d` is \\d'\\. m-AFC assumes the signal and distractor distributions
share an SD, so there is no `sdratio` parameter, and \\d'\\ coincides
with the balanced index \\d_a\\ that
[`sdt_yn()`](https://popov-lab.github.io/bmm/dev/reference/sdt_yn.md)
reports once its `sdratio` is estimated.

At `m = 2` the two models agree: 2AFC proportion correct equals the area
under the yes/no ROC (Green's theorem), and for Gaussian noise both give
\\P_c = \Phi(d/\sqrt{2})\\. The same observer therefore yields the same
`d` whether it is measured by a yes/no ROC or by 2AFC accuracy, which is
the property that makes \\d_a\\ the right common scale for the SDT
family. For `dist = "normal"` this holds even when
[`sdt_yn()`](https://popov-lab.github.io/bmm/dev/reference/sdt_yn.md)
estimates unequal variance. For the other distributions it holds when
[`sdt_yn()`](https://popov-lab.github.io/bmm/dev/reference/sdt_yn.md)
keeps `sdratio` at its default, because only the Gaussian \\d_a\\ is
exactly the AUC-equivalent index (see the sensitivity section of
[`sdt_yn()`](https://popov-lab.github.io/bmm/dev/reference/sdt_yn.md)).

## References

Green, D. M., & Swets, J. A. (1966). *Signal detection theory and
psychophysics*. Wiley.

DeCarlo, L. T. (2012). On a signal detection approach to m-alternative
forced choice with bias, with maximum likelihood and Bayesian approaches
to estimation. *Journal of Mathematical Psychology*, *56*(3), 196–207.
[doi:10.1016/j.jmp.2012.02.004](https://doi.org/10.1016/j.jmp.2012.02.004)

## Examples

``` r
if (FALSE) { # \dontrun{
dat <- data.frame(id = 1:20, n_trials = 200L)
dat$n_correct <- rsdt_mafc(nrow(dat), dat$n_trials, m = 4,
                           d = rnorm(20, 1.5, 0.4))

model <- sdt_mafc(
  response = "n_correct",
  n_trials = "n_trials",
  m = 4
)

fit <- bmm(
  formula = bmf(d ~ 1),
  data = dat,
  model = model,
  cores = 4,
  backend = "cmdstanr"
)
} # }
```
