# Psychometric Function Model

Estimates the psychometric function, the probability of a positive
response as a function of stimulus intensity, with its midpoint, width,
guess rate and lapse rate, from trial-level binary responses or from
response counts per stimulus level.

## Usage

``` r
psychometric_function(
  response,
  intensity,
  n_trials = NULL,
  sigmoid = c("normal", "logistic", "gumbel_min", "gumbel_max", "weibull", "lognormal",
    "loglogistic"),
  guess = NULL,
  lapse = NULL,
  ...
)
```

## Arguments

- response:

  The name of the variable in the data with the responses: 0/1 (or
  logical) per trial when `n_trials` is `NULL`, and otherwise the number
  of positive responses out of `n_trials`.

- intensity:

  The name of the variable in the data with the stimulus intensity of
  each row.

- n_trials:

  The name of the variable in the data with the number of trials per
  row, for data aggregated into counts. `NULL` (default) reads
  `response` as one binary response per row.

- sigmoid:

  The sigmoid \\F\\: one of `"normal"` (default), `"logistic"`,
  `"gumbel_min"`, `"gumbel_max"`, `"weibull"`, `"lognormal"` or
  `"loglogistic"`. See the Sigmoids section.

- guess:

  The guess rate. `NULL` (default) estimates it; a number in \[0, 1)
  fixes it; the name of a variable in the data fixes it per row.

- lapse:

  The lapse rate. `NULL` (default) estimates it; a number in \[0, 1)
  fixes it; the name of a variable in the data fixes it per row.

- ...:

  used internally for testing, ignore it

## Value

An object of class `bmmodel`

## Details

- **Domain:** Perception

- **Task:** Detection, recognition and confidence judgments

- **Name:** Psychometric Function

- **Citation:**

  - Wichmann, F. A., & Hill, N. J. (2001). The psychometric function: I.
    Fitting, sampling, and goodness of fit. Perception & Psychophysics,
    63(8), 1293-1313. https://doi.org/10.3758/BF03194544

  - Schütt, H. H., Harmeling, S., Macke, J. H., & Wichmann, F. A.
    (2016). Painfree and accurate Bayesian estimation of psychometric
    functions for (potentially) overdispersed data. Vision Research,
    122, 105-123. https://doi.org/10.1016/j.visres.2016.02.002

- **Requirements:**

  Provide data with one row per trial or per stimulus level:

&nbsp;

- Response (response): 0/1 per trial, or the number of positive
  responses out of n_trials

- Stimulus intensity (intensity): the level of the manipulated stimulus
  dimension

- Number of trials (n_trials, optional): trials per row for counts

&nbsp;

- **Parameters:**

  - `midpoint`: Stimulus intensity at which the sigmoid F is 0.5,
    halfway between the asymptotes, in the units of the intensity

  - `width`: Width of the sigmoid: the distance on the intensity axis
    from F = 0.05 to F = 0.95; log link

  - `guess`: Guess rate: the lower asymptote, the response probability
    far below the midpoint; logit link

  - `lapse`: Lapse rate: the proportion of trials on which the observer
    does not attend and responds at the guess rate; the upper asymptote
    is 1 - lapse \* (1 - guess); logit link

- **Fixed parameters:**

- **Default parameter links:**

  - midpoint = identity; width = log; guess = logit; lapse = logit

- **Default priors:**

  - `midpoint`:

    - `main`: normal(0.5, 0.5)

    - `effects`: normal(0, 0.25)

    - `sd`: exponential(4)

  - `width`:

    - `main`: normal(-0.693, 1)

    - `effects`: normal(0, 0.5)

    - `sd`: exponential(2)

  - `guess`:

    - `main`: normal(-3, 1.3)

    - `effects`: normal(0, 1)

    - `sd`: exponential(1)

  - `lapse`:

    - `main`: normal(-3, 1.3)

    - `effects`: normal(0, 1)

    - `sd`: exponential(1)

The default links listed above are those of the default sigmoid;
`midpoint` takes a log link for `"weibull"`, `"lognormal"` and
`"loglogistic"`.

The probability of a positive response ("yes", "correct", "longer") at
intensity \\x\\ is \$\$\psi(x) = \gamma + (1 - \gamma)(1 - \lambda)
F(x),\$\$ where \\\gamma\\ is `guess`, \\\lambda\\ is `lapse` and \\F\\
the sigmoid. A lapse is a trial on which the observer does not attend
and so responds at the guess rate, which keeps \\\psi\\ increasing for
every value of the two rates.

## Midpoint and width

The sigmoid is placed by `midpoint`, the intensity at which \\F = 0.5\\,
and scaled by `width`, the distance on the intensity axis between \\F =
0.05\\ and \\F = 0.95\\ (Schütt et al., 2016). Both mean the same thing
for every sigmoid, so estimates are comparable across sigmoids.

`midpoint` is not the threshold at a target performance level: because
of the guess and lapse rates, \\F = 0.5\\ is \\\psi = 0.75\\ in 2AFC
without lapses, but \\\psi = 0.5\\ in a yes/no task. Use
[`psychometric_threshold()`](https://venpopov.com/bmm/dev/reference/psychometric_threshold.md)
for the intensity at any performance level, with its posterior, and for
the slope there.

## Sigmoids

- `"normal"`: cumulative normal (probit), on the intensity axis

- `"logistic"`: logistic, on the intensity axis

- `"gumbel_min"`: \\1 - \exp(-\exp(z))\\, skewed towards low intensities
  (psignifit's `"gumbel"`)

- `"gumbel_max"`: \\\exp(-\exp(-z))\\, skewed towards high intensities
  (psignifit's `"rgumbel"`)

- `"weibull"`: `"gumbel_min"` on log intensity, the Weibull function

- `"lognormal"`, `"loglogistic"`: `"normal"` and `"logistic"` on log
  intensity

The last three read the intensity on the log scale, so it must be
positive. Their `midpoint` is still reported in the units of the
intensity, through a log link, while their `width` is in log units.

## Guess and lapse rates

Each rate is either estimated, with a logit link and its own formula, or
fixed by the design. Fix them in the call, on the probability scale, not
in the formula: `guess = 0.5` for 2AFC, `guess = 1 / m` for m-AFC,
`guess = "chance"` to read a per-row guess rate from a column (for
designs that mix 2AFC and 4AFC trials), `lapse = 0` for no lapses. A
fixed value of zero is allowed. Leave `guess` estimated for yes/no
tasks, where it is the false-alarm rate far below the midpoint.

The lapse rate of Wichmann and Hill (2001), and of psignifit and
Palamedes, is the gap between the upper asymptote and 1, which is
`lapse * (1 - guess)` here: in 2AFC, half the `lapse` this model
reports.
[`psychometric_threshold()`](https://venpopov.com/bmm/dev/reference/psychometric_threshold.md)
returns both.

## Default priors

Stimulus units are arbitrary, so the priors on `midpoint` and `width`
are scaled to the range of the intensities in the data (on the log scale
for the log-scale sigmoids), as psignifit does. With \\c\\ the centre
and \\R\\ the range, the `midpoint` intercept gets `normal(c, R / 2)`,
its effects `normal(0, R / 4)` and its random-effect SDs
`exponential(4 / R)`; the log of `width` gets `normal(log(R / 2), 1)`.
The priors listed above are these for a range of 0 to 1;
[`default_prior()`](https://venpopov.com/bmm/dev/reference/default_prior.bmmformula.md)
shows the ones for your data.
[`update()`](https://rdrr.io/r/stats/update.html) with `newdata` alone
keeps the priors of the original fit, and says so when the new
intensities shift the centre by more than a quarter of the old range or
change the range by more than a factor of two. Passing the formula again
as `formula.` scales the priors to the new data; `prior` sets your own.

## References

Schütt, H. H., Harmeling, S., Macke, J. H., & Wichmann, F. A. (2016).
Painfree and accurate Bayesian estimation of psychometric functions for
(potentially) overdispersed data. *Vision Research*, *122*, 105–123.
[doi:10.1016/j.visres.2016.02.002](https://doi.org/10.1016/j.visres.2016.02.002)

Wichmann, F. A., & Hill, N. J. (2001). The psychometric function: I.
Fitting, sampling, and goodness of fit. *Perception & Psychophysics*,
*63*(8), 1293–1313.
[doi:10.3758/BF03194544](https://doi.org/10.3758/BF03194544)

## See also

[`dpsychometric_function()`](https://venpopov.com/bmm/dev/reference/psychometric_function_dist.md),
[`rpsychometric_function()`](https://venpopov.com/bmm/dev/reference/psychometric_function_dist.md),
[`psychometric_threshold()`](https://venpopov.com/bmm/dev/reference/psychometric_threshold.md),
and [`sdt_mafc()`](https://venpopov.com/bmm/dev/reference/sdt_mafc.md)
for m-AFC accuracy at a single stimulus level.

## Examples

``` r
if (FALSE) { # \dontrun{
# 2AFC contrast detection, 40 trials per level, 10 observers
dat <- expand.grid(contrast = c(0.01, 0.02, 0.04, 0.08, 0.16), id = 1:10)
dat$n_trials <- 40L
dat$n_correct <- rpsychometric_function(nrow(dat), dat$contrast, midpoint = 0.04,
                                        width = 2, guess = 0.5, lapse = 0.02,
                                        n_trials = dat$n_trials, sigmoid = "weibull")

model <- psychometric_function(response = "n_correct", intensity = "contrast",
                               n_trials = "n_trials", sigmoid = "weibull",
                               guess = 0.5)

fit <- bmm(
  formula = bmf(midpoint ~ 1 + (1 | id), width ~ 1, lapse ~ 1),
  data = dat,
  model = model,
  cores = 4,
  backend = "cmdstanr"
)

# contrast at 75% correct, and the slope there
psychometric_threshold(fit, p = 0.75)
} # }
```
