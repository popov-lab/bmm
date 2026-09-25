# Yes/No Signal Detection Theory Model

Estimates sensitivity (`d`, which is \\d'\\ under equal variance and
\\d_a\\ once `sdratio` is estimated) and response bias (`criterion`)
from yes/no detection or old/new recognition counts.

## Usage

``` r
sdt_yn(
  response,
  stimulus,
  n_trials,
  dist = c("normal", "gumbel_min", "gumbel_max", "logistic"),
  links = NULL,
  ...
)
```

## Arguments

- response:

  The name of the variable in the dataset containing the count of
  "old"/"signal" responses for each cell.

- stimulus:

  The name of the variable in the dataset coding the stimulus type: 0
  (noise/new) and 1 (signal/old). Logical, and factor or character
  columns holding "0"/"1", are coerced automatically; anything else
  (e.g. "noise"/"signal" labels) must be recoded by hand, since bmm
  cannot guess which level is the signal.
  [`dsdt_yn()`](https://venpopov.com/bmm/dev/reference/sdt_yn_dist.md)
  and
  [`rsdt_yn()`](https://venpopov.com/bmm/dev/reference/sdt_yn_dist.md)
  take the same column but are stricter, accepting only numeric or
  logical input.

- n_trials:

  The name of the variable in the dataset containing the total number of
  trials for each cell. It may differ from cell to cell. Trial-level
  data also works: keep one row per trial, with the column named by
  `response` holding 0 or 1 and the column named by `n_trials` holding
  `1`s.

- dist:

  The noise distribution assumed for the latent evidence variable, given
  here by its cumulative distribution function. One of:

  - "normal" (default): Gaussian SDT, \\\Phi(x)\\

  - "gumbel_min": smallest-extreme-value SDT, \\1 - \exp(-\exp(x))\\
    (complementary log-log)

  - "gumbel_max": largest-extreme-value SDT, \\\exp(-\exp(-x))\\
    (log-log, as in
    [`evd::pgumbel`](https://rdrr.io/pkg/evd/man/gumbel.html))

  - "logistic": logistic SDT, \\1 / (1 + \exp(-x))\\

- links:

  A named list of link functions for the parameters, one entry per
  parameter you want to change, e.g. `links = list(d = "log")`. Only `d`
  and `criterion` can be set. `sdratio` keeps its log link, because the
  model's default of an equal-variance SD ratio is stored as the 0 that
  the log link maps to 1; read on any other link that same 0 is a
  different ratio, and on an identity link it is a ratio of zero that
  makes the signal trials' evidence infinite.

- ...:

  used internally for testing, ignore it

## Value

An object of class `bmmodel`

## Details

- **Domain:** Perception & Recognition Memory

- **Task:** Yes/No Detection or Old/New Recognition

- **Name:** Signal Detection Theory (Yes/No)

- **Citation:**

  - Green, D. M., & Swets, J. A. (1966). Signal detection theory and
    psychophysics. Wiley.

- **Requirements:**

  Provide pre-aggregated data with the following columns:

&nbsp;

- Response counts (response): number of 'old'/'signal' responses

- Stimulus type (stimulus): 0 = noise, 1 = signal

- Number of trials (n_trials): total trials per cell

&nbsp;

- **Parameters:**

  - `d`: d_a sensitivity (= d' when sdratio is fixed): the distance
    between the signal and noise distributions in units of their
    root-mean-square SD

  - `criterion`: Response bias: location of decision boundary

  - `sdratio`: SD ratio signal/noise, log link (0 = equal SDs):
    summary() prints the log of this ratio, not the ratio itself, so
    exp() a posterior value to read it

- **Fixed parameters:**

  - `sdratio` = 0

- **Default parameter links:**

  - d = identity; criterion = identity; sdratio = log

- **Default priors:**

  - `d`:

    - `main`: normal(1, 1)

    - `effects`: normal(0, 0.5)

    - `sd`: exponential(1)

  - `criterion`:

    - `main`: normal(0, 1.5)

    - `effects`: normal(0, 0.5)

    - `sd`: exponential(2)

  - `sdratio`:

    - `main`: normal(0, 0.5)

    - `effects`: normal(0, 0.3)

    - `sd`: exponential(2)

## Which sensitivity measure `d` is

**`d` is the familiar \\d'\\** whenever `sdratio` keeps its default of 0
(an SD ratio of 1, equal variance), which is every fit that does not
give `sdratio` a formula. The rest of this section only matters once you
estimate `sdratio`.

When the signal and noise distributions have different widths, their
separation only becomes dimensionless after choosing an SD to divide by.
`bmm` then reports \\d_a\\, the separation divided by the
root-mean-square of the two SDs: \$\$d_a = \sqrt{2}\\\delta / \sqrt{1 +
r^2},\$\$ where \\\delta\\ is the separation in noise-SD units and \\r\\
is the SD ratio, `exp(sdratio)`. This weights the two distributions
equally, and it is the measure Simpson and Fitter (1973), Macmillan and
Creelman (2005), and Mickes et al. (2007) recommend under unequal
variance. The classical noise-standardized index is \\d_N = \delta = d_a
\sqrt{(1 + r^2)/2}\\, so a published \\d'\\ from an unequal-variance
analysis is larger than `d` when \\r \> 1\\: by 13% at \\r = 1.25\\ and
by 33% at \\r = 1.6\\.

`d` is \\d_a\\ rather than \\d_N\\ because only \\d_a\\ is comparable
across conditions or subjects that differ in `sdratio`: two conditions
that are equally discriminable can show a large, confidently estimated
difference in \\d_N\\. When `sdratio` is estimated but constant across
the conditions you compare, the two indices differ by one common factor
and give the same contrasts up to scale.

**Units of the other parameters.** `criterion` is *not* rescaled. It is
the location of the decision boundary relative to the midpoint between
the two distributions, in noise-SD units, so under unequal variance `d`
(in root-mean-square SD units) and `criterion` (in noise-SD units) are
on different scales, and a ratio such as `criterion / d` mixes them.

**Extreme-value distributions.** For `dist = "normal"`, \\d_a\\ is also
the AUC-equivalent index, \\d_a = \sqrt{2}\\\Phi^{-1}(\mathrm{AUC})\\,
so it carries the same information as 2AFC accuracy. For `"gumbel_min"`
and `"gumbel_max"` that identity holds only under equal variance. With
`sdratio` estimated, \\d_a\\ keeps its balanced geometry but drifts away
from the AUC-equivalent index as the SD ratio moves away from 1 (for
`"gumbel_min"` at \\\delta = 1.5\\ and \\r = 2\\, \\d_a\\ is 31%
larger), so compare such fits on the AUC rather than on `d`.

Because `d` is a short name, a column called `d` in your data that is
also used as a predictor will collide with this parameter;
[`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md) warns when that
happens.

## Identifying `sdratio`

`sdratio` needs a design that supplies more than one operating point. A
single (hit, false-alarm) pair is two numbers for three unknowns, so
when every parameter is intercept-only with no random effects,
`sdratio ~ 1` returns its prior and `d` is pulled along the resulting
ridge: sampling converges, `Rhat` is fine, and the profile likelihood
over `sdratio` is flat to 1e-12.
[`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md) warns for that
one design, which is the only shape that provably cannot work.

Any linear predictor that moves the operating point along the ROC
supplies what is missing, and it need not sit on `criterion`: a
sensitivity manipulation (`d ~ 0 + condition`, a study-time or strength
manipulation with bias held constant) identifies `sdratio` just as a
criterion manipulation does, and so does between-subject variation
entering through a random effect such as `criterion ~ 1 + (1 | id)` —
though, like any predictor, only in proportion to how far it actually
moves the operating point: a random effect with little between-subject
spread carries little information and will not trigger the warning
above, because that warning counts formula terms, not how much they move
the design.

A criterion manipulation is still the cleanest design, because it traces
the ROC at fixed sensitivity: give `criterion` a predictor that shifts
the decision boundary — a base-rate, payoff, or confidence manipulation
— as in `criterion ~ 0 + condition`; see
[broeder_schuetz_2009_e3](https://venpopov.com/bmm/dev/reference/broeder_schuetz_2009_e3.md).
Leaving `sdratio` at its default is always identified.

## Reading `sdratio` and carrying it to [`dsdt_yn()`](https://venpopov.com/bmm/dev/reference/sdt_yn_dist.md)/[`rsdt_yn()`](https://venpopov.com/bmm/dev/reference/sdt_yn_dist.md)

As in every `bmm` model, the parameters the model *estimates* are on
their link scale, while the distribution functions take their arguments
on the *natural* scale. `sdratio` has a log link, so
[`summary()`](https://rdrr.io/r/base/summary.html) reports \\\log r\\
whereas
[`dsdt_yn()`](https://venpopov.com/bmm/dev/reference/sdt_yn_dist.md) and
[`rsdt_yn()`](https://venpopov.com/bmm/dev/reference/sdt_yn_dist.md)
expect the ratio \\r\\ itself (their default is 1, equal variance).
Exponentiate before carrying a posterior value across:

    r <- exp(as_draws_matrix(fit)[, "b_sdratio_Intercept"])

A posterior mean of `sdratio = 0.22` is a ratio of `exp(0.22) = 1.25`.
Passing `0.22` straight to
[`rsdt_yn()`](https://venpopov.com/bmm/dev/reference/sdt_yn_dist.md)
instead asks for a signal distribution 4.5 times *narrower* than the
noise — a legal value that raises no error, and the one mistake worth
checking for in a posterior predictive check written by hand. `d` and
`criterion` have identity links, so they carry across unchanged.

The **zROC slope** reported in the recognition-memory literature is the
reciprocal of that ratio, `1 / exp(sdratio)`, so a `sdratio` posterior
mean of 0.375 is a zROC slope of 0.69.

The same log link applies going the other way: a constant you supply
yourself, whether as `bmf(sdratio = )` or through a hand-written
`brms::set_prior(..., dpar = "sdratio")`, is read on it too.
`bmf(sdratio = 1)` does not fix a ratio of 1 — it fixes `exp(1) = 2.72`,
and [`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md) raises no
warning; a fixed ratio of 1.25 needs `bmf(sdratio = log(1.25))`. Both
`sdratio` defaults
[`default_prior()`](https://venpopov.com/bmm/dev/reference/default_prior.bmmformula.md)
reports — `normal(0, 0.5)` on the intercept and `exponential(2)` on the
random-effect SDs — are on that same scale, unannotated.

## Terms used on this page

- **noise-standardized axis** — latent evidence expressed in units of
  the noise distribution's SD. `criterion` always lives on this axis.

- \\\delta\\ (**separation**) — the distance between the signal and
  noise means, on that axis.

- \\d_N\\ — the classical \\d'\\: \\\delta\\ itself, i.e. the separation
  in noise-SD units.

- \\d_a\\ — the separation divided by the root-mean-square of the two
  SDs, which is what `d` reports. \\d_a = d_N\\ under equal variance.

- **operating point** — one (false-alarm rate, hit rate) pair, i.e. one
  point of an ROC curve. One condition gives one operating point.

- **AUC** — the area under that ROC curve, equivalently the probability
  that a random signal trial yields more evidence than a random noise
  trial. Obtain it from the posterior with `pnorm(d / sqrt(2))` for
  `dist = "normal"`.

- `main` / `effects` / `sd` — the keys of the default priors shown in
  the model description above: `main` is the prior on the intercept,
  `effects` the prior on regression coefficients, and `sd` the prior on
  the standard deviations of the parameter's random effects.

## References

Green, D. M., & Swets, J. A. (1966). *Signal detection theory and
psychophysics*. Wiley.

DeCarlo, L. T. (1998). Signal detection theory and generalized linear
models. *Psychological Methods*, *3*(2), 186–205.
[doi:10.1037/1082-989X.3.2.186](https://doi.org/10.1037/1082-989X.3.2.186)

Simpson, A. J., & Fitter, M. J. (1973). What is the best index of
detectability? *Psychological Bulletin*, *80*(6), 481–488.
[doi:10.1037/h0035203](https://doi.org/10.1037/h0035203)

Macmillan, N. A., & Creelman, C. D. (2005). *Detection theory: A user's
guide* (2nd ed.). Erlbaum.

Mickes, L., Wixted, J. T., & Wais, P. E. (2007). A direct test of the
unequal-variance signal detection model of recognition memory.
*Psychonomic Bulletin & Review*, *14*(5), 858–865.
[doi:10.3758/BF03194112](https://doi.org/10.3758/BF03194112)

## See also

[`sdt_d()`](https://venpopov.com/bmm/dev/reference/SDTdist.md) and
[`sdt_criterion()`](https://venpopov.com/bmm/dev/reference/SDTdist.md)
compute the `d` and `criterion` of this model in closed form from a
single pair of observed hit and false-alarm rates, without fitting: use
them for a quick check of a fitted value, and this model when you need a
hierarchical or condition-wise estimate, or unequal variance.

## Examples

``` r
if (FALSE) { # \dontrun{
dat <- expand.grid(id = 1:20, stimulus = c(0L, 1L))
dat$n_trials <- 100L
dat$n_old <- rsdt_yn(nrow(dat), dat$n_trials, dat$stimulus,
                     d = 1.5, criterion = 0.2)

model <- sdt_yn(
  response = "n_old",
  stimulus = "stimulus",
  n_trials = "n_trials"
)

fit <- bmm(
  formula = bmf(d ~ 1, criterion ~ 1),
  data = dat,
  model = model,
  cores = 4,
  backend = "cmdstanr"
)

# Sensitivity and bias per participant
fit_re <- bmm(
  formula = bmf(d ~ 1 + (1 | id), criterion ~ 1 + (1 | id)),
  data = dat,
  model = model,
  cores = 4,
  backend = "cmdstanr"
)

# Unequal-variance yes/no SDT. sdratio needs more than one operating point:
# on the single-condition `dat` above it would not be identified.
# `model` already names the columns this dataset uses.
fit_uv <- bmm(
  formula = bmf(d ~ 1, criterion ~ 0 + condition, sdratio ~ 1),
  data = broeder_schuetz_2009_e3,
  model = model,
  cores = 4,
  backend = "cmdstanr"
)

# Simulating from a fitted unequal-variance model: the model reports
# sdratio on its log link, rsdt_yn() takes the ratio itself.
sdratio_posterior <- 0.375
rsdt_yn(2, 100L, c(0L, 1L),
        d = 1.3, criterion = 0.1, sdratio = exp(sdratio_posterior))
} # }
```
