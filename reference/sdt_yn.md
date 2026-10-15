# Yes/No Signal Detection Theory Model

Estimates sensitivity (`d`, which is \\d'\\ under equal variance and
\\d_a\\ once `sdratio` is estimated) and the decision criterion from
yes/no detection or old/new recognition counts. The criterion is
estimated either as `criterion`, measured from the mean of the noise
distribution, or as `bias`, measured from the midpoint between the noise
and signal means; the formula names the one to estimate.

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
  [`dsdt_yn()`](https://venpopov.com/bmm/reference/sdt_yn_dist.md) and
  [`rsdt_yn()`](https://venpopov.com/bmm/reference/sdt_yn_dist.md) take
  the same column but are stricter, accepting only numeric or logical
  input.

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
  parameter you want to change, e.g. `links = list(d = "log")`. Only
  `d`, `criterion` and `bias` can be set. `sdratio` keeps its log link,
  because the model's default of an equal-variance SD ratio is stored as
  the 0 that the log link maps to 1; read on any other link that same 0
  is a different ratio, and on an identity link it is a ratio of zero
  that makes the signal trials' evidence infinite.

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

  - `criterion`: Decision criterion, measured from the mean of the noise
    distribution in noise-SD units; the false-alarm rate depends on it
    alone

  - `bias`: Alternative to criterion: the decision criterion measured
    from the midpoint between the noise and signal means (0 = neutral,
    positive = conservative). Use one of the two

  - `sdratio`: SD ratio signal/noise, log link (0 = equal SDs):
    summary() prints the log of this ratio, not the ratio itself, so
    exp() a posterior value to read it

- **Fixed parameters:**

  - `sdratio` = 0

- **Default parameter links:**

  - d = identity; criterion = identity; bias = identity; sdratio = log

- **Default priors:**

  - `d`:

    - `main`: normal(1, 1)

    - `effects`: normal(0, 0.5)

    - `sd`: exponential(1)

  - `criterion`:

    - `main`: normal(0.5, 1.5)

    - `effects`: normal(0, 0.5)

    - `sd`: exponential(1)

  - `bias`:

    - `main`: normal(0, 1.5)

    - `effects`: normal(0, 0.5)

    - `sd`: exponential(2)

  - `sdratio`:

    - `main`: normal(0, 0.3)

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
[`sdt_sensitivity()`](https://venpopov.com/bmm/reference/sdt_sensitivity.md)
returns \\d_N\\ and the signal-standardized \\d_S\\ from a fitted model,
and [`summary()`](https://rdrr.io/r/base/summary.html) notes when `d` is
\\d_a\\.

`d` is \\d_a\\ rather than \\d_N\\ because only \\d_a\\ is comparable
across conditions or subjects that differ in `sdratio`: two conditions
that are equally discriminable can show a large, confidently estimated
difference in \\d_N\\. When `sdratio` is estimated but constant across
the conditions you compare, the two indices differ by one common factor
and give the same contrasts up to scale.

**Units of the other parameters.** `criterion` and `bias` are *not*
rescaled: both are in noise-SD units, so under unequal variance `d` (in
root-mean-square SD units) and the criterion are on different scales,
and a ratio such as `bias / d` mixes them.

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
[`bmm()`](https://venpopov.com/bmm/reference/bmm.md) warns when that
happens.

## Choosing between `criterion` and `bias`

Both parameters place the same decision boundary on the
noise-standardized evidence axis, and they differ only in the point they
measure it from:

- `criterion` is the distance from the mean of the noise distribution.
  The false-alarm rate depends on it alone: \\FA = 1 - F(criterion)\\.

- `bias` is the distance from the midpoint between the noise and signal
  means: 0 is the neutral point, positive values a conservative bias
  toward "new"/"noise". Under equal variance it is the index usually
  written \\c\\.

They are related by \\bias = criterion - \delta / 2\\, where \\\delta =
d \sqrt{(1 + r^2)/2}\\ is the separation in noise-SD units. Write one of
them in the formula; a formula that names neither estimates
`criterion ~ 1`, and one that names both is refused.
[`sdt_criterion()`](https://venpopov.com/bmm/reference/SDTdist.md) and
[`sdt_bias()`](https://venpopov.com/bmm/reference/SDTdist.md) compute
the two from a single hit and false-alarm rate.

**The choice matters whenever `d` or `sdratio` differs between rows that
share a criterion term**, because a linear formula holds constant the
parameter you write. `criterion ~ 1` with `d ~ 1 + strength` says the
observer keeps one boundary while targets get stronger; `bias ~ 1` with
the same `d` formula says the boundary moves up by half of every change
in separation. When the criterion and `d` share condition-wise (factor)
predictors and `sdratio` does not vary within those cells, or share
predictors while `d` keeps its identity link and `sdratio` is the same
in every row, the two give the same likelihood and their coefficients
differ by \\\delta / 2\\. Otherwise they are different models. Even with
the same likelihood, the posteriors differ through the default priors:
`criterion` takes `normal(0.5, 1.5)` for its intercept, or for each cell
under `~ 0 + condition`, and `exponential(1)` for its random-effect SDs,
while `bias` takes `normal(0, 1.5)` and `exponential(2)`.

Use `criterion` when the lures, and the observer's boundary, do not
depend on the manipulation of `d`: the classic recognition design in
which weak and strong targets are tested against one shared set of lures
under the same criteria. Its lure rows have no strength level, and with
`criterion` they need none, because the false-alarm rate does not
involve `d`. Use `bias` when you want to hold, or compare, the
observer's position relative to the two distributions, e.g. a bias
contrast between groups that differ in sensitivity.

`bias` is measured from the midpoint, which is the point where the two
densities are equal – the optimal boundary for equal base rates and
payoffs – only for the symmetric distributions (`"normal"`,
`"logistic"`) under equal variance. Under unequal variance, or for the
Gumbel distributions, `bias = 0` is the midpoint, not the optimum.

## Identifying `sdratio`

`sdratio` needs a design that supplies more than one operating point,
and whether it gets one depends on which criterion parameter the formula
holds constant.

**With `criterion`**, the boundary must move among rows that share the
same `d`: give `criterion` a predictor that shifts it while `d` stays
put, such as a base-rate, payoff, or confidence manipulation (`d ~ 1`,
`criterion ~ 0 + condition`; see
[broeder_schuetz_2009_e3](https://venpopov.com/bmm/reference/broeder_schuetz_2009_e3.md)).
A criterion that is the same in every row is fixed by the false-alarm
rate, and each hit rate then fixes only \\(\delta - criterion) / r\\, so
`sdratio` trades against `d`. A change in `sdratio` is matched by
shifting the separation of every cell by one common amount and rescaling
it by one common factor, and whenever the formula for `d` can make that
change, `sdratio` returns its prior. A factor can, under the identity
and the log link alike (`d ~ 0 + condition` or `d ~ 1 + condition`),
because it gives every cell a value of its own. An intercept with a
numeric predictor can under the identity link (`d ~ 1 + strength`),
because the intercept shifts and the slope rescales. `sdratio` is
identified only when the formula constrains `d` across cells beyond
that. A numeric predictor without an intercept (`d ~ 0 + strength`)
holds `d` proportional to `strength`. It identifies `sdratio` when
`strength` takes at least two nonzero values and the criterion lies
clearly away from the noise mean (for `"normal"`, a false-alarm rate
clearly away from .5); at a criterion of 0 it does not. The estimate
then rests on the assumption that `d` is zero where `strength` is zero,
not on a second operating point in the data. A numeric predictor under a
log link on `d` barely identifies `sdratio`, and in practice `sdratio`
returns its prior. [`bmm()`](https://venpopov.com/bmm/reference/bmm.md)
warns for every criterion that is the same in every row, so it also
warns for `d ~ 0 + strength`. Giving `d` and `criterion` the same
condition-wise predictors (`d ~ 0 + condition`,
`criterion ~ 0 + condition`) fails too, because each condition then has
two rates for two parameters of its own, and no warning catches it. A
random effect such as `criterion ~ 1 + (1 | id)` moves the boundary
between subjects and does not trigger the warning, but it informs
`sdratio` only weakly, through the assumption that the subject effects
on `d` and on the criterion are independent.

**With `bias`**, any predictor that moves the operating point supplies
what is missing, including a sensitivity manipulation such as
`d ~ 0 + condition` with `bias ~ 1`: holding the bias constant moves the
criterion with the separation, and that movement traces the ROC. The
identification then rests on that assumption about the observer. The
design in which every parameter is intercept-only with no random effects
cannot work: a single (hit, false-alarm) pair is two numbers for three
unknowns, and [`bmm()`](https://venpopov.com/bmm/reference/bmm.md) warns
for it. Neither can `bias` and `d` with the same condition-wise
predictors, for the reason given above, and that design raises no
warning.

Sampling converges and `Rhat` is fine in every unidentified case, so the
warnings are the only signal. They count formula terms, not how far the
terms move the operating point, so a predictor with a small effect
carries little information without a warning. Leaving `sdratio` at its
default is always identified.

## Reading `sdratio` and carrying it to [`dsdt_yn()`](https://venpopov.com/bmm/reference/sdt_yn_dist.md)/[`rsdt_yn()`](https://venpopov.com/bmm/reference/sdt_yn_dist.md)

As in every `bmm` model, the parameters the model *estimates* are on
their link scale, while the distribution functions take their arguments
on the *natural* scale. `sdratio` has a log link, so
[`summary()`](https://rdrr.io/r/base/summary.html) reports \\\log r\\
whereas [`dsdt_yn()`](https://venpopov.com/bmm/reference/sdt_yn_dist.md)
and [`rsdt_yn()`](https://venpopov.com/bmm/reference/sdt_yn_dist.md)
expect the ratio \\r\\ itself (their default is 1, equal variance).
Exponentiate before carrying a posterior value across:

    r <- exp(as_draws_matrix(fit)[, "b_sdratio_Intercept"])

A posterior mean of `sdratio = 0.22` is a ratio of `exp(0.22) = 1.25`.
Passing `0.22` straight to
[`rsdt_yn()`](https://venpopov.com/bmm/reference/sdt_yn_dist.md) instead
asks for a signal distribution 4.5 times *narrower* than the noise — a
legal value that raises no error, and the one mistake worth checking for
in a posterior predictive check written by hand. `d`, `criterion` and
`bias` have identity links, so they carry across unchanged.

The **zROC slope** reported in the recognition-memory literature is the
reciprocal of that ratio, `1 / exp(sdratio)`, so a `sdratio` posterior
mean of 0.375 is a zROC slope of 0.69.

The same log link applies going the other way: a constant you supply
yourself, whether as `bmf(sdratio = )` or through a hand-written
`brms::set_prior(..., dpar = "sdratio")`, is read on it too.
`bmf(sdratio = 1)` does not fix a ratio of 1 — it fixes `exp(1) = 2.72`,
and [`bmm()`](https://venpopov.com/bmm/reference/bmm.md) raises no
warning; a fixed ratio of 1.25 needs `bmf(sdratio = log(1.25))`. Both
`sdratio` defaults
[`default_prior()`](https://venpopov.com/bmm/reference/default_prior.bmmformula.md)
reports — `normal(0, 0.3)` on the intercept and `exponential(2)` on the
random-effect SDs — are on that same scale, unannotated.

## Terms used on this page

- **noise-standardized axis** — latent evidence expressed in units of
  the noise distribution's SD, with its origin at the noise mean.
  `criterion` is a position on this axis; `bias` is a distance on it,
  measured from the midpoint between the means.

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

[`sdt_d()`](https://venpopov.com/bmm/reference/SDTdist.md),
[`sdt_criterion()`](https://venpopov.com/bmm/reference/SDTdist.md) and
[`sdt_bias()`](https://venpopov.com/bmm/reference/SDTdist.md) compute
the `d`, `criterion` and `bias` of this model in closed form from a
single pair of observed hit and false-alarm rates, without fitting: use
them for a quick check of a fitted value, and this model when you need a
hierarchical or condition-wise estimate, or unequal variance.

## Examples

``` r
if (FALSE) { # \dontrun{
dat <- expand.grid(id = 1:20, stimulus = c(0L, 1L))
dat$n_trials <- 100L
dat$n_old <- rsdt_yn(nrow(dat), dat$n_trials, dat$stimulus,
                     d = 1.5, criterion = 0.9)

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

# The same observer, with the boundary measured from the midpoint
fit_bias <- bmm(
  formula = bmf(d ~ 1 + (1 | id), bias ~ 1 + (1 | id)),
  data = dat,
  model = model,
  cores = 4,
  backend = "cmdstanr"
)

# Weak and strong targets against one shared set of lures, one criterion per
# subject: `strong` is 1 for strong targets and 0 for weak targets and lures
# (the lure rows' value does not matter, their false-alarm rate depends on
# the criterion alone)
lure_dat <- expand.grid(id = 1:20, cell = c("lure", "weak", "strong"))
lure_dat$stimulus <- as.integer(lure_dat$cell != "lure")
lure_dat$strong <- as.integer(lure_dat$cell == "strong")
lure_dat$n_trials <- 100L
lure_dat$n_old <- rsdt_yn(nrow(lure_dat), lure_dat$n_trials,
                          lure_dat$stimulus, d = 1 + 0.8 * lure_dat$strong,
                          criterion = 0.7)
fit_strength <- bmm(
  formula = bmf(d ~ 1 + strong + (1 | id), criterion ~ 1 + (1 | id)),
  data = lure_dat,
  model = model,
  cores = 4,
  backend = "cmdstanr"
)

# Simulating from a fitted unequal-variance model: the model reports
# sdratio on its log link, rsdt_yn() takes the ratio itself.
sdratio_posterior <- 0.375
rsdt_yn(2, 100L, c(0L, 1L),
        d = 1.3, criterion = 0.75, sdratio = exp(sdratio_posterior))
} # }
```
