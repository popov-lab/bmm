# Confidence Rating Signal Detection Theory Model

Confidence Rating Signal Detection Theory Model

## Usage

``` r
sdt_rating(
  response,
  stimulus,
  dist = c("normal", "gumbel_min", "gumbel_max", "logistic"),
  threshold_type = c("parsimonious", "equidistant", "log_distance", "log_ratio",
    "softmax"),
  version = c("standard", "dpsdt", "metad"),
  links = NULL,
  ...
)
```

## Arguments

- response:

  A character vector of K column names containing response counts per
  rating category, ordered from "definitely noise" to "definitely
  signal".

- stimulus:

  The name of the variable coding the stimulus type. Must be coded as 0
  (noise/new) and 1 (signal/old).

- dist:

  The distribution assumed for the latent evidence, given here by its
  cumulative distribution function. One of:

  - "normal" (default): Gaussian SDT, \\\Phi(x)\\

  - "gumbel_min": smallest-extreme-value SDT, \\1 - \exp(-\exp(x))\\
    (complementary log-log)

  - "gumbel_max": largest-extreme-value SDT, \\\exp(-\exp(-x))\\
    (log-log, as in
    [`evd::pgumbel`](https://rdrr.io/pkg/evd/man/gumbel.html))

  - "logistic": logistic SDT, \\1 / (1 + \exp(-x))\\

- threshold_type:

  Character. Threshold parameterization:

  - "parsimonious" (default): 2 parameters (criterion + spacing).
    Thresholds follow logit-spaced canonical positions (Selker et al.,
    2019).

  - "equidistant": 2 parameters (criterion + spacing). Thresholds are
    equally spaced, exp(spacing) apart.

  - "log_distance": K-2 parameters. `delta<i>` is the log width of the
    interval between thresholds i and i + 1, so the intervals are free
    and the ordering is guaranteed (Paulewicz & Blaut, 2022).

  - "log_ratio": K-2 parameters (Paulewicz & Blaut, 2022; the odd-K form
    is bmm's). One interval is the spread, `exp(delta)`: for even K the
    interval just above `criterion`, for odd K the middle category. The
    first interval on the other side (even K) or on each side (odd K) is
    a ratio times the spread, and every further interval a ratio times
    the first interval on its side. The deltas are therefore not
    exchangeable; `model$parameters` names the role of each.

  - "softmax": K-2 parameters. A shared spacing parameter sets the mean
    interval width, exp(spacing), while K-3 delta parameters share the
    total width out over the intervals through a softmax (each delta is
    the log ratio of its interval to the last one).

- version:

  Character. The latent-process version of the rating model. One of
  `"standard"` (default, a single familiarity SDT process), `"dpsdt"`
  (dual-process SDT with recollection parameters `Ro`/`Rn`), or
  `"metad"` (meta-d' with a metacognitive efficiency parameter
  `logmratio`, the log M-ratio). All three use the same
  confidence-rating response interface. See Details.

- links:

  A named list of link functions for the parameters, e.g.
  `links = list(d = "log")`. Only `d` and `criterion` can be set, to
  `"identity"`, `"log"`, `"softplus"`, `"logit"` or `"probit"`.
  `sdratio` and the threshold parameters keep their identity links,
  because the model reads each of them through
  [`exp()`](https://rdrr.io/r/base/Log.html) and fixes `sdratio` at 0
  for equal variance. The same holds for `Ro` and `Rn`, read through
  `inv_logit()` and fixed off by default, and for `logmratio`, read
  through [`exp()`](https://rdrr.io/r/base/Log.html).

- ...:

  used internally for testing, ignore it

## Value

An object of class `bmmodel`

## Details

Three versions share the confidence-rating response interface, selected
with `version`:

### Version: `standard` (default)

- **Domain:** Perception & Recognition Memory

- **Task:** Signal/Noise or Old/New Recognition

- **Name:** Signal Detection Theory (Confidence Rating)

- **Citation:**

  - Green, D. M., & Swets, J. A. (1966). Signal detection theory and
    psychophysics. Wiley.

- **Version:** standard

- **Requirements:**

  Provide pre-aggregated data with the following columns:

&nbsp;

- Response counts: one column per rating category (K columns)

- Stimulus type (stimulus): 0 = noise, 1 = signal Categories should be
  ordered: 1 = 'definitely noise' to K = 'definitely signal'

&nbsp;

- **Parameters:**

  - `d`: Sensitivity: d' under equal variance (the default). When
    sdratio is estimated, d is d_a, the distance between the signal and
    noise distributions in units of their root-mean-square SD

  - `criterion`: Response bias, on the noise-standardized axis: the
    middle threshold (the old/new boundary) for an even number of
    categories, the centre of the middle category for an odd number

  - `spacing`: Threshold spacing: controls distance between adjacent
    thresholds (exp(spacing) ensures positive spacing)

  - `sdratio`: Log SD ratio: the log of the signal-to-noise standard
    deviation ratio, so exp(sdratio) is the ratio itself and 0 means
    equal variance

- **Fixed parameters:**

  - `sdratio` = 0

- **Default parameter links:**

  - d = identity; criterion = identity; spacing = identity; sdratio =
    identity

- **Default priors:**

  - `d`:

    - `main`: normal(1, 1)

    - `effects`: normal(0, 0.5)

    - `sd`: exponential(1)

  - `criterion`:

    - `main`: normal(0, 1.5)

    - `effects`: normal(0, 0.5)

    - `sd`: exponential(2)

  - `spacing`:

    - `main`: normal(0, 0.5)

    - `effects`: normal(0, 0.3)

    - `sd`: exponential(2)

  - `sdratio`:

    - `main`: normal(0, 0.3)

    - `effects`: normal(0, 0.3)

    - `sd`: exponential(2)

### Version: `dpsdt`

Dual-process SDT (Yonelinas, 1994): a familiarity SDT process plus an
all-or-none recollection process that loads the most-confident category.
`Ro` is recollection of old items (loads the most-confident "signal"
category), the only recollection term of the classic model; `Rn` is
recall-to-reject of new items (loads the most-confident "noise"
category), the two-sided extension that Yonelinas (2024) applies to
visual working memory. `inv_logit(Ro)`/`inv_logit(Rn)` are the
recollection probabilities. Both are fixed off by default (recovering
`standard`); add `Ro ~ 1` for the one-sided model and `Ro ~ 1, Rn ~ 1`
for the two-sided model. `d` describes the familiarity distributions
only – the observed ROC is a mixture of familiarity and recollection, so
it is not the discriminability of that mixture.
[`summary()`](https://rdrr.io/r/base/summary.html) reports `Ro` and `Rn`
on the logit scale;
[`latent_sdt()`](https://popov-lab.github.io/bmm/dev/reference/latent_sdt.md)
returns them as probabilities (attribute `extra`). Two cautions. Freeing
`sdratio` alongside `Ro` is weakly identified from a single ROC: in
simulation the two are correlated at about -0.8 in the posterior and
`Ro` is pulled down while `exp(sdratio)` is pulled above 1; keep the
familiarity process equal-variance unless the design separates them. And
a recollection probability near zero is reported as the tail of its
prior: on the probability scale the interval cannot include 0, so the
test of "no recollection" is a comparison with the fit that leaves the
parameter fixed off, not the interval.

- **Domain:** Perception & Recognition Memory

- **Task:** Signal/Noise or Old/New Recognition

- **Name:** Signal Detection Theory (Confidence Rating)

- **Citation:**

  - Yonelinas, A. P. (1994). Receiver-operating characteristics in
    recognition memory: Evidence for a dual-process model. Journal of
    Experimental Psychology: Learning, Memory, and Cognition, 20(6),
    1341-1354. https://doi.org/10.1037/0278-7393.20.6.1341

  - Yonelinas, A. P. (2024). The role of recollection and familiarity in
    visual working memory: A mixture of threshold and signal detection
    processes. Psychological Review, 131(2), 321-348.
    https://doi.org/10.1037/rev0000432

- **Version:** dpsdt

- **Requirements:**

  Provide pre-aggregated data with the following columns:

&nbsp;

- Response counts: one column per rating category (K columns)

- Stimulus type (stimulus): 0 = noise, 1 = signal Categories should be
  ordered: 1 = 'definitely noise' to K = 'definitely signal'

&nbsp;

- **Parameters:**

  - `d`: Sensitivity: d' under equal variance (the default). When
    sdratio is estimated, d is d_a, the distance between the signal and
    noise distributions in units of their root-mean-square SD

  - `criterion`: Response bias, on the noise-standardized axis: the
    middle threshold (the old/new boundary) for an even number of
    categories, the centre of the middle category for an odd number

  - `spacing`: Threshold spacing: controls distance between adjacent
    thresholds (exp(spacing) ensures positive spacing)

  - `sdratio`: Log SD ratio: the log of the signal-to-noise standard
    deviation ratio, so exp(sdratio) is the ratio itself and 0 means
    equal variance

  - `Ro`: Recollection of old items: probability inv_logit(Ro) that an
    old item is recollected as old, loading the most-confident 'signal'
    category

  - `Rn`: Recollection of new items: probability inv_logit(Rn) that a
    new item is recall-rejected, loading the most-confident 'noise'
    category

- **Fixed parameters:**

  - `sdratio` = 0

  - `Ro` = -100

  - `Rn` = -100

- **Default parameter links:**

  - d = identity; criterion = identity; spacing = identity; sdratio =
    identity; Ro = identity; Rn = identity

- **Default priors:**

  - `d`:

    - `main`: normal(1, 1)

    - `effects`: normal(0, 0.5)

    - `sd`: exponential(1)

  - `criterion`:

    - `main`: normal(0, 1.5)

    - `effects`: normal(0, 0.5)

    - `sd`: exponential(2)

  - `spacing`:

    - `main`: normal(0, 0.5)

    - `effects`: normal(0, 0.3)

    - `sd`: exponential(2)

  - `sdratio`:

    - `main`: normal(0, 0.3)

    - `effects`: normal(0, 0.3)

    - `sd`: exponential(2)

  - `Ro`:

    - `main`: normal(0, 1)

    - `effects`: normal(0, 0.5)

    - `sd`: exponential(1)

  - `Rn`:

    - `main`: normal(0, 1)

    - `effects`: normal(0, 0.5)

    - `sd`: exponential(1)

### Version: `metad`

Meta-d' (Maniscalco & Lau, 2012): a type-2 metacognitive sensitivity
governs confidence-threshold placement, with the total "old"/"new"
response rates held to what type-1 `d` implies. Rather than estimating
meta-d' directly, the model estimates `logmratio`, the log M-ratio
\\\log(\mathrm{meta\text{-}d'}/d')\\, and recovers meta-d' as
`exp(logmratio) * d`. The M-ratio is the field-standard measure of
metacognitive efficiency (Maniscalco & Lau, 2012; Fleming, 2017):
estimating it on the log scale keeps meta-d' positive, regularizes it
toward `d`, and anchors the ideal point (meta-d' = `d`, perfect
metacognition) at `logmratio = 0`, which recovers `standard`. Type-1 and
type-2 sensitivity share one scale (\\d'\\, or \\d_a\\ when `sdratio` is
estimated), so the M-ratio is unaffected by `sdratio`. The type-1
boundary is `criterion`, the middle threshold, so this version needs an
even number of rating categories: with an odd number the middle category
straddles the boundary (see "Where `criterion` sits"). The type-1
criterion is held at the same location in the meta-d' space as in the
type-1 model, as the HMeta-d toolbox's model equations do (Fleming,
2017, Appendix). Maximum-likelihood meta-d' instead constrains meta-c' =
c' (Maniscalco & Lau, 2014), so its estimates differ from bmm's when the
criterion is far from the midpoint. Extract the M-ratio posterior with
[`mratio()`](https://popov-lab.github.io/bmm/dev/reference/mratio.md).

- **Domain:** Perception & Recognition Memory

- **Task:** Signal/Noise or Old/New Recognition

- **Name:** Signal Detection Theory (Confidence Rating)

- **Citation:**

  - Maniscalco, B., & Lau, H. (2012). A signal detection theoretic
    approach for estimating metacognitive sensitivity from confidence
    ratings. Consciousness and Cognition, 21(1), 422-430.
    https://doi.org/10.1016/j.concog.2011.09.021

- **Version:** metad

- **Requirements:**

  Provide pre-aggregated data with the following columns:

&nbsp;

- Response counts: one column per rating category (K columns)

- Stimulus type (stimulus): 0 = noise, 1 = signal Categories should be
  ordered: 1 = 'definitely noise' to K = 'definitely signal'

&nbsp;

- **Parameters:**

  - `d`: Sensitivity: d' under equal variance (the default). When
    sdratio is estimated, d is d_a, the distance between the signal and
    noise distributions in units of their root-mean-square SD

  - `criterion`: Response bias, on the noise-standardized axis: the
    middle threshold (the old/new boundary) for an even number of
    categories, the centre of the middle category for an odd number

  - `spacing`: Threshold spacing: controls distance between adjacent
    thresholds (exp(spacing) ensures positive spacing)

  - `sdratio`: Log SD ratio: the log of the signal-to-noise standard
    deviation ratio, so exp(sdratio) is the ratio itself and 0 means
    equal variance

  - `logmratio`: Log M-ratio, log(meta-d/d): metacognitive efficiency. 0
    is ideal metacognition (meta-d = d), negative is inefficiency,
    positive is hyper-efficiency. meta-d is recovered as exp(logmratio)
    \* d, on the same scale as d

- **Fixed parameters:**

  - `sdratio` = 0

- **Default parameter links:**

  - d = identity; criterion = identity; spacing = identity; sdratio =
    identity; logmratio = identity

- **Default priors:**

  - `d`:

    - `main`: normal(1, 1)

    - `effects`: normal(0, 0.5)

    - `sd`: exponential(1)

  - `criterion`:

    - `main`: normal(0, 1.5)

    - `effects`: normal(0, 0.5)

    - `sd`: exponential(2)

  - `spacing`:

    - `main`: normal(0, 0.5)

    - `effects`: normal(0, 0.3)

    - `sd`: exponential(2)

  - `sdratio`:

    - `main`: normal(0, 0.3)

    - `effects`: normal(0, 0.3)

    - `sd`: exponential(2)

  - `logmratio`:

    - `main`: normal(0, 0.5)

    - `effects`: normal(0, 0.3)

    - `sd`: exponential(2)

By default, the model assumes equal variance (sdratio fixed to 0). To
estimate unequal variance, add `sdratio ~ 1` (or `sdratio ~ predictors`)
to the formula.

## Sensitivity is on the same scale as [`sdt_yn()`](https://popov-lab.github.io/bmm/dev/reference/sdt_yn.md)

`d` is \\d'\\ whenever `sdratio` stays fixed at 0, which is every fit
that does not give `sdratio` a formula. With `sdratio` estimated, `d` is
the balanced index \\d_a\\ that
[`sdt_yn()`](https://popov-lab.github.io/bmm/dev/reference/sdt_yn.md)
reports: the separation between the signal and noise distributions
divided by the root-mean-square of their SDs. Unlike the
noise-standardized \\d'\\, it remains comparable across conditions that
differ in `sdratio`; see the sensitivity section of
[`sdt_yn()`](https://popov-lab.github.io/bmm/dev/reference/sdt_yn.md)
for the reasoning, for how far the two indices lie apart, and for the
caveat that under the Gumbel distributions \\d_a\\ is not the
AUC-equivalent index once `sdratio` is estimated.

The `criterion` and the confidence thresholds are **not** rescaled. They
stay on the noise-standardized axis, so under unequal variance they and
`d` are in different units.

The `dpsdt` and `metad` versions inherit the same convention: `d` there
describes the familiarity (type-1) distributions, and meta-d' is scaled
the same way, which leaves the M-ratio invariant to `sdratio`.

## Where `criterion` sits

Every `threshold_type` places `criterion` the same way. With an even
number of categories K it is the middle threshold, the boundary between
the K/2 "noise" categories and the K/2 "signal" categories. With an odd
number there is no such boundary – the middle category straddles it – so
`criterion` is the centre of that category, and the two thresholds
around it sit half an interval below and above. A shift in `criterion`
therefore moves the whole threshold set, and the threshold parameters
(`spacing`, `delta`) describe its shape around that point.

## Identifying `sdratio`

Unlike
[`sdt_yn()`](https://popov-lab.github.io/bmm/dev/reference/sdt_yn.md), a
rating design identifies `sdratio` from a single condition: the K - 1
thresholds give K - 1 operating points per condition, which trace the
ROC; under `dist = "normal"` its z-transform is a line with slope
`1 / exp(sdratio)`. `sdratio ~ 1` on a one-condition dataset is
identified from K = 3 categories, where it uses every degree of freedom
(the fit is saturated), and leaves the zROC testable from K = 4; the
identification caveats on the
[`sdt_yn()`](https://popov-lab.github.io/bmm/dev/reference/sdt_yn.md)
page do not carry over.

## References

Yonelinas, A. P. (1994). Receiver-operating characteristics in
recognition memory: Evidence for a dual-process model. *Journal of
Experimental Psychology: Learning, Memory, and Cognition*, *20*(6),
1341–1354.

Yonelinas, A. P. (2024). The role of recollection and familiarity in
visual working memory: A mixture of threshold and signal detection
processes. *Psychological Review*, *131*(2), 321–348.
[doi:10.1037/rev0000432](https://doi.org/10.1037/rev0000432)

Maniscalco, B., & Lau, H. (2012). A signal detection theoretic approach
for estimating metacognitive sensitivity from confidence ratings.
*Consciousness and Cognition*, *21*(1), 422–430.

Maniscalco, B., & Lau, H. (2014). Signal detection theory analysis of
type 1 and type 2 data: meta-d', response-specific meta-d', and the
unequal variance SDT model. In S. M. Fleming & C. D. Frith (Eds.), *The
cognitive neuroscience of metacognition* (pp. 25–66). Springer.

Fleming, S. M. (2017). HMeta-d: hierarchical Bayesian estimation of
metacognitive efficiency from confidence ratings. *Neuroscience of
Consciousness*, *2017*(1), nix007.
[doi:10.1093/nc/nix007](https://doi.org/10.1093/nc/nix007)

Green, D. M., & Swets, J. A. (1966). *Signal detection theory and
psychophysics*. Wiley.

Selker, R., van den Bergh, D., Criss, A. H., & Wagenmakers, E.-J.
(2019). Parsimonious estimation of signal detection models from
confidence ratings. *Behavior Research Methods*, *51*(5), 1953–1967.
[doi:10.3758/s13428-019-01231-3](https://doi.org/10.3758/s13428-019-01231-3)

Paulewicz, B., & Blaut, A. (2022). The general causal cumulative model
of ordinal response. *PsyArXiv*.
[doi:10.31234/osf.io/e7a3x](https://doi.org/10.31234/osf.io/e7a3x)

## Examples

``` r
if (FALSE) { # \dontrun{
# EV-SDT rating model
dat <- expand.grid(id = 1:20, stimulus = c(0L, 1L))
dat <- cbind(dat, rsdt_rating(nrow(dat), 200, dat$stimulus,
                              d = 1.5, thresholds = c(-0.5, 0, 0.5)))

model <- sdt_rating(
  response = c("r1", "r2", "r3", "r4"),
  stimulus = "stimulus"
)

fit <- bmm(
  formula = bmf(d ~ 1, criterion ~ 1, spacing ~ 1),
  data = dat,
  model = model,
  cores = 4,
  backend = "cmdstanr"
)

# UV-SDT: add sdratio to the formula
fit_uv <- bmm(
  formula = bmf(d ~ 1, criterion ~ 1, spacing ~ 1, sdratio ~ 1),
  data = dat,
  model = model,
  cores = 4,
  backend = "cmdstanr"
)

# Dual-process SDT: free recollection of old items (one-sided) or both
# old and new items (two-sided) via the formula
model_dp <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus", version = "dpsdt")
fit_dp <- bmm(bmf(d ~ 1, criterion ~ 1, spacing ~ 1, Ro ~ 1, Rn ~ 1),
              data = dat, model = model_dp, backend = "cmdstanr")

# Meta-d': estimate metacognitive efficiency (log M-ratio)
model_md <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus", version = "metad")
fit_md <- bmm(bmf(d ~ 1, criterion ~ 1, spacing ~ 1, logmratio ~ 1),
              data = dat, model = model_md, backend = "cmdstanr")
mratio(fit_md) # posterior M-ratio (meta-d'/d')
} # }
```
