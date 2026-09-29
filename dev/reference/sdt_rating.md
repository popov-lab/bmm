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

- links:

  A named list of link functions for the parameters, e.g.
  `links = list(d = "log")`. Only `d` and `criterion` can be set, to
  `"identity"`, `"log"`, `"softplus"`, `"logit"` or `"probit"`.
  `sdratio` and the threshold parameters keep their identity links,
  because the model reads each of them through
  [`exp()`](https://rdrr.io/r/base/Log.html) and fixes `sdratio` at 0
  for equal variance.

- ...:

  used internally for testing, ignore it

## Value

An object of class `bmmodel`

## Details

- **Domain:** Perception & Recognition Memory

- **Task:** Signal/Noise or Old/New Recognition

- **Name:** Signal Detection Theory (Confidence Rating)

- **Citation:**

  - Green, D. M., & Swets, J. A. (1966). Signal detection theory and
    psychophysics. Wiley.

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
} # }
```
