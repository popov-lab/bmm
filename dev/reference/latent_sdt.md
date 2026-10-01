# Model-implied latent decision-variable distributions

Computes the noise and signal latent decision-variable densities implied
by a fitted
[`sdt_yn()`](https://popov-lab.github.io/bmm/dev/reference/sdt_yn.md) or
[`sdt_rating()`](https://popov-lab.github.io/bmm/dev/reference/sdt_rating.md)
model, together with the response criterion (binary) or the K-1
confidence thresholds (rating). This is the canonical signal detection
picture: two evidence distributions separated by the sensitivity and cut
by one or more decision boundaries.

## Usage

``` r
latent_sdt(
  fit,
  conditions = NULL,
  n_grid = 200,
  probs = c(0.025, 0.975),
  collapse = NULL,
  show_competitors = FALSE,
  ...
)
```

## Arguments

- fit:

  A `bmmfit` object returned by
  [`bmm()`](https://popov-lab.github.io/bmm/dev/reference/bmm.md) from
  an SDT model.

- conditions:

  Optional data frame of predictor values at which to evaluate the
  model. Column names must match predictor variables used in the
  formula. If `NULL` (default), unique predictor combinations are
  derived from the data.

- n_grid:

  Integer. Number of points on the evidence-axis grid at which each
  density is evaluated (default 200).

- probs:

  Numeric vector of length 2. Lower and upper quantiles for the credible
  band on the criterion/threshold locations (default `c(0.025, 0.975)`).

- collapse:

  Control of which predictor dimensions are collapsed into one panel.
  `NULL` (default) auto-detects: predictors of the criterion/thresholds
  only (and the set size) are collapsed; predictors of `d`/`sdratio` are
  faceted. Pass a character vector of column names to force those
  columns to collapse, or `FALSE` to facet every dimension (one panel
  per combination).

- show_competitors:

  Logical (default `FALSE`). For `sdt_mafc`/`sdt_ranking` only,
  additionally overlay the density of the maximum of the `m - 1`
  distractor samples (one curve per set size) – the "effective
  competitor" the target must beat – which shifts rightward as `m` grows
  and visualises why accuracy falls with set size. The overlay shows the
  set sizes in the fitted data, whatever `conditions` requests. Ignored
  for `sdt_yn`/`sdt_rating`.

- ...:

  Additional arguments passed to
  [`brms::posterior_linpred()`](https://mc-stan.org/rstantools/reference/posterior_linpred.html),
  such as `draw_ids` to use a subset of the posterior draws. `ndraws` is
  refused: each parameter is drawn by its own call, so random subsets
  would not match across parameters.

## Value

A data frame of class `"bmm_sdt_latent"` with columns `x` (the evidence
axis), `density`, `distribution` (`"noise"` or `"signal"`), and any
faceting (density) condition columns. A `lines` attribute holds the
decision-boundary positions (columns `position`, `lower`, `upper`,
`marker`, `level`, plus condition columns), or `NULL` for the
criterion-free `sdt_mafc`/`sdt_ranking` models. When
`show_competitors = TRUE`, a `competitors` attribute holds the
max-of-distractors densities. An `extra` attribute carries the
response-process parameters of the
[`sdt_rating()`](https://popov-lab.github.io/bmm/dev/reference/sdt_rating.md)
versions as a posterior summary (`parameter`, `mean`, `lower`, `upper`,
plus condition columns): for `version = "dpsdt"`, `Ro` and `Rn` on the
probability scale; for `version = "metad"`, the M-ratio (`mratio`) and
meta-d' (`metad`); `NULL` for `version = "standard"`. The object also
carries `probs`, `model_class`, `dist`, `is_rating` and `conditions`
(the conditions of the density panels).

## Details

Densities are evaluated at the posterior-mean parameters of each
condition, on the centred evidence axis the model uses internally: noise
at `-sep/2` with unit SD and signal at `+sep/2` with SD `exp(sdratio)`,
where `sep = d * sqrt((1 + exp(sdratio)^2) / 2)` is the separation in
noise-SD units. `d` is \\d_a\\, measured in root-mean-square SD units,
so `sep` equals `d` whenever `sdratio` is fixed at 0. The decision
boundaries additionally carry a posterior credible band on their
location. The area beyond a boundary equals the corresponding
hit/false-alarm rate for every noise distribution: the likelihood
evaluates the survivor function of the evidence distribution at the
boundary, which is exactly that area.

Available for all four SDT models.
[`sdt_yn()`](https://popov-lab.github.io/bmm/dev/reference/sdt_yn.md)
draws the criterion and
[`sdt_rating()`](https://popov-lab.github.io/bmm/dev/reference/sdt_rating.md)
the K-1 confidence thresholds as boundary lines.
[`sdt_mafc()`](https://popov-lab.github.io/bmm/dev/reference/sdt_mafc.md)
and
[`sdt_ranking()`](https://popov-lab.github.io/bmm/dev/reference/sdt_ranking.md)
have no response criterion – the decision is a max/rank rule over the
`m` alternatives – so only the densities are drawn (no boundary lines);
there the noise density represents each of the `m - 1` distractor
alternatives and the signal density the target.

Because the densities depend only on `d`/`sdratio`, predictors that vary
the criterion/thresholds *only* (e.g. a base-rate manipulation, as in
[broeder_schuetz_2009_e3](https://popov-lab.github.io/bmm/dev/reference/broeder_schuetz_2009_e3.md))
leave the densities unchanged: by default they are collapsed into a
single panel with their several boundaries overlaid and colour-coded,
rather than shown as repeated identical panels. Predictors that vary
`d`/`sdratio` produce distinct density panels (faceted). For
`sdt_mafc`/`sdt_ranking` the set size is likewise collapsed unless it
predicts `d`. Use `collapse` to override the auto-detection.

## See also

[`roc_sdt()`](https://popov-lab.github.io/bmm/dev/reference/roc_sdt.md),
[`plot.bmm_sdt_latent()`](https://popov-lab.github.io/bmm/dev/reference/plot.bmm_sdt_latent.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Three base-rate conditions shift the criterion, which identifies sdratio
dat <- expand.grid(id = 1:20, stimulus = c(0L, 1L),
                   condition = c("liberal", "neutral", "strict"))
dat$n_trials <- 100L
criteria <- c(liberal = -0.5, neutral = 0, strict = 0.5)
dat$n_old <- rsdt_yn(nrow(dat), dat$n_trials, dat$stimulus, d = 1.5,
                     criterion = criteria[as.character(dat$condition)],
                     sdratio = 1.3)

fit <- bmm(
  formula = bmf(d ~ 1, criterion ~ 0 + condition, sdratio ~ 1),
  data = dat,
  model = sdt_yn(response = "n_old", stimulus = "stimulus",
                 n_trials = "n_trials"),
  cores = 4,
  backend = "cmdstanr"
)

latent <- latent_sdt(fit)
latent
plot(latent)
} # }
```
