# Area under the ROC curve from a fitted SDT model

Computes the posterior area under the ROC curve (AUC). For Gaussian and
Gumbel (min or max) equal-variance binary SDT the AUC is available in
closed form from the `d` draws; otherwise it is obtained by trapezoidal
integration of the model-implied curve (for rating fits, the curve swept
from the posterior of `d` and `sdratio`, lifted by `Ro`/`Rn` for
`version = "dpsdt"`, so the AUC is the area under the predicted mixture
ROC and not \\\Phi(d/\sqrt{2})\\; for `version = "metad"` it is the
type-1 curve). The returned AUC is always the area under the full curve,
not the trapezoid of the discrete operating points or the K-1 rating
thresholds; for binary multi-criteria fits it is one value per curve.

## Usage

``` r
auc_sdt(
  fit,
  conditions = NULL,
  probs = c(0.025, 0.975),
  criterion_points = NULL,
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

- probs:

  Numeric vector of length 2. Quantiles for the credible interval
  (default `c(0.025, 0.975)`).

- criterion_points:

  Optional. Control of the binary multi-criteria behaviour. `NULL`
  (default) auto-detects predictors that vary the criterion only (not
  `d`/`sdratio`) and treats their levels as operating points on one
  curve. Pass a character vector of column names to force that
  classification, or `FALSE` to disable it (one separate curve per
  predictor combination). Ignored for rating models.

- ...:

  Additional arguments passed to
  [`brms::posterior_linpred()`](https://mc-stan.org/rstantools/reference/posterior_linpred.html),
  such as `draw_ids` to use a subset of the posterior draws. `ndraws` is
  refused: each parameter is drawn by its own call, so random subsets
  would not match across parameters.

## Value

A data frame of class `"bmm_sdt_auc"` with columns `AUC`, `.draw`, and
any condition columns, plus a `summary` attribute (`AUC_mean`,
`AUC_lower`, `AUC_upper`), `model_class`, `dist` and `conditions`
attributes.

## Details

The closed form is used only when `sdratio` is fixed at 0, where `d`
(which is \\d_a\\) equals \\d'\\; every fit with `sdratio` estimated or
fixed away from 0 takes the numerical route, so the AUC is invariant to
the sensitivity parameterization.

Analytical formulas (equal variance, where \\d_a = d'\\): normal EV-SDT
\\AUC = \Phi(d'/\sqrt{2})\\; Gumbel-min and Gumbel-max EV-SDT \\AUC =
\mathrm{logistic}(g')\\.

## See also

[`roc_sdt()`](https://popov-lab.github.io/bmm/dev/reference/roc_sdt.md),
[`plot.bmm_sdt_auc()`](https://popov-lab.github.io/bmm/dev/reference/plot.bmm_sdt_auc.md)

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

auc <- auc_sdt(fit)
auc
plot(auc)
} # }
```
