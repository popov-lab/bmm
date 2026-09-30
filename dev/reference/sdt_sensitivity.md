# Sensitivity on the noise, signal, or root-mean-square scale

Re-expresses the posterior sensitivity of a fitted SDT model against a
different reference standard deviation. As long as `sdratio` is fixed at
its default, `d` is the familiar \\d'\\ and all three scales below are
the same number. Once `sdratio` is estimated, `d` is \\d_a\\, which
measures the separation of the two evidence distributions in units of
their root-mean-square SD, and the same separation can also be read
against the noise SD (\\d_N\\, the classical \\d'\\) or against the
signal SD (\\d_S\\). This function returns any of the three as posterior
draws, so contrasts and intervals can be computed on whichever scale a
literature reports.

## Usage

``` r
sdt_sensitivity(
  fit,
  measure = c("da", "dn", "ds"),
  conditions = NULL,
  probs = c(0.025, 0.975),
  ...
)
```

## Arguments

- fit:

  A `bmmfit` object returned by
  [`bmm()`](https://popov-lab.github.io/bmm/dev/reference/bmm.md) from
  an SDT model.

- measure:

  Character vector naming the scales to return: `"da"` (root-mean-square
  SD, the estimated parameter), `"dn"` (noise SD), and/or `"ds"` (signal
  SD). Defaults to all three.

- conditions:

  Optional data frame of predictor values at which to evaluate the
  model. Column names must match predictor variables used in the
  formula. If `NULL` (default), unique predictor combinations are
  derived from the data.

- probs:

  Numeric vector of length 2. Lower and upper quantiles for the credible
  interval (default `c(0.025, 0.975)`).

- ...:

  Additional arguments passed to
  [`brms::posterior_linpred()`](https://mc-stan.org/rstantools/reference/posterior_linpred.html),
  such as `draw_ids` to use a subset of the posterior draws. `ndraws` is
  refused: each parameter is drawn by its own call, so random subsets
  would not match across parameters.

## Value

A data frame of class `"bmm_sdt_sensitivity"` with columns `measure`,
`value`, `.draw`, and any condition columns. The object carries a
`summary` attribute (`measure`, `mean`, `lower`, `upper`, plus condition
columns), and `probs`, `model_class`, `dist` and `conditions`
attributes.

## Details

The model places the noise distribution at \\-\delta/2\\ with SD
\\\sigma_N = 1\\ and the signal distribution at \\+\delta/2\\ with SD
\\\sigma_S = \exp(\mathrm{sdratio})\\, so \\\delta\\ is the separation
in noise-SD units. Dividing that separation by each reference SD gives

\$\$d_N = \delta / \sigma_N, \quad d_S = \delta / \sigma_S, \quad d_a =
\delta / \sqrt{(\sigma_N^2 + \sigma_S^2)/2}.\$\$

The estimated parameter is \\d_a\\, hence \\\delta = d_a \sqrt{(1 +
\sigma_S^2)/2}\\ and

\$\$d_N = d_a \sqrt{(1 + \sigma_S^2)/2}, \qquad d_S = d_a \sqrt{(1 +
\sigma_S^2)/2} \\/\\ \sigma_S.\$\$

All three coincide when `sdratio` is 0 (equal variance), which is the
case for
[`sdt_mafc()`](https://popov-lab.github.io/bmm/dev/reference/sdt_mafc.md)
and for any fit that leaves `sdratio` at its default of 0. The
conversion is applied draw by draw, so the returned intervals propagate
the joint posterior uncertainty in `d` and `sdratio` rather than
combining point estimates.

Only \\d_a\\ is invariant to which distribution is treated as the
reference. \\d_N\\ and \\d_S\\ are not comparable across conditions that
differ in `sdratio`: two conditions with identical discriminability can
show a large, confidently estimated \\d_N\\ difference. Prefer \\d_a\\
for contrasts, and use \\d_N\\/\\d_S\\ for comparison with published
values. When `sdratio` is estimated but constant across the conditions
being compared, all three differ by one common factor per draw and give
the same contrasts up to scale.

The `criterion` and rating thresholds are not converted: they stay on
the noise-SD axis whichever sensitivity scale you report. Under the
Gumbel distributions with `sdratio` estimated, none of the three is the
AUC-equivalent index; use
[`auc_sdt()`](https://popov-lab.github.io/bmm/dev/reference/auc_sdt.md)
to compare discriminability on the probability scale.

## See also

[`auc_sdt()`](https://popov-lab.github.io/bmm/dev/reference/auc_sdt.md),
[`roc_sdt()`](https://popov-lab.github.io/bmm/dev/reference/roc_sdt.md),
[`latent_sdt()`](https://popov-lab.github.io/bmm/dev/reference/latent_sdt.md)

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

# d_a (the fitted d), d_N and d_S, converted draw by draw
sdt_sensitivity(fit)
sdt_sensitivity(fit, measure = "dn", draw_ids = 1:500)
} # }
```
