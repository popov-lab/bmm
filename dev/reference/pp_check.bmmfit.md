# Posterior predictive check for bmmfit objects

For models where brms provides `pp_check` support, this method delegates
to
[`brms::pp_check()`](https://paulbuerkner.com/brms/reference/pp_check.brmsfit.html).
For models with multinomial families (e.g., the m3 model), brms's
`pp_check` is unavailable; this method dispatches to a model-specific
visualisation instead.

## Usage

``` r
# S3 method for class 'bmmfit'
pp_check(object, type = "dens_overlay", ndraws = NULL, group = NULL, ...)
```

## Arguments

- object:

  A `bmmfit` object returned by
  [`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md).

- type:

  Character. Type of pp_check (default `"dens_overlay"`). For
  non-multinomial models, passed to
  [`brms::pp_check()`](https://paulbuerkner.com/brms/reference/pp_check.brmsfit.html).
  When `group` is specified, the grouped variant (e.g.,
  `"dens_overlay_grouped"`) is auto-selected if available. Multinomial
  models produce a response proportion profile regardless of the value
  supplied.

- ndraws:

  Integer. Number of posterior draws. Defaults to `100` for multinomial
  models; otherwise passed to
  [`brms::pp_check()`](https://paulbuerkner.com/brms/reference/pp_check.brmsfit.html).

- group:

  Character. Optional grouping variable for faceting. For
  non-multinomial models, passed to
  [`brms::pp_check()`](https://paulbuerkner.com/brms/reference/pp_check.brmsfit.html);
  when specified, the grouped variant of `type` (e.g.,
  `"dens_overlay_grouped"`) is auto-selected if available. For
  multinomial models, facets by the named predictor.

- ...:

  Additional arguments forwarded to
  [`brms::pp_check()`](https://paulbuerkner.com/brms/reference/pp_check.brmsfit.html)
  (non-multinomial) or to
  [`brms::posterior_predict()`](https://paulbuerkner.com/brms/reference/posterior_predict.brmsfit.html)
  (multinomial). For multinomial models, `probs` (numeric vector of
  length 2, default `c(0.025, 0.975)`) controls the credible interval.
  Both model types accept `re_formula` (e.g., `re_formula = NA` to
  predict at the population level, excluding random effects).

## Value

For multinomial models, a `ggplot2` object. For other models, the result
of
[`brms::pp_check()`](https://paulbuerkner.com/brms/reference/pp_check.brmsfit.html).

## Details

For **multinomial models**, the plot mirrors the bayesplot `ppc_bars`
style: observed proportions are shown as bars and posterior predictive
medians with credible intervals are shown as point-ranges, using the
bayesplot default colour scheme and theme.

## See also

[`brms::pp_check()`](https://paulbuerkner.com/brms/reference/pp_check.brmsfit.html)
