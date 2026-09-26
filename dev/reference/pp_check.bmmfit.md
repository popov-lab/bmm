# Posterior predictive check for bmmfit objects

For models where brms provides `pp_check` support, this method delegates
to
[`brms::pp_check()`](https://mc-stan.org/bayesplot/reference/pp_check.html).
For models with multinomial families (e.g., the m3 model), brms's
`pp_check` is unavailable; this method dispatches to a model-specific
visualisation instead.

## Usage

``` r
# S3 method for class 'bmmfit'
pp_check(
  object,
  type = NULL,
  ndraws = NULL,
  group = NULL,
  resp_var = NULL,
  ...
)
```

## Arguments

- object:

  A `bmmfit` object returned by
  [`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md).

- type:

  Character. Type of pp_check. When `NULL` (default), resolves to
  `"dens_overlay"`, or to the selected observable's default type when
  `resp_var` is specified. When `group` is specified, the grouped
  variant (e.g., `"dens_overlay_grouped"`) is auto-selected if
  available. Multinomial models produce a response proportion profile
  regardless of the value supplied. With `resp_var`,
  `type = "bars_binned"` is also available: it bins a continuous
  statistic like a histogram, with bars for the observed number of
  observations per bin and points with intervals for the predicted
  number. It is the default for the
  [`ezdm()`](https://venpopov.com/bmm/dev/reference/ezdm.md) accuracy
  check.

- ndraws:

  Integer. Number of posterior draws. Defaults to `100` for multinomial
  models and `10` when `resp_var` is specified; otherwise passed to
  [`brms::pp_check()`](https://mc-stan.org/bayesplot/reference/pp_check.html).

- group:

  Character. Optional grouping variable for faceting. For
  non-multinomial models, passed to
  [`brms::pp_check()`](https://mc-stan.org/bayesplot/reference/pp_check.html);
  when specified, the grouped variant of `type` (e.g.,
  `"dens_overlay_grouped"`) is auto-selected if available. For
  multinomial models, facets by the named predictor.

- resp_var:

  Character. For models that declare several observables, the name of
  the observable to check, or `"all"` for a panel of all available
  checks built from one shared simulation. See
  [`pp_check_vars()`](https://venpopov.com/bmm/dev/reference/pp_check_vars.md)
  for the options of a fitted model. The default `NULL` checks the
  primary response via
  [`brms::pp_check()`](https://mc-stan.org/bayesplot/reference/pp_check.html).
  For the RT models, passing `negative_rt = TRUE` (a
  [`brms::posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html)
  argument) is redirected to `resp_var = "signed_rt"`, so that observed
  and predicted response times are both signed by the response.

- ...:

  Additional arguments. Without `resp_var`, forwarded to
  [`brms::pp_check()`](https://mc-stan.org/bayesplot/reference/pp_check.html),
  or for multinomial models to
  [`brms::posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html)
  (`probs`, a numeric vector of length 2 with default `c(0.025, 0.975)`,
  sets the credible interval). With `resp_var`, `draw_ids` and
  `re_formula` go to
  [`brms::prepare_predictions()`](https://paulbuerkner.com/brms/reference/prepare_predictions.html)
  and the rest to the `bayesplot::ppc_*` function.
  `type = "bars_binned"` takes `breaks` (bin edges that cover the
  observed and predicted values), `prob` (interval width, default `0.9`)
  and `freq` (`FALSE` for proportions instead of counts).
  `re_formula = NA` predicts at the population level on every path.

## Value

For multinomial models or when `resp_var` is specified, a `ggplot2`
object (a `bayesplot_grid` for `resp_var = "all"`). For other models,
the result of
[`brms::pp_check()`](https://mc-stan.org/bayesplot/reference/pp_check.html).

## Details

For **multinomial models**, the plot mirrors the bayesplot `ppc_bars`
style: observed proportions are shown as bars and posterior predictive
medians with credible intervals are shown as point-ranges, using the
bayesplot default colour scheme and theme.

Some models describe several observables jointly (e.g.
[`ddm()`](https://venpopov.com/bmm/dev/reference/ddm.md): response times
*and* responses;
[`ezdm()`](https://venpopov.com/bmm/dev/reference/ezdm.md): mean RT, RT
variance *and* accuracy), but brms only ever checks the primary
response. For these models the `resp_var` argument selects which
observable to check;
[`pp_check_vars()`](https://venpopov.com/bmm/dev/reference/pp_check_vars.md)
lists the available checks. All selected observables are drawn from
**one** joint posterior predictive simulation, so `resp_var = "all"`
panels are mutually consistent.

Some observables are undefined for some cells — an
[`ezdm()`](https://venpopov.com/bmm/dev/reference/ezdm.md) boundary has
no mean response time when fewer than two responses reach it.
Observations whose *observed* value is undefined are dropped from the
check; undefined values in the *simulated* replicates are absorbed by
dropping those posterior draws instead, so the number of observations
checked does not depend on `ndraws`. Dropped draws are not missing at
random — they are draws whose parameters made a boundary sparse — so the
retained predictive is mildly conditioned; both reductions are reported
with a warning. With `resp_var = "all"` one reduction is shared by every
panel, so the panels are computed on the same observations and draws.

## See also

[`brms::pp_check()`](https://mc-stan.org/bayesplot/reference/pp_check.html),
[`pp_check_vars()`](https://venpopov.com/bmm/dev/reference/pp_check_vars.md)
