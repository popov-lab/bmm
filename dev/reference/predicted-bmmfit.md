# Predictions of a model parameter for easystats and marginaleffects

With these methods, insight and marginaleffects, and through them
`modelbased::estimate_means()`, `modelbased::estimate_contrasts()`,
`modelbased::estimate_relation()` and the see plots of their results,
predict a model parameter on its native scale, as
[`native_parameters()`](https://venpopov.com/bmm/dev/reference/native_parameters.md)
computes it: with the inverse link applied draw by draw and the mixture
weights softmaxed jointly, so that `predict = "thetat"` returns the
probability of a target response.

## Usage

``` r
# S3 method for class 'bmmfit'
get_predicted(
  x,
  data = NULL,
  predict = NULL,
  iterations = NULL,
  ci = NULL,
  ci_method = NULL,
  include_random = "default",
  ...
)

# S3 method for class 'bmmfit'
get_predict(model, newdata, type = "response", ...)
```

## Arguments

- x, model:

  A `bmmfit` object (created by
  [`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md)).

- data, newdata:

  The data to predict for. `get_predicted()` also takes `newdata` in
  `...`, and predicts for the data of the fit if both are `NULL`.

- predict:

  The name of a model parameter, or one of the prediction types of
  `insight::get_predicted()`. A `type` in `...` replaces it.

- iterations:

  The number of posterior draws to predict from, all of them if `NULL`.

- ci:

  The probability mass of the credible intervals.

- ci_method:

  `"quantile"` (the default), `"eti"` or `"hdi"`.

- include_random:

  `FALSE` or `NA` for population-level predictions, a formula for the
  random effects to include, or `"default"` or `TRUE` for all of them if
  the data contain the grouping factors.

- ...:

  For a model parameter: `draw_ids`, `allow_new_levels` and
  `sample_new_levels`, passed to
  [`native_parameters()`](https://venpopov.com/bmm/dev/reference/native_parameters.md);
  for `get_predicted()` also `newdata`, `type` and
  `centrality_function`, and for `get_predict()` `dpar`, `ndraws` and
  `re_formula` (or `re.form`). For any other prediction, passed to the
  methods for brms fits.

- type:

  The prediction type of marginaleffects, used only when `dpar` is not a
  model parameter.

## Value

`get_predicted()`: a `get_predicted` object, as described in
`insight::get_predicted()`. `get_predict()`: a `data.table` of the
estimates with the draws in the `posterior_draws` attribute, as
marginaleffects expects.

## Details

`predict` (or `dpar` for `get_predict()`, which is how modelbased passes
it on) names the model parameter. Without it, `get_predicted()` predicts
the first parameter the model estimates and says which. The prediction
types of insight (`"expectation"`, `"link"`, `"prediction"`, ...) and of
marginaleffects are left to the methods for brms fits; the expected
response is not defined for most bmm models.

Random effects enter the predictions only if every grouping factor is a
column of the data to predict for and holds only levels of the fit, as
insight does for brms fits, or new levels with
`allow_new_levels = TRUE`. Otherwise the predictions are for the
population level, with a warning if you asked for
`include_random = TRUE`. A grouping factor that is `NA`, as in the grids
of modelbased, always means the population level. For a mixed model,
`modelbased::estimate_means()` averages the predictions over the
participants in its grid, so its means differ from those of
`native_parameters(re_formula = NA)`.

The point prediction of `get_predicted()` is the `centrality_function`
of the draws, [`base::mean`](https://rdrr.io/r/base/mean.html) by
default, as for brms fits.

## See also

[`native_parameters()`](https://venpopov.com/bmm/dev/reference/native_parameters.md),
[easystats-bmmfit](https://venpopov.com/bmm/dev/reference/easystats-bmmfit.md)

## Examples

``` r
if (FALSE) { # \dontrun{
insight::get_predicted(fit, predict = "c")
modelbased::estimate_means(fit, by = "set_size", predict = "c")
modelbased::estimate_contrasts(fit, contrast = "set_size", predict = "c")
} # }
```
