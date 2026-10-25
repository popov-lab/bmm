# easystats support for bmmfit objects

These S3 methods let the insight package, and through it the easystats
packages parameters, bayestestR and modelbased, read a bmm fit as the
model you specified with
[`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md): its
parameters, their formulas, the response variables and the data, rather
than the brms model that bmm builds from them. With them, for example,
`parameters::model_parameters()` and `bayestestR::describe_posterior()`
summarise each model parameter in its own component, and for models such
as [`sdm()`](https://venpopov.com/bmm/dev/reference/sdm.md),
`modelbased::estimate_means(fit, by = "set_size", predict = "c")`
returns the estimates of a model parameter on its native scale.

## Usage

``` r
# S3 method for class 'bmmfit'
find_formula(x, verbose = TRUE, ...)

# S3 method for class 'bmmfit'
find_predictors(
  x,
  effects = "fixed",
  component = "all",
  flatten = FALSE,
  verbose = TRUE,
  ...
)

# S3 method for class 'bmmfit'
find_response(x, combine = TRUE, ...)

# S3 method for class 'bmmfit'
get_data(
  x,
  effects = "all",
  component = "all",
  source = "frame",
  verbose = FALSE,
  ...
)

# S3 method for class 'bmmfit'
find_parameters(
  x,
  effects = "all",
  component = "all",
  flatten = FALSE,
  parameters = NULL,
  ...
)

# S3 method for class 'bmmfit'
get_priors(x, verbose = TRUE, ...)

# S3 method for class 'bmmfit'
model_info(x, response = NULL, ...)

# S3 method for class 'bmmfit'
link_function(x, ...)

# S3 method for class 'bmmfit'
link_inverse(x, ...)
```

## Arguments

- x:

  A `bmmfit` object (created by
  [`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md)).

- ...:

  Ignored, for consistency with the generics.

- effects:

  `"fixed"` for the population-level part of the model, `"random"` for
  the group-level part or `"all"` for both: coefficients for
  `find_parameters()`, predictors for `find_predictors()`, columns for
  `get_data()`. `find_parameters()` also takes `"full"`,
  `"random_variance"` and `"grouplevel"`, as the methods of insight for
  brms fits do.

- component:

  `"all"`, `"conditional"`, `"location"`, `"distributional"` or
  `"auxiliary"`, which all select every model parameter, or for
  `find_parameters()` the name of a model parameter; ignored by
  `get_data()`.

- flatten:

  Should a character vector be returned instead of a list?

- combine, source, verbose:

  Ignored, for consistency with the generics.

- parameters:

  A regular expression the names of the returned parameters must match.

- response:

  Ignored, for consistency with the generic.

## Value

See the insight generics, e.g. `insight::find_formula()` and
`insight::find_parameters()`.

## Details

Each method reshapes what the fit already stores:

- `find_formula()` returns one formula per model parameter, and per
  helper parameter of a non-linear formula, with the population-level
  and group-level parts as separate elements (`<parameter>` and
  `<parameter>_random`), plus a `conditional` formula of the response
  variables on all population-level predictors.

- `find_predictors()` returns the predictors of these formulas: the
  population-level predictors for `effects = "fixed"`, the grouping
  factors for `"random"` and both for `"all"`. Every model parameter
  counts as both a location and a distributional parameter, so each
  `component` returns the same variables.

- `find_response()` returns the names of the response columns, as
  `get_data()` returns them.

- `get_data()` returns the data the model was fitted to as you passed
  them: without the columns bmm adds for its own use or brms computes
  from the formula (such as `ID:session` or `log(x)`), and with the
  response categories of a model fitted to counts (e.g.
  [`m3()`](https://venpopov.com/bmm/dev/reference/m3.md)) as columns of
  their own. The model frame stored in the fit is always the source;
  `source` is ignored.

- `find_parameters()` returns the names of the posterior draws in one
  element per model parameter (`<parameter>` for the population-level
  coefficients, `<parameter>_random` for the group-level ones), without
  the parameters that the brms family needs but the model does not have,
  such as the fixed `mu` of the response-time models. `component`
  selects a model parameter by name; `"all"`, `"conditional"`,
  `"location"`, `"distributional"` and `"auxiliary"` select all of them.
  `insight::get_parameters()` and `insight::clean_parameters()` build on
  it.

- `get_priors()` returns one row per population-level coefficient,
  group-level standard deviation and correlation, with the prior that
  applies to it. A parameter fixed to a constant has the distribution
  `"constant"` and the constant as its location.

- `model_info()` describes the model as not linear and as having no
  dispersion parameter; the family is the one of the brms model.

- `link_function()` and `link_inverse()` return the link function of the
  first parameter that has one and its inverse. bmm models have one link
  per parameter, so for anything else use `model$links`.

## See also

`insight::find_formula()`, `insight::find_parameters()`,
`insight::get_priors()`

## Examples

``` r
if (FALSE) { # \dontrun{
insight::find_formula(fit)
parameters::model_parameters(fit)
bayestestR::describe_posterior(fit)
modelbased::estimate_means(fit, by = "set_size", predict = "c")
} # }
```
