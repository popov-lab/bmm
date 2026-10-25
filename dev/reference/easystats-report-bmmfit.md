# report and performance support for bmmfit objects

These S3 methods let `report::report()` describe a bmm fit as the model
you specified with
[`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md), and let
`performance::model_performance()` return only the indices that are
defined for bmm models. They build on the insight methods described in
[easystats-bmmfit](https://venpopov.com/bmm/dev/reference/easystats-bmmfit.md).

## Usage

``` r
# S3 method for class 'bmmfit'
format_model(x)

# S3 method for class 'bmmfit'
report_model(x, table = NULL, ...)

# S3 method for class 'bmmfit'
report_random(x, ...)

# S3 method for class 'bmmfit'
report_priors(x, ...)

# S3 method for class 'bmmfit'
report_info(x, parameters = NULL, ...)

# S3 method for class 'bmmfit'
model_performance(model, metrics = c("LOOIC", "WAIC"), verbose = TRUE, ...)
```

## Arguments

- x, model:

  A `bmmfit` object (created by
  [`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md)).

- table:

  Ignored, for consistency with the generic.

- ...:

  Passed on to the methods of report and performance for brms fits.

- parameters:

  The result of `report::report_parameters()` for `x`, or `NULL` to
  compute it.

- metrics:

  The indices to return: `"LOOIC"`, `"WAIC"`, or `"all"` and `"common"`
  for both.

- verbose:

  Should a message name the requested indices that are left out?

## Value

See the generics, e.g. `report::report_model()` and
`performance::model_performance()`.

## Details

- `format_model()` names the model and its version, e.g. "Bayesian
  Signal Discrimination Model (SDM) by Oberauer (2023) (simple)".

- `report_model()` describes the model in one sentence, with all
  response variables of models such as
  [`ezdm()`](https://venpopov.com/bmm/dev/reference/ezdm.md) or
  [`m3()`](https://venpopov.com/bmm/dev/reference/m3.md) in one list,
  the formula of each model parameter as
  [`summary()`](https://rdrr.io/r/base/summary.html) prints it, and the
  grouping factors of the random effects.

- `report_priors()` gives the text of
  `prior_info(fit, format = "text")`: the priors in the names of the
  model parameters, with a note that they apply on the sampling (link)
  scale. The short version of the report, `summary(report(fit))`, leaves
  the priors out.

- `report_info()` keeps the paragraph of report that explains the
  parameter summaries, and replaces its sentence on convergence (the
  thresholds for R-hat and ESS) with the sentence of
  `report_convergence(fit, format = "text")`: the largest R-hat, the
  smallest bulk and tail ESS and the number of divergent transitions,
  without verdicts.

- `model_performance()` returns LOOIC and WAIC (with ELPD and the
  standard errors). R2, RMSE and the other indices of performance need
  the expected response or residuals, which bmm models do not define;
  requested ones are left out with a message, and when none of the
  requested indices is defined, LOOIC and WAIC are returned.

The parameter summaries, with their R-hat and ESS, come from report
unchanged. report lists no parameter estimates for a model without
predictors, such as an
[`mpt()`](https://venpopov.com/bmm/dev/reference/mpt.md) or
[`ezdm()`](https://venpopov.com/bmm/dev/reference/ezdm.md) fit with
intercepts only; use [`summary()`](https://rdrr.io/r/base/summary.html)
for those.
[`report_methods()`](https://venpopov.com/bmm/dev/reference/report_methods.md)
drafts a methods paragraph from the same facts in bmm's own words.

## See also

[`prior_info()`](https://venpopov.com/bmm/dev/reference/prior_info.md),
[`report_methods()`](https://venpopov.com/bmm/dev/reference/report_methods.md),
[`report_convergence()`](https://venpopov.com/bmm/dev/reference/report_convergence.md),
[easystats-bmmfit](https://venpopov.com/bmm/dev/reference/easystats-bmmfit.md)

## Examples

``` r
if (FALSE) { # \dontrun{
report::report(fit)
performance::model_performance(fit)
} # }
```
