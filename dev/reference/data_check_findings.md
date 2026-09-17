# Generic S3 method for model-specific pre-fit data findings

Called by
[`bmm_data_check()`](https://venpopov.com/bmm/dev/reference/bmm_data_check.md)
to collect model-specific diagnostics of common data coding mistakes.
Like
[`check_data()`](https://venpopov.com/bmm/dev/reference/check_data.md),
methods dispatch over the model class chain from general to specific,
and every method should concatenate its own findings with
[`NextMethod()`](https://rdrr.io/r/base/UseMethod.html). Each finding is
a list with elements `severity` ("warning" or "note") and `message`.
Unlike
[`check_data()`](https://venpopov.com/bmm/dev/reference/check_data.md),
these diagnostics never throw - they only describe likely problems, so a
report can always be produced.

Defining a method is entirely optional: the default method returns an
empty list, and
[`bmm_data_check()`](https://venpopov.com/bmm/dev/reference/bmm_data_check.md)
produces its full generic report for models without any method. Add one
only when a model has a common data mistake that the hard checks in
[`check_data()`](https://venpopov.com/bmm/dev/reference/check_data.md)
deliberately tolerate - a mistake
[`check_data()`](https://venpopov.com/bmm/dev/reference/check_data.md)
already warns about needs no method, because
[`bmm_data_check()`](https://venpopov.com/bmm/dev/reference/bmm_data_check.md)
captures that warning and reports it under "Hard checks".

## Usage

``` r
data_check_findings(model, data, formula)
```

## Arguments

- model:

  A `bmmodel` object

- data:

  The user supplied data.frame, before any transformations by
  [`check_data()`](https://venpopov.com/bmm/dev/reference/check_data.md)

- formula:

  The user supplied `bmmformula`

## Value

A list of findings (possibly empty), each built with
[`data_check_finding()`](https://venpopov.com/bmm/dev/reference/data_check_finding.md)

## See also

[`data_check_finding()`](https://venpopov.com/bmm/dev/reference/data_check_finding.md),
[`bmm_data_check()`](https://venpopov.com/bmm/dev/reference/bmm_data_check.md)

## Examples

``` r
# a method for a model whose check_data() tolerates unnormalized weights
data_check_findings.my_model <- function(model, data, formula) {
  weights <- data[[model$other_vars$weights]]
  findings <- list()
  if (!isTRUE(all.equal(sum(weights), 1))) {
    findings <- list(data_check_finding(
      "warning",
      glue::glue("The weights sum to {sum(weights)} rather than 1.")
    ))
  }
  c(findings, NextMethod())
}
```
