# Stan data of a fitted `bmm` model

The prediction methods of `brms`, such as
[`brms::log_lik()`](https://mc-stan.org/rstantools/reference/log_lik.html),
[`brms::posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html),
[`brms::posterior_epred()`](https://mc-stan.org/rstantools/reference/posterior_epred.html),
[`brms::predict.brmsfit()`](https://paulbuerkner.com/brms/reference/predict.brmsfit.html),
[`brms::fitted.brmsfit()`](https://paulbuerkner.com/brms/reference/fitted.brmsfit.html)
and [`brms::loo()`](https://mc-stan.org/loo/reference/loo.html), turn
`newdata` into Stan data through
[`standata()`](https://venpopov.com/bmm/dev/reference/standata.bmmformula.md).
For a `bmmfit`, `newdata` can be in the layout of the data you gave to
[`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md): the columns
`bmm` built from those data when it fitted the model are built again.
Rows of the fit's model frame, `fit$data`, are used as they are.

## Usage

``` r
# S3 method for class 'bmmfit'
standata(object, newdata = NULL, ...)
```

## Arguments

- object:

  A `bmmfit` object

- newdata:

  An optional data frame, in the layout of the data given to
  [`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md) or a subset
  of the rows of `object$data`. Its rows keep their order. When
  `newdata` is in the layout given to
  [`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md), the columns
  `bmm` builds itself are built again, and a column of `newdata` with
  one of those names is replaced.

- ...:

  Further arguments passed to
  [`brms::standata()`](https://paulbuerkner.com/brms/reference/standata.html)

## Value

A named list of objects containing the data to evaluate the model with
Stan.

## See also

[`brms::standata()`](https://paulbuerkner.com/brms/reference/standata.html),
[`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md)
