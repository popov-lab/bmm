# Create a summary of a fitted model represented by a `bmmfit` object

Create a summary of a fitted model represented by a `bmmfit` object

## Usage

``` r
# S3 method for class 'bmmfit'
summary(
  object,
  priors = FALSE,
  prob = 0.95,
  robust = FALSE,
  mc_se = FALSE,
  ...,
  backend = "bmm"
)
```

## Arguments

- object:

  An object of class `brmsfit`.

- priors:

  Logical; Indicating if priors should be included in the summary.
  Default is `FALSE`.

- prob:

  A value between 0 and 1 indicating the desired probability to be
  covered by the uncertainty intervals. The default is 0.95.

- robust:

  If `FALSE` (the default) the mean is used as the measure of central
  tendency and the standard deviation as the measure of variability. If
  `TRUE`, the median and the median absolute deviation (MAD) are applied
  instead.

- mc_se:

  Logical; Indicating if the uncertainty in `Estimate` caused by the
  MCMC sampling should be shown in the summary. Defaults to `FALSE`.

- ...:

  Other potential arguments

- backend:

  Choose whether to display the *bmm* summary method (default), or to
  display the *brms* summary method.

## Value

A list of class `bmmsummary` containing the summary of the model's
parameters, the model formula, the model, and the data used to fit the
model.

## Note

You can turn off the color output by setting the option
options(bmm.color_summary = FALSE) or bmm_options(color_summary = FALSE)

## See also

[`summary.brmsfit`](https://paulbuerkner.com/brms/reference/summary.brmsfit.html)

## Examples

``` r
if (FALSE) { # isTRUE(Sys.getenv("BMM_EXAMPLES"))
# generate artificial data from the Signal Discrimination Model
dat <- data.frame(y = rsdm(2000))

# define formula
ff <- bmmformula(c ~ 1, kappa ~ 1)

# fit the model
fit <- bmm(
  formula = ff,
  data = dat,
  model = sdm(resp_error = "y"),
  cores = 4,
  backend = "cmdstanr"
)

# summary of the model
summary(fit)
}
```
