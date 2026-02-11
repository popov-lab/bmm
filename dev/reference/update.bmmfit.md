# Update a bmm model

Update an existing bmm mode. This function calls
[`brms::update.brmsfit()`](https://paulbuerkner.com/brms/reference/update.brmsfit.html),
but it applies the necessary bmm postprocessing to the model object
before and after the update.

## Usage

``` r
# S3 method for class 'bmmfit'
update(object, formula., newdata = NULL, recompile = NULL, ...)
```

## Arguments

- object:

  An object of class `bmmfit`

- formula.:

  A
  [`bmmformula()`](https://venpopov.com/bmm/dev/reference/bmmformula.md).
  If missing, the original formula is used. Currently you have to
  specify a full `bmmformula`

- newdata:

  An optional data frame containing the variables in the model

- recompile:

  Logical, indicating whether the Stan model should be recompiled. If
  NULL (the default), update tries to figure out internally, if
  recompilation is necessary. Setting it to FALSE will cause all Stan
  code changing arguments to be ignored.

- ...:

  Further arguments passed to
  [`brms::update.brmsfit()`](https://paulbuerkner.com/brms/reference/update.brmsfit.html)

## Value

An updated `bmmfit` object refit to the new data and/or formula

## Details

When updating a brmsfit created with the cmdstanr backend in a different
R session, a recompilation will be triggered because by default,
cmdstanr writes the model executable to a temporary directory. To avoid
that, set option "cmdstanr_write_stan_file_dir" to a nontemporary path
of your choice before creating the original bmmfit.

For more information and examples, see
[`brms::update.brmsfit()`](https://paulbuerkner.com/brms/reference/update.brmsfit.html)

## Examples

``` r
if (FALSE) { # isTRUE(Sys.getenv("BMM_EXAMPLES"))
# generate artificial data from the Signal Discrimination Model
# generate artificial data from the Signal Discrimination Model
dat <- data.frame(y = rsdm(2000))

# define formula
ff <- bmf(c ~ 1, kappa ~ 1)

# fit the model
fit <- bmm(
  formula = ff,
  data = dat,
  model = sdm(resp_error = "y"),
  cores = 4,
  backend = "cmdstanr"
)

# update the model
fit <- update(fit, newdata = data.frame(y = rsdm(2000, kappa = 5)))
}
```
