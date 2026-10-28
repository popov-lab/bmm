# Extract information from a brmsfit object

Extract information from a brmsfit object

## Usage

``` r
fit_info(fit, what)
```

## Arguments

- fit:

  A brmsfit object, or a list of brmsfit objects

- what:

  String. What to return:

  - "time" for the sampling time per chain

  - "time_mean" for the mean sampling time

  - "inv_metric" for the adapted inverse metric per chain

  - "step_size" for the adapted step size per chain

## Value

Depends on `what` and the class of `fit`. For `brmsfit` objects,
information about the single fit is returned. For `brmsfit_list`
objects, a list or data.frame with the information for each fit is
returned.

- "time": A data.frame with the warmup, sampling and total time per
  chain

- "time_mean": A named numeric vector with the mean warmup, sampling and
  total time across chains

- "inv_metric": A list with one numeric vector per chain (the diagonal),
  or one matrix per chain for `metric = "dense_e"`

- "step_size": A list with one number per chain

Every chain of the fit is reported, also in a fit combined from several
fits with
[`brms::combine_models()`](https://paulbuerkner.com/brms/reference/combine_models.html)
or extended with
[`add_chains()`](https://venpopov.com/bmm/dev/reference/add_chains.md).
`NULL` is returned for what the sampler did not record: anything for a
fit without a stanfit or without chains, "time" and "time_mean" for
variational and Laplace runs, "inv_metric" and "step_size" for any
algorithm other than NUTS, and "inv_metric" for `metric = "unit_e"`. In
a fit combined from several fits, a unit_e chain next to chains with a
metric has a `NULL` element in the "inv_metric" list. Values read from
Stan's printed adaptation output (rstan fits, and fits combined from
several fits) have six significant digits by default. A cmdstanr unit_e
chain carries no printed step size, so in a fit combined from several
fits its step size is that of its last draw.

## Examples

``` r
if (FALSE) { # identical(Sys.getenv("BMM_EXAMPLES"), "true")
fit <- bmm(
  formula = bmmformula(c ~ 1, kappa ~ 1),
  data = data.frame(y = rsdm(1000)),
  model = sdm(resp_error = "y")
)

fit_info(fit, "time")
}
```
