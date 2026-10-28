# Sample more chains for a fitted model

`add_chains()` samples more chains of a model fitted with
[`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md) and returns one
fit with the chains of `fit` followed by the new ones. Use it when a fit
needs more draws, for example because its effective sample sizes are too
small, and the chains it has should be kept.

## Usage

``` r
add_chains(fit, chains = 1, seed = NULL, ...)
```

## Arguments

- fit:

  A `bmmfit` object sampled with MCMC (`algorithm = "sampling"`) by the
  rstan or cmdstanr backend.

- chains:

  The number of chains to add, a positive whole number.

- seed:

  The seed of the new chains. If `NULL` (the default), a seed is drawn
  from R's random number generator. bmm draws the default initial values
  from R's random number generator too, so call
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before
  `add_chains()` to make the new chains reproducible, with or without
  `seed`. A seed that a chain of `fit` was sampled with is refused,
  because the new chains could repeat that chain's draws.

- ...:

  Further arguments passed to
  [`update.bmmfit()`](https://venpopov.com/bmm/dev/reference/update.bmmfit.md),
  such as `cores`, `refresh` or `control`. `threads` works only for a
  fit that was sampled with threading, because a threaded Stan program
  differs from an unthreaded one. Arguments that would change the model,
  the data, the draws per chain or what is sampled (`formula.`,
  `newdata`, `data`, `data2`, `model`, `prior`, `sample_prior`, `iter`,
  `warmup`, `thin`, `algorithm`, `backend` and `normalize`) are refused:
  use [`update()`](https://rdrr.io/r/stats/update.html) for those and
  sample all chains anew. `file` and `file_compress` are refused too;
  save the result with
  [`saveRDS()`](https://rdrr.io/r/base/readRDS.html).

## Value

A `bmmfit` object with the chains of `fit` followed by the new ones,
combined with
[`brms::combine_models()`](https://paulbuerkner.com/brms/reference/combine_models.html).
Everything else is taken from `fit`, except criteria added with
[`brms::add_criterion()`](https://paulbuerkner.com/brms/reference/add_criterion.html):
they were computed from the chains of `fit` alone and are removed, so
add them again. The result has no `file`, so
[`brms::add_criterion()`](https://paulbuerkner.com/brms/reference/add_criterion.html)
does not write it over the file of `fit`. `version` and the other fields
that describe the fit are those of `fit`; the new chains ran the same
Stan program.
[`fit_info()`](https://venpopov.com/bmm/dev/reference/fit_info.md)
reports the time, step size and inverse metric of every chain.

## Details

The new chains are sampled with
[`update()`](https://rdrr.io/r/stats/update.html), which recompiles the
Stan model only when it needs to, for example for a cmdstanr fit loaded
in a new R session. If the Stan program that bmm now writes for the
model differs from the one `fit` was sampled with, as it can after an
update of bmm or brms, the new chains would not sample the same
posterior, and `add_chains()` stops after sampling them. Refit the model
with [`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md) in that
case.

## See also

[`fit_info()`](https://venpopov.com/bmm/dev/reference/fit_info.md),
[`update.bmmfit()`](https://venpopov.com/bmm/dev/reference/update.bmmfit.md)

## Examples

``` r
if (FALSE) { # identical(Sys.getenv("BMM_EXAMPLES"), "true")
fit <- bmm(
  formula = bmmformula(c ~ 1, kappa ~ 1),
  data = data.frame(y = rsdm(1000)),
  model = sdm(resp_error = "y"),
  chains = 2
)

# two more chains, four in all
fit <- add_chains(fit, chains = 2)
fit_info(fit, "time")
}
```
