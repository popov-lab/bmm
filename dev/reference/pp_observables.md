# Declare the observables of a bmm model for posterior predictive checks

`pp_observables()` returns the model's observable declaration used by
[`pp_check.bmmfit()`](https://venpopov.com/bmm/dev/reference/pp_check.bmmfit.md)
when `resp_var` is specified, or `NULL` for models that delegate fully
to
[`brms::pp_check()`](https://mc-stan.org/bayesplot/reference/pp_check.html).
`pp_simulate()` draws all observables jointly from the posterior
predictive distribution.

## Usage

``` r
pp_observables(model)

pp_simulate(model, prep)
```

## Arguments

- model:

  A `bmmodel` object.

- prep:

  A `brmsprep` object from
  [`brms::prepare_predictions()`](https://paulbuerkner.com/brms/reference/prepare_predictions.html).

## Value

`pp_observables()` returns `NULL` for a model that delegates fully to
[`brms::pp_check()`](https://mc-stan.org/bayesplot/reference/pp_check.html),
or a list with elements `observed` (a named character vector mapping
observable names to brms standata slots) and `checks` (a named list of
check definitions, each with a `compute` closure, a `label` and a
default bayesplot `type`). `pp_simulate()` returns a named list of
`ndraws` x `nobs` matrices, one per simulated observable.

## Details

A `pp_observables()` method returns `list(observed, checks)`:

- `observed`: named character vector mapping observable names to slots
  of the brms standata (`"Y"`, `"vreal1"`, `"vint1"`, `"trials"`,
  `"dec"`). The observable mapped to `"Y"` is the default check.

- `checks`: named list of entries built by the internal
  `.pp_observable()` constructor. Each `compute` closure receives a
  named list keyed by `names(observed)` and must be elementwise, so the
  identical closure produces `y` from length-N vectors and `yrep` from
  ndraws x N matrices.

A `pp_simulate()` method returns a named list of ndraws x nobs matrices
drawn jointly, typically through the internal `.pp_simulate_joint()`
helper around the model's `r*()` function. Simulating observables
independently would break their joint distribution (e.g. rt and response
under the DDM). Names not in `observed` are ignored; declared
observables that are not simulated (design quantities such as trial
counts) are filled in from the data.

Register exactly one method per model at the most general class level
where the declaration is identical across versions.

These two generics are exported so that model methods defined outside
bmm can be registered against them, but they are an internal developer
interface documented for bmm's own model authors and carry no stability
guarantee across releases.
