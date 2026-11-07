# Generic S3 method for configuring the model to be fit by brms

Called by bmm() to automatically construct the model formula, with its
family, and any stanvars for the model specified by the user. Default
priors are not built here but by
[`configure_prior()`](https://venpopov.com/bmm/dev/reference/configure_prior.md)
from the model's `default_priors`; a model that needs priors beyond
those defines its own `configure_prior.<model>()` method (e.g. `mpt`).
It will call the appropriate configure_model.\* functions based on the
list of classes defined in the .model\_\* functions. Each model has one
method, for its model class, or one per version class when the versions
are configured differently (e.g. `imm`, `ezdm`, `cswald`).

## Usage

``` r
configure_model(model, data, formula)
```

## Arguments

- model:

  A model list object returned from check_model()

- data:

  The user supplied data.frame containing the data to be checked

- formula:

  The user supplied formula

## Value

A named list built with `nlist()`, containing:

- formula: An object of class `brmsformula`, the constructed model
  formula with the brms family attached as `formula$family`

- data: the user supplied data.frame, preprocessed by check_data

- stanvars: (optional) An object of class `stanvars`, when the model
  adds Stan code (a custom family, or Stan functions called from the
  formula). See
  [`brms::custom_family()`](https://paulbuerkner.com/brms/reference/custom_family.html)
  and
  [`brms::stanvar()`](https://paulbuerkner.com/brms/reference/stanvar.html)
  for more details.

## Details

A bare bones configure_model.\* method should look like this:

    configure_model.newmodel <- function(model, data, formula) {

       # preprocessing - e.g. extract arguments from data check, construct new variables
       <preprocessing code>

       # construct the formula
       formula <- bmf2bf(model, formula)

       # attach the family to the formula
       formula$family <- <code for new family>

       # return the list (add stanvars when the model adds Stan code: a custom
       # family, or Stan functions called from the formula)
       nlist(formula, data)
    }

## Examples

``` r
if (FALSE) { # identical(Sys.getenv("BMM_EXAMPLES"), "true")
configure_model.mixture3p <- function(model, data, formula) {
  # the number of non-target covariates the family needs is fixed by the data
  n_nt <- attr(data, "max_set_size") - 1

  formula <- bmf2bf(model, formula)
  formula$family <- .circmix_custom_family(
    model,
    family = "mixture3p_simple",
    weight_parameters = c("thetat", "thetant"),
    vint = TRUE, n_vreal = n_nt,
    log_lik = log_lik_mixture3p_simple,
    posterior_predict = posterior_predict_mixture3p_simple,
    posterior_epred = posterior_epred_mixture3p_simple
  )

  nlist(
    formula, data,
    stanvars = .circmix_model_stanvars(
      model, formula$family, "mixture3p_functions.stan",
      vint = "ss", vreal = list(nt = n_nt)
    )
  )
}
}
```
