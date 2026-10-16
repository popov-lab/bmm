# Generic S3 method for configuring the model to be fit by brms

Called by bmm() to automatically construct the model formula, with its
family, and any stanvars for the model specified by the user. Default
priors are not built here but by
[`configure_prior()`](https://venpopov.com/bmm/dev/reference/configure_prior.md)
from the model's `default_priors`; a model that needs priors beyond
those defines its own `configure_prior.<model>()` method (e.g.
`mixture3p`, `mpt`). It will call the appropriate configure_model.\*
functions based on the list of classes defined in the .model\_\*
functions. Each model has one method, for its model class, or one per
version class when the versions are configured differently (e.g. `imm`,
`ezdm`, `cswald`).

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
  # retrieve arguments from the data check
  max_set_size <- attr(data, "max_set_size")
  lure_idx <- attr(data, "lure_idx_vars")
  nt_features <- model$other_vars$nt_features
  set_size_var <- model$other_vars$set_size

  # construct initial brms formula
  formula <- bmf2bf(model, formula) +
    brms::lf(kappa2 ~ 1) +
    brms::lf(mu2 ~ 1) +
    brms::nlf(theta1 ~ thetat) +
    brms::nlf(kappa1 ~ kappa)

  # additional internal terms for the mixture model formula
  kappa_nts <- paste0("kappa", 3:(max_set_size + 1))
  theta_nts <- paste0("theta", 3:(max_set_size + 1))
  mu_nts <- paste0("mu", 3:(max_set_size + 1))

  for (i in 1:(max_set_size - 1)) {
    formula <- formula +
      glue_nlf("{kappa_nts[i]} ~ kappa") +
      glue_nlf(
        "{theta_nts[i]} ~ {lure_idx[i]} * (thetant + log(inv_ss)) + ",
        "(1 - {lure_idx[i]}) * (-100)"
      ) +
      glue_nlf("{mu_nts[i]} ~ {nt_features[i]}")
  }

  # define mixture family
  vm_list <- lapply(1:(max_set_size + 1), function(x) brms::von_mises(link = "identity"))
  vm_list$order <- "none"
  formula$family <- brms::do_call(brms::mixture, vm_list)

  nlist(formula, data)
}
}
```
