# Generic S3 method for configuring the default prior for a bmmodel

Called by bmm() to automatically construct the priors for a given model,
data and formula, and combine it with the prior given by the user. The
first method executed is configure_prior.bmmodel, which will build the
prior based on information from the model object such as
fixed_parameters, default_priors, etc. Thus it is important to define
these values in the model object. The function will also recognize if
the user has specified that some parameters should be fixed to a
constant and put the appropriate constant priors. Any additional priors
that a developer wants to specify, which are not based on information in
the model object, can be defined in the configure_prior.\* method for
the model. See configure_prior.imm_full for an example.

## Usage

``` r
configure_prior(model, data, formula, user_prior, ...)
```

## Arguments

- model:

  A `bmmodel` object

- data:

  A data.frame containing the data used in the model

- formula:

  A `brmsformula` object returned from configure_model()

- user_prior:

  A `brmsprior` object given by the user as an argument to bmm()

- ...:

  Additional arguments passed to the method

## Value

A `brmsprior` object containing the default priors for the model

## Examples

``` r
if (FALSE) { # isTRUE(Sys.getenv("BMM_EXAMPLES"))
configure_prior.mixture3p <- function(model, data, formula, user_prior, ...) {
  # if there is set_size 1 in the data, set constant prior over thetant for set_size1
  prior <- brms::empty_prior()
  set_size_var <- model$other_vars$set_size
  prior_cond <- any(data$ss_numeric == 1) && !is.numeric(data[[set_size_var]])

  thetant_preds <- rhs_vars(formula$pforms$thetant)
  if (prior_cond && set_size_var %in% thetant_preds) {
    prior <- prior + brms::prior_("constant(-100)",
      class = "b",
      coef = paste0(set_size_var, 1),
      nlpar = "thetant"
    )
  }
  # check if there is a random effect on theetant that include set_size as predictor
  bterms <- brms::brmsterms(formula$pforms$thetant)
  re_terms <- bterms$dpars$mu$re
  if (!is.null(re_terms)) {
    for (i in 1:nrow(re_terms)) {
      group <- re_terms$group[[i]]
      form <- re_terms$form[[i]]
      thetant_preds <- rhs_vars(form)

      if (prior_cond && set_size_var %in% thetant_preds) {
        prior <- prior + brms::prior_("constant(1e-8)",
          class = "sd",
          coef = paste0(set_size_var, 1),
          group = group,
          nlpar = "thetant"
        )
      }
    }
  }

  prior
}
}
```
