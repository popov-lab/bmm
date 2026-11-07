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
the model, as `mpt` does; the example shows the form such a method
takes.

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
if (FALSE) { # identical(Sys.getenv("BMM_EXAMPLES"), "true")
configure_prior.mymodel <- function(model, data, formula, user_prior, ...) {
  # a prior that depends on the data, which the model object cannot declare
  if (!"condition" %in% names(data)) {
    return(NULL)
  }
  brms::prior_("normal(0, 0.5)", class = "b", dpar = "kappa")
}
}
```
