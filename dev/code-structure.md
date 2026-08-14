# NA

## Main user entry point - `bmm()`

The main function for fitting models is
[`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md). This function
is the main entry point for users to fit models. It is set-up to be
independent of the specific models that are implemented in the package.

``` r

bmm <- function(formula, data, model,
                prior = NULL,
                sort_data = getOption('bmm.sort_data', "check"),
                silent = getOption('bmm.silent', 1),
                backend = getOption('brms.backend', NULL),
                file = NULL, file_compress = TRUE,
                file_refit = getOption('bmm.file_refit', FALSE), ...) {
  deprecated_args(...)
  dots <- list(...)

  x <- read_bmmfit(file, file_refit)
  if (!is.null(x)) return(x)

  # set temporary global options and return modified arguments for brms
  configure_opts <- nlist(sort_data, silent, backend, parallel = dots$parallel,
                          cores = dots$cores)
  opts <- configure_options(configure_opts)
  dots$parallel <- NULL

  # check model, formula and data, and transform data if necessary
  user_formula <- formula
  model <- check_model(model, data, formula)
  data <- check_data(model, data, formula)
  formula <- check_formula(model, data, formula)

  # generate the model specification to pass to brms later
  config_args <- configure_model(model, data, formula)

  # configure the default prior and combine with user-specified prior
  prior <- configure_prior(model, data, config_args$formula, prior)

  # estimate the model
  fit_args <- combine_args(nlist(config_args, opts, dots, prior))
  fit <- call_brm(fit_args)

  # model post-processing
  fit <- postprocess_brm(model, fit, fit_args = fit_args, user_formula = user_formula,
                         configure_opts = configure_opts)

  # save the fitted model object if !is.null
  save_bmmfit(fit, file, compress = file_compress)
}
```

It calls several subroutines, implemented as generic S3 methods, to:

- `configure_options()` - to configure local options for fitting, such
  as parallel sampling,
- [`check_model()`](https://venpopov.com/bmm/dev/reference/check_model.md) -
  check if the model exists
- [`check_formula()`](https://venpopov.com/bmm/dev/reference/check_formula.md) -
  check if the formula is specified correctly and transform it to a
  brmsformula
- [`check_data()`](https://venpopov.com/bmm/dev/reference/check_data.md) -
  check whether the data contains all necessary information
- [`configure_model()`](https://venpopov.com/bmm/dev/reference/configure_model.md) -
  configures the model called for fitting
- [`configure_prior()`](https://venpopov.com/bmm/dev/reference/configure_prior.md) -
  sets the default priors for the model and combines them with the user
  prior
- `call_brm()` - fit the model using the `brm()` function from the
  `brms` package
- [`postprocess_brm()`](https://venpopov.com/bmm/dev/reference/postprocess_brm.md) -
  to post-process the fitted model

In addition, it also tests if the specified `bmmodel` has already been
estimated and saved to a `file`. This is done via the `read_bmmfit` and
`save_bmmfit` functions.

## Models

All models in the package are defined as S3 classes and follow a strict
template. This allows us to implement general methods for handling model
fitting, data checking, and post-processing. Each model has an internal
function that defines the model and its parameters, and a user-facing
alias. For a complete example model file and an explanation, see
\[Section @sec-example-model\]. The general model template looks like
this:

``` r

.model_my_new_model <- function(resp_var1 = NULL, required_args1 = NULL, 
                                required_arg2 = NULL, links = NULL, version = NULL,
                                call = NULL, ...) {
  out <- structure(
    list(
      resp_vars = nlist(resp_error),
      other_vars = nlist(),
      domain = "",
      task = "",
      name = "",
      version = "",
      citation = "",
      requirements = "",
      parameters = list(),
      links = list(),
      fixed_parameters = list(),
      default_priors = list(),
      version = version,
      void_mu = FALSE
    ),
    class = c("bmmodel", "my_new_model"),
    call = call
  )
  out$links[names(links)] <- links
  out
}
```

Each model is accompanied by a user-facing alias, the documentation of
which is generated automatically based on the info list in the model
definition.

``` r

# user facing alias
# information in the title and details sections will be filled in
# automatically based on the information in the .model_modelname()$info
#' @title `r .model_my_new_model()name`
#' @name Model Name#' @details `r model_info(.model_my_new_model())`
#' @param resp_var1 A description of the response variable
#' @param required_arg1 A description of the required argument
#' @param required_arg2 A description of the required argument
#' @param ... used internally for testing, ignore it
#' @return An object of class `bmmmodel`
#' @export
#' @examples
#' \dontrun{
#' # put a full example here (see 'R/bmm_model_mixture3p.R' for an example)
#' }
my_new_model <- function(resp_var1, required_arg1, required_arg2, 
                         links = NULL, version = NULL, ...) {
  call <- match.call()
  stop_missing_args()
  .model_my_new_model(resp_var1 = resp_var1, required_arg1 = required_arg1,
                      required_arg2 = required_arg2, links = links, version = version,
                      call = call, ...)
}
```

Then users can fit the model using the
[`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md) function, and
the model will be automatically recognized and handled by the package:

``` r

fit <- bmm(formula = my_bmmformula, 
           data = my_data, 
           model = my_new_model(resp_var1, required_arg1, required_arg2))
```

## S3 methods

The package uses S3 methods to handle different models. This means that
the same function can behave differently depending on the class of the
object it is called with. For example, the `configure_model(model)`
function called by
[`fit_model()`](https://venpopov.com/bmm/dev/reference/bmm.md), is
generally defined as:

``` r

configure_model <- function(model) {
   UseMethod('configure_model')
}
```

and it will call a function `configure_model.modelname()` that is
specified for each model. The same is true for other functions, such as
[`check_data()`](https://venpopov.com/bmm/dev/reference/check_data.md),
[`postprocess_brm()`](https://venpopov.com/bmm/dev/reference/postprocess_brm.md),
and
[`check_formula()`](https://venpopov.com/bmm/dev/reference/check_formula.md).
This allows us to add new models without having to edit the main fitting
function, [`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md).

## File organization

The `bmm` package is organized into several files. The main files are:

### `R/bmm.R`

It contains the main function for fitting models,
[`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md). This function
is the main entry point for users to fit models. It is set-up to be
independent of the specific models that are implemented in the package.

To add new models, you do not have to edit this file. The functions
above are generic S3 methods, and they will automatically recognize new
models if you add appropriate methods for them (see section [Adding new
models](#adding-new-models-to-bmm)).

### `R/helpers-*.R`

`R/helpers-data.R`, `R/helpers-parameters.R`, `R/helpers-postprocess.R`,
`R/helpers-model.R`, and `R/helpers-prior.R`

These files define the main generic S3 methods for checking data,
post-processing the fitted model, configuring the model, checking the
model formula, and combining priors. They contain the default methods
for these functions, which are called by
[`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md) if no specific
method is defined for a model. If you want to add a new model, you will
need to add specific methods for these functions for your model. *You do
not need to edit these files to add a new model.*

### `R/bmmformula.R`

This file contains the definition of the `bmmformula` class, which is
used to represent the formula for the model. It contains the
[`bmmformula()`](https://venpopov.com/bmm/dev/reference/bmmformula.md)
function and its alias
[`bmf()`](https://venpopov.com/bmm/dev/reference/bmmformula.md), which
is used to create a new formula object.

In addition, it contains the definition of the `bmf2bf` S3 method that
is used to convert a `bmmformula` object into a `brms_formula`\`object.
This is necessary, as `brmsformula` objects are required to include the
response variable in the first formula line. In contrast `bmmformula`
objects only contain formulas predicting the paramaeters of a
`bmmmodel`. The `bmf2bf` S3 method is used to perform this conversion
and add the first formula line and including the response variable in
the `brmsformula` created during model configuration.

### `R/model_*.R`

Each model and it’s methods is defined in a separate file. For example,
the 3-parameter mixture model is defined in `model_mixture3p.R`. This
file contains the internal function that defines the model and its
parameters, and the specific methods for the generic S3 functions. Your
new model will exist in a file like this. The name of the file should be
`model_name_of_your_model.R`. You don’t have to add this file manually -
see section [Adding new models](#adding-new-models-to-bmm).

### `R/distributions.R`

This file contains the definition of the custom distributions that are
used in the package. It specifies the density, random number generation,
and probability functions for the custom distributions. If your model
requires a custom distribution, you will need to add it to this file.
These are not used during model fitting, but can be used to generate
data from the model, and to plot the model fit.

### `R/utils.R`, `R/brms-misc.R`, `R/restructure.R`, `R/summary.R`, `R/update.R`

Various utility functions.

### `inst/stan_chunks/`

This directory contains the Stan chunks that are passed to the
[`brms::stanvar()`](https://paulbuerkner.com/brms/reference/stanvar.html)
function. These are used to define the custom distributions that are
used in the package. If you add a new custom distribution, you will need
to add a new Stan chunk to this directory. Each model has several files,
one for each corresponding stanvar block.
