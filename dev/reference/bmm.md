# Fit Bayesian Measurement Models

Fit a Bayesian measurement model using **brms** as a backend interface
to Stan.

## Usage

``` r
bmm(
  formula,
  data,
  model,
  prior = NULL,
  sort_data = getOption("bmm.sort_data", "check"),
  silent = getOption("bmm.silent", 1),
  backend = getOption("brms.backend", NULL),
  file = NULL,
  file_compress = TRUE,
  file_refit = getOption("bmm.file_refit", FALSE),
  ...
)

fit_model(
  formula,
  data,
  model,
  prior = NULL,
  sort_data = getOption("bmm.sort_data", "check"),
  silent = getOption("bmm.silent", 1),
  backend = getOption("brms.backend", NULL),
  ...
)
```

## Arguments

- formula:

  An object of class `bmmformula`. A symbolic description of the model
  to be fitted.

- data:

  An object of class data.frame, containing data of all variables used
  in the model. The names of the variables must match the variable names
  passed to the `bmmodel` object for required argurments.

- model:

  A description of the model to be fitted. This is a call to a `bmmodel`
  such as
  [`mixture3p()`](https://venpopov.com/bmm/dev/reference/mixture3p.md)
  function. Every model function has a number of required arguments
  which need to be specified within the function call. Call
  [`supported_models()`](https://venpopov.com/bmm/dev/reference/supported_models.md)
  to see the list of supported models and their required arguments

- prior:

  One or more `brmsprior` objects created by
  [`brms::set_prior()`](https://paulbuerkner.com/brms/reference/set_prior.html)
  or related functions and combined using the c method or the +
  operator. See also
  [`default_prior()`](https://venpopov.com/bmm/dev/reference/default_prior.bmmformula.md)
  for more help. Not necessary for the default model fitting, but you
  can provide prior constraints to model parameters

- sort_data:

  Logical. If TRUE, the data will be sorted by the predictor variables
  for faster sampling. If FALSE, the data will not be sorted, but
  sampling will be slower. If "check" (the default), `bmm()` will check
  if the data is sorted, and ask you via a console prompt if it should
  be sorted. You can set the default value for this option using global
  `options(bmm.sort_data = TRUE/FALSE/"check"))` or via
  `bmm_options(sort_data)`

- silent:

  Verbosity level between 0 and 2. If 1 (the default), most of the
  informational messages of compiler and sampler are suppressed. If 2,
  even more messages are suppressed. The actual sampling progress is
  still printed. Set refresh = 0 to turn this off as well. If using
  backend = "rstan" you can also set open_progress = FALSE to prevent
  opening additional progress bars.

- backend:

  Character. The backend to use for fitting the model. Can be "rstan" or
  "cmdstanr". If NULL (the default), "cmdstanr" will be used if the
  cmdstanr package is installed, otherwise "rstan" will be used. You can
  set the default backend using global
  `options(brms.backend = "rstan"/"cmdstanr")`

- file:

  Either `NULL` or a character string. If a string, the fitted model
  object is saved via [saveRDS](https://rdrr.io/r/base/readRDS.html) in
  a file named after the string. The .rds extension is added
  automatically. If the file already exists, `bmm` will load and return
  the saved model object. Unless you specify the `file_refit` argument
  as well, the existing files won't be overwritten, you have to manually
  remove the file in order to refit and save the model under an existing
  file name. The file name is stored in the `bmmfit` object for later
  usage. If the directory of the file does not exist, it will be
  created.

- file_compress:

  Logical or a character string, specifying one of the compression
  algorithms supported by [saveRDS](https://rdrr.io/r/base/readRDS.html)
  when saving the fitted model object.

- file_refit:

  Logical or character string. Modifies when the fit stored via the
  `file` argument is re-used. Can be set globally for the current R
  session via the `"bmm.file_refit"` option (see
  [options](https://rdrr.io/r/base/options.html)). If `TRUE` or
  "always", the model is fitted again. If `FALSE` or "never" (the
  default), the model saved under the name specified in `file` will be
  re-used. Note that unlike in `brms`, there is no "on_change" option

- ...:

  Further arguments passed to
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html) or
  Stan. See the description of
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html) for
  more details

## Value

An object of class brmsfit which contains the posterior draws along with
many other useful information about the model. Use methods(class =
"brmsfit") for an overview on available methods.

## Supported Models

The following models are supported:

- cswald(rt, response, links, version)

- ddm(rt, response, links)

- ezdm(mean_rt, var_rt, n_upper, n_trials, links, version)

- imm(resp_error, nt_features, nt_distances, set_size, regex, version)

- m3(resp_cats, num_options, choice_rule, version)

- mixture2p(resp_error)

- mixture3p(resp_error, nt_features, set_size, regex)

- sdm(resp_error, version)

Type ?modelname to get information about a specific model, e.g. ?imm

## bmmformula syntax

see [this online
article](https://venpopov.com/bmm/articles/bmm_bmmformula.html) for a
detailed description of the syntax and how it differs from the syntax
for **brmsformula**

## Default priors, Stan code and Stan data

For more information about the default priors in **bmm** and about who
to extract the Stan code and data generated by bmm and \#' brms, see
[the online
article](https://venpopov.com/bmm/articles/bmm_extract_info.html).

## Miscellaneous

Type [`help(package=bmm)`](https://venpopov.github.io/bmm/reference) for
a full list of available help topics.

**fit_model()** is a deprecated alias for **bmm()**.

## References

Frischkorn, G. T., & Popov, V. (2023). A tutorial for estimating mixture
models for visual working memory tasks in brms: Introducing the Bayesian
Measurement Modeling (bmm) package for R.
https://doi.org/10.31234/osf.io/umt57

## See also

[`supported_models()`](https://venpopov.com/bmm/dev/reference/supported_models.md),
[`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html),
[default_prior()](https://venpopov.com/bmm/dev/reference/default_prior.bmmformula.md),
[`bmmformula()`](https://venpopov.com/bmm/dev/reference/bmmformula.md),
[stancode()](https://venpopov.com/bmm/dev/reference/stancode.bmmformula.md),
[standata()](https://venpopov.com/bmm/dev/reference/standata.bmmformula.md)

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
}
```
