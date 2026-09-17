# Get Default priors for Measurement Models specified in BMM

Obtain the default priors for a Bayesian multilevel measurement model,
as well as information for which parameters priors can be specified.
Given the `model`, the `data` and the `formula` for the model, this
function will return the default priors that would be used to estimate
the model. Additionally, it will return all model parameters that have
no prior specified (flat priors). This can help to get an idea about
which priors need to be specified and also know which priors were used
if no user-specified priors were passed to the
[`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md) function.

The default priors in `bmm` tend to be more informative than the default
priors in `brms`, as we use domain knowledge to specify the priors.

## Usage

``` r
# S3 method for class 'bmmformula'
default_prior(object, data, model, formula = object, ...)
```

## Arguments

- object:

  A `bmmformula` object

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

- formula:

  An object of class `bmmformula`. A symbolic description of the model
  to be fitted.

- ...:

  Further arguments passed to
  [`brms::default_prior()`](https://paulbuerkner.com/brms/reference/default_prior.html)

## Value

A data.frame with columns specifying the `prior`, the `class`, the
`coef` and `group` for each of the priors specified. Separate rows
contain the information on the parameters (or parameter classes) for
which priors can be specified.

## See also

[`supported_models()`](https://venpopov.com/bmm/dev/reference/supported_models.md),
[`brms::default_prior()`](https://paulbuerkner.com/brms/reference/default_prior.html)

## Examples

``` r
default_prior(bmf(c ~ 1, kappa ~ 1),
  data = oberauer_lin_2017,
  model = sdm(resp_error = "dev_rad")
)
#>                     prior     class coef group resp  dpar nlpar   lb   ub tag
#>     student_t(5, 2, 0.75) Intercept                     c       <NA> <NA>    
#>  student_t(5, 1.75, 0.75) Intercept                 kappa       <NA> <NA>    
#>               constant(0) Intercept                             <NA> <NA>    
#>  source
#>    user
#>    user
#>    user
```
