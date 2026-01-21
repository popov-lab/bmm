# Generic S3 method for reverting any postprocessing of the fitted brm model

Called by update.bmmfit() to automatically revert some of the
postprocessing depending on the model type. It will call the appropriate
revert_postprocess_brm.\* methods based on the list of classes defined
in the .model\_\* functions. For models with several classes listed, it
will call the functions in the order they are listed. For example, for
the sdm model, the postprocessing involves setting the link function for
the c parameter to "log", because it was coded manually in the stan
code, but it was specified as "identity" in the brms custom family.
However, during the update process, the link function should be set back
to "identity". Only use this if you have a specific reason to revert the
postprocessing (if otherwise the update method would produce incorrect
results).

## Usage

``` r
revert_postprocess_brm(model, fit, ...)
```

## Arguments

- model:

  A model list object returned from check_model()

- fit:

  the fitted brm model

- ...:

  Additional arguments passed to the method

## Value

An object of class brmsfit, with any necessary postprocessing applied

## Examples

``` r
if (FALSE) { # isTRUE(Sys.getenv("BMM_EXAMPLES"))
fit <- readRDS("my_saved_fit.rds")
postprocessed_fit <- prostprocess_brm(fit)
reverted_fit <- revert_postprocess_brm(postprocessed_fit)
}
```
