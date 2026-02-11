# Generic S3 method for checking data based on model type

Called by [`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md) to
automatically perform checks on the data depending on the model type. It
will call the appropriate check_data methods based on the list of
classes defined in the .model\_\* functions. For models with several
classes listed, it will call the functions in the order they are listed.
Thus, any operations that are common to a group of models should be
defined in the appropriate check_data.\* function, where \\ corresponds
to the shared class. For example, for the .model_imm_abc model, this
corresponds to the following order of check_data.\* functions:
check_data() -\> check_data.circular(), check_data.non_targets() the
output of the final function is returned to bmm().

## Usage

``` r
check_data(model, data, formula)
```

## Arguments

- model:

  A model list object returned from check_model()

- data:

  The user supplied data.frame containing the data to be checked

- formula:

  The user supplied formula

## Value

A data.frame with the same number of rows as the input data, but with
additional columns added as necessary, any necessary transformations
applied, and attributes added to the data.frame for later use. If you
need to reuse variables created by the check_data.\* functions in
subsequent stages (e.g. in configure_model()), you can store and access
them using the attr() function.

## Examples

``` r
data <- oberauer_lin_2017
model <- sdmSimple(resp_error = "dev_rad")
#> Warning: The function `sdmSimple()` is deprecated. Please use `sdm()` instead.
formula <- bmf(c ~ 1, kappa ~ 1)
checked_data <- check_data(model, data, formula)
```
