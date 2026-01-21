# Get Activation Functions for different M3 versions

This function generates the activation functions for different versions
of the Memory Measurement Model (m3) implemented in the `bmm` package.
If no `bmmodel` object is passed then it will print the available model
versions.

## Usage

``` r
construct_m3_act_funs(model = NULL, warnings = TRUE)
```

## Arguments

- model:

  A bmmodel object that specifies the M3 model for which the activation
  functions should be generated. If no model is passed the available M3
  versions will be printed to the console.

- warnings:

  Logical flag to indicate if information about the generated model
  formulas should be printed when the function is called.

## Value

A bmmformula object with the activation functions for the m3 version
specified in the model object. The activation functions use the names of
the response categories specified in the model object.

## Examples

``` r
model <- m3(
 resp_cats = c("correct","other", "npl"),
 num_options = c(1, 4, 5),
 version = "ss"
)

construct_m3_act_funs(model, warnings = FALSE)
#> correct ~ b + a + c
#> other ~ b + a
#> npl ~ b
```
