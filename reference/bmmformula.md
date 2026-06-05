# Create formula for predicting parameters of a `bmmodel`

This function is used to specify the formulas predicting the different
parameters of a `bmmodel`.

## Usage

``` r
bmmformula(...)

bmf(...)
```

## Arguments

- ...:

  Formulas for predicting a `bmmodel` parameter. Each formula for a
  parameter should be specified as a separate argument, separated by
  commas

## Value

A list of formulas for each parameters being predicted

## General formula structure

The formula argument accepts formulas of the following syntax:

    parameter ~ fixed_effects + (random_effects | grouping_variable)

`bmm` formulas are built on `brms` formulas and function in nearly the
same way, so you can use most of the `brms` formula syntax. The main
differences is that in `bmm` formulas, the response variable is not
specified in the formula. Instead, each parameter of the model is
explicitly specified as the left-hand side of the formula. In `brms`,
the response variable is always specified as the left-hand side of the
first formula, which implicitly means that any predictors in the first
formula are predictors of the `mu` parameter of the model. In general,
measurement models do not all have a `mu` parameter, therefore it is
more straigthforward to explicitely predict each parameter of the model.

For example, in the following `brms` formula for the drift diffusion
model, the first line corresponds to the drift rate parameter, but this
is not explicitely stated.

          brmsformula(rt | dec(response) ~ condition + (condition | id),
                      bs ~ 1 + (1 | id),
                      ndt ~ 1 + (1 | id),
                      bias ~ 1 + (1 | id))

In `bmm`, the same formula would be written as:

          bmmformula(drift ~ condition + (condition | id),
                     bs ~ 1 + (1 | id),
                     ndt ~ 1 + (1 | id),
                     bias ~ 1 + (1 | id))

and the rt and response variables would be specified in the model
argument of the [`bmm()`](https://venpopov.com/bmm/reference/bmm.md)
function.

Aside from that, the `bmm` formula syntax is the same as the `brms`
formula syntax. For more information on the `brms` formula syntax, see
[`brms::brmsformula()`](https://paulbuerkner.com/brms/reference/brmsformula.html).

You can also use the `bmf()` function as a shorthand for `bmmformula()`.

You can also set parameters to a constant value by using par = value
syntax:

          bmf(drift ~ condition + (condition | id),
              bs ~ 1 + (1 | id),
              ndt ~ 1 + (1 | id),
              bias = 0.5)

in which case the bias parameter will not be estimated but it will be
fixed to 0.5

## Examples

``` r
imm_formula <- bmmformula(
  c ~ 0 + set_size + (0 + set_size | id),
  a ~ 1,
  kappa ~ 0 + set_size + (0 + set_size | id)
)

# or use the shorter alias 'bmf'
imm_formula2 <- bmf(
  c ~ 0 + set_size + (0 + set_size | id),
  a ~ 1,
  kappa ~ 0 + set_size + (0 + set_size | id)
)
identical(imm_formula, imm_formula2)
#> [1] TRUE
```
