# Convert `bmmformula` objects to `brmsformula` objects

Called by
[`configure_model()`](https://venpopov.com/bmm/dev/reference/configure_model.md)
inside [`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md) to
convert the `bmmformula` into a `brmsformula` based on information in
the model object. It will call the appropriate bmf2bf.\\ methods based
on the classes defined in the model\_\\ function.

## Usage

``` r
bmf2bf(model, formula)
```

## Arguments

- model:

  The model object defining one of the supported \`bmmodels“

- formula:

  The `bmmformula` that should be converted to a `brmsformula`

## Value

A `brmsformula` defining the response variables and the additional
parameter formulas for the specified `bmmodel`

## Examples

``` r
model <- mixture2p(resp_error = "error")

formula <- bmmformula(
  thetat ~ 0 + set_size + (0 + set_size | id),
  kappa ~ 1 + (1 | id)
)

brms_formula <- bmf2bf(model, formula)
```
