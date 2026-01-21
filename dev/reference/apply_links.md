# Apply link functions for parameters in a formula or bmmformula

This function applies the specified link functions in the list of
`links` to the `formula` or `bmmformula` that is passed to it. This
function is mostly used internally for configuring `bmmodels`.

## Usage

``` r
apply_links(formula, links = nlist())
```

## Arguments

- formula:

  A `formula` or `bmmformula` that the links should be applied to

- links:

  A list of `links` that should be applied to the formula. Each element
  in this list should be named using the parameter labels the links
  should be applied for and contain a character variable specifying the
  link to be applied. Currently implemented links are: "log", "logit",
  "probit", and "identity".

## Value

The `formula` or `bmmformula` the links have been applied to

## Examples

``` r
# specify a bmmformula
form <- bmf(x ~ a + c, kappa ~ 1, a ~ 1, c ~ 1)
links <- list(a = "log", c = "logit")

apply_links(form, links)
#> x ~ exp(a) + inv_logit(c)
#> kappa ~ 1
#> a ~ 1
#> c ~ 1
```
