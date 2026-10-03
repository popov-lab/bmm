# Get parameter information for a bmm model

Returns a data frame with information about the model parameters: their
descriptions, whether they are fixed (and to which value), and their
link functions. Use
[`prior_info()`](https://popov-lab.github.io/bmm/dev/reference/prior_info.md)
for the priors of a fitted model.

## Usage

``` r
parameter_info(x, ...)

# S3 method for class 'bmmodel'
parameter_info(x, formula = NULL, ...)

# S3 method for class 'bmmfit'
parameter_info(x, ...)
```

## Arguments

- x:

  A `bmmodel` object (e.g., `sdm(resp_error = "y")`) or a `bmmfit`
  object (a fitted model returned by
  [`bmm`](https://popov-lab.github.io/bmm/dev/reference/bmm.md)).

- ...:

  Additional arguments (currently unused).

- formula:

  An optional `bmmformula` object. Only relevant for M3 custom models,
  where additional parameters are discovered from the formula. Ignored
  for all other models.

## Value

A data frame of class `bmm_parameters` with one row per parameter and
columns: `parameter`, `description`, `fixed`, `value`, and `link`.

## Examples

``` r
# For an unfitted model
parameter_info(sdm(resp_error = "y"))
#> Model:  Signal Discrimination Model (SDM) by Oberauer (2023) 
#> 
#>  parameter description                                        fixed value
#>  mu        Location parameter of the SDM distribution (in ... yes   0    
#>  c         Memory strength parameter of the SDM distribution  no    --   
#>  kappa     Precision parameter of the SDM distribution        no    --   
#>  link    
#>  tan_half
#>  log     
#>  log     

# For an M3 model
parameter_info(m3(
  resp_cats = c("corr", "other", "npl"),
  num_options = c(1, 4, 5),
  version = "ss"
))
#> Model:  The Multinomial / Memory Measurement Model 
#> 
#>  parameter description                                        fixed value
#>  b         Background activation. Added to each response c... yes   0    
#>  c         Context activation. Added to the item cued to b... no    --   
#>  a         General activation. Added to all items that wer... no    --   
#>  link    
#>  identity
#>  identity
#>  identity
```
