# Generic S3 method for model-specific lines in `print()` of a model

Called by [`print()`](https://rdrr.io/r/base/print.html) on a `bmmodel`
after the lines all models share and before the closing pointer to
[`parameter_info()`](https://venpopov.com/bmm/reference/parameter_info.md).
The default method prints nothing.

## Usage

``` r
print_model_details(model, ...)
```

## Arguments

- model:

  A `bmmodel` object

- ...:

  Unused

## Value

Called for its printed output; returns `NULL` invisibly
