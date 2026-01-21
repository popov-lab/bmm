# Restructure Old `bmmfit` Objects

Restructure old `bmmfit` objects to work with the latest bmm version.
This function is called internally when applying post-processing
methods.

## Usage

``` r
# S3 method for class 'bmmfit'
restructure(x, ...)
```

## Arguments

- x:

  An object of class `bmmfit`.

- ...:

  Currently ignored.

## Value

A `bmmfit` object compatible with the latest version of bmm and brms.

## Examples

``` r
if (FALSE) { # file.exists("bmmfit_old.rds")
# Load an old bmmfit object
old_fit <- readRDS("bmmfit_old.rds")
new_fit <- restructure(old_fit)
}
```
