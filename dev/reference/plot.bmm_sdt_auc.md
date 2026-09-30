# Plot the posterior AUC distribution from a fitted SDT model

Plot the posterior AUC distribution from a fitted SDT model

## Usage

``` r
# S3 method for class 'bmm_sdt_auc'
plot(x, condition_col = NULL, ...)
```

## Arguments

- x:

  A `"bmm_sdt_auc"` object from
  [`auc_sdt()`](https://popov-lab.github.io/bmm/dev/reference/auc_sdt.md).

- condition_col:

  Optional character. Condition column for colour. If `NULL` (default),
  auto-detected.

- ...:

  Ignored.

## Value

A `ggplot2` object.

## See also

[`auc_sdt()`](https://popov-lab.github.io/bmm/dev/reference/auc_sdt.md),
[`roc_sdt()`](https://popov-lab.github.io/bmm/dev/reference/roc_sdt.md)
