# Plot the latent decision-variable distributions of an SDT model

Draws the model-implied noise and signal evidence densities on the
latent decision axis. For
[`sdt_yn()`](https://popov-lab.github.io/bmm/dev/reference/sdt_yn.md)/[`sdt_rating()`](https://popov-lab.github.io/bmm/dev/reference/sdt_rating.md)
the response criterion or the K-1 confidence thresholds are added as
dashed vertical lines with a shaded credible band; when several
boundaries are shown together (e.g. base-rate criteria or the rating
thresholds) they are colour-coded.
[`sdt_mafc()`](https://popov-lab.github.io/bmm/dev/reference/sdt_mafc.md)
and
[`sdt_ranking()`](https://popov-lab.github.io/bmm/dev/reference/sdt_ranking.md)
have no boundary; with `show_competitors = TRUE` in
[`latent_sdt()`](https://popov-lab.github.io/bmm/dev/reference/latent_sdt.md)
the max-of-distractors densities are overlaid as dashed lines, one per
set size. Distinct `d`/`sdratio` conditions are faceted.

## Usage

``` r
# S3 method for class 'bmm_sdt_latent'
plot(x, condition_col = NULL, line_alpha = 0.9, ...)
```

## Arguments

- x:

  A `"bmm_sdt_latent"` object from
  [`latent_sdt()`](https://popov-lab.github.io/bmm/dev/reference/latent_sdt.md).

- condition_col:

  Optional character. Condition column(s) for faceting. If `NULL`
  (default), all condition columns are used.

- line_alpha:

  Numeric. Opacity of the criterion/threshold lines (default `0.9`).

- ...:

  Ignored.

## Value

A `ggplot2` object.

## See also

[`latent_sdt()`](https://popov-lab.github.io/bmm/dev/reference/latent_sdt.md),
[`roc_sdt()`](https://popov-lab.github.io/bmm/dev/reference/roc_sdt.md)
