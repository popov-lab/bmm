# Plot a model-implied SDT ROC curve

Shows the smooth model-implied ROC as a credible-band ribbon, with the
model-implied operating points overlaid as colour-coded markers with
crosshair error bars (uncertainty in both directions). For **binary**
multi-criteria fits the points are the several criterion levels; for
**rating** models they are the K-1 confidence thresholds (labelled
`t1`..`t(K-1)`). Pass `observed` to overlay empirical points from
[`roc_observed()`](https://popov-lab.github.io/bmm/dev/reference/roc_observed.md).

## Usage

``` r
# S3 method for class 'bmm_sdt_roc'
plot(
  x,
  observed = NULL,
  condition_col = NULL,
  add_diagonal = TRUE,
  scale = c("probability", "quantile", "z"),
  ribbon_alpha = 0.25,
  point_size = 2.5,
  ...
)
```

## Arguments

- x:

  A `"bmm_sdt_roc"` object from
  [`roc_sdt()`](https://popov-lab.github.io/bmm/dev/reference/roc_sdt.md).

- observed:

  Optional `"bmm_sdt_roc_observed"` object from
  [`roc_observed()`](https://popov-lab.github.io/bmm/dev/reference/roc_observed.md)
  to overlay as empirical points.

- condition_col:

  Optional character. Condition column for colour/faceting. If `NULL`
  (default), auto-detected.

- add_diagonal:

  Logical. Draw the chance-level diagonal (default `TRUE`).

- scale:

  Either `"probability"` (default) for the usual hit vs. false-
  alarm-rate axes, or `"quantile"` (alias `"z"`) for the
  quantile-transformed axes (the inverse CDF of the fitted noise
  distribution).

- ribbon_alpha:

  Numeric. Transparency of the credible band (default `0.25`).

- point_size:

  Numeric. Size of operating-point markers (default `2.5`).

- ...:

  Ignored.

## Value

A `ggplot2` object. The credible band (ribbon) is the pointwise
posterior interval of the hit rate at each criterion value, plotted at
that criterion's posterior-mean false-alarm rate; it does not include
uncertainty in the false-alarm rate.

## Details

With `scale = "quantile"` (or its alias `scale = "z"`) the rates are
read on the distribution's quantile axis via the transform
`-qf(1 - rate)`, where `qf` is the inverse CDF of the fitted noise
distribution, and the axis label names that transform (`z` for normal,
`logit` for logistic, `loglog`/`cloglog` for the Gumbel distributions).
On this axis the binary model ROC is a straight line with slope
`1 / exp(sdratio)` and intercept
`d * sqrt((1 + exp(sdratio)^2) / 2) / exp(sdratio)` – the separation in
noise-SD units over the signal SD, since `d` is \\d_a\\ – so a slope
below 1 is the unequal-variance signature (signal SD \> noise SD), and
departures of the observed points from a straight line diagnose misfit.
This linearity holds for all four distributions. The (0,0) and (1,1)
endpoints map to infinity and are dropped on the transformed scale.

## See also

[`roc_sdt()`](https://popov-lab.github.io/bmm/dev/reference/roc_sdt.md),
[`roc_observed()`](https://popov-lab.github.io/bmm/dev/reference/roc_observed.md),
[`auc_sdt()`](https://popov-lab.github.io/bmm/dev/reference/auc_sdt.md)
