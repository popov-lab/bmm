# Rating SDT category log-probability (multinomial logit)

R companion to the Stan `sdt_rating_logmu` function. It returns
`log(p_k)` for rating category `cat`, which the
[`sdt_rating()`](https://popov-lab.github.io/bmm/dev/reference/sdt_rating.md)
multinomial formula uses as the category logit (so `softmax` recovers
the SDT category probabilities). `brms` evaluates the non-linear formula
in R for `posterior_predict()` and `posterior_epred()`, so this function
must be on the search path; it is exported for that reason and is not
called directly.

## Usage

``` r
sdt_rating_logmu(
  cat,
  K,
  dist,
  thresh,
  d,
  criterion,
  spacing,
  sdratio,
  stimulus,
  ...
)

sdt_dpsdt_logmu(
  cat,
  K,
  dist,
  thresh,
  d,
  criterion,
  spacing,
  sdratio,
  Ro,
  Rn,
  stimulus,
  ...
)

sdt_metad_logmu(
  cat,
  K,
  dist,
  thresh,
  d,
  criterion,
  spacing,
  sdratio,
  logmratio,
  stimulus,
  ...
)
```

## Arguments

- cat:

  Integer rating category index.

- K:

  Integer number of rating categories.

- dist:

  Integer noise-distribution id (see the `.sdt_dists` registry).

- thresh:

  Integer threshold-parameterization id.

- d, criterion, spacing, sdratio:

  Model parameters (draws-by-observation matrices supplied by brms). `d`
  is d', or d_a when sdratio is not 0; `spacing` is `0` for threshold
  types without it.

- stimulus:

  Stimulus covariate (0 = noise, 1 = signal).

- ...:

  Threshold `delta` parameters, when the threshold type uses them.

- Ro, Rn:

  Recollection parameters of the `dpsdt` version on the logit
  (linear-predictor) scale; `inv_logit(Ro)`/`inv_logit(Rn)` are the
  recollection and recall-to-reject probabilities.

- logmratio:

  Log M-ratio for the `metad` version; meta-d' is recovered as
  `exp(logmratio) * d`, on the same scale as `d`.

## Value

`log(p_cat)`, matching the shape of `d`.
