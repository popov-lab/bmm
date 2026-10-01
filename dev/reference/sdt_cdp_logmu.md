# CDP category log-probability (multinomial logit)

R companion to the Stan `sdt_cdp_logmu` function. It returns
`log(p_cat)` for response category `cat` (unnormalized, exactly like the
Stan function – the category probabilities sum to 1 analytically and
`softmax` absorbs the shared constant), which the
[`sdt_cdp()`](https://popov-lab.github.io/bmm/dev/reference/sdt_cdp.md)
multinomial formula uses as the category logit. `brms` evaluates the
non-linear formula in R for `posterior_predict()` and
`posterior_epred()`, so this function must be on the search path; it is
exported for that reason and is not called directly. Vectorized over the
draws-by-observations matrices brms supplies.

## Usage

``` r
sdt_cdp_logmu(
  cat,
  n_new,
  n_old,
  thresh,
  has_guess,
  dfam,
  drec,
  criterion,
  spacing,
  rcrit,
  sigmar,
  rho,
  kcrit,
  stimulus,
  ...
)
```

## Arguments

- cat:

  Integer response-category index (canonical order: new, then
  guess/know/remember blocks).

- n_new, n_old:

  Integer numbers of "new" and "old" confidence levels.

- thresh:

  Integer threshold-parameterization id (1 parsimonious, 2 equidistant,
  3 log_distance).

- has_guess:

  Integer flag (1 if the Know/Guess split is active).

- dfam, drec, criterion, spacing, rcrit, sigmar, rho, kcrit:

  Model parameters (draws-by-observation matrices supplied by brms).
  `kcrit` is ignored when `has_guess` is 0. `spacing` is the fixed
  literal 0 for log_distance, which uses the per-distance deltas in
  `...` instead.

- stimulus:

  Stimulus covariate (0 = new/lure, 1 = old/target).

- ...:

  For `log_distance`, the per-distance threshold parameters (`delta1`,
  ...) supplied by brms in contiguous order; otherwise unused.

## Value

`log(p_cat)`, matching the shape of `dfam`.
