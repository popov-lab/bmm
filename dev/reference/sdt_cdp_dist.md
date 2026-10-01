# Distribution functions for Continuous Dual-Process SDT (CDP)

Density and random generation for the continuous dual-process signal
detection theory model (Wixted & Mickes, 2010). Two correlated
continuous dimensions, Familiarity (F) and Recollection (R), generate
the aggregate strength S = F + R that drives old/new confidence;
Remember/Know judgments split "old" responses on R, and an optional
Know/Guess split uses F. The response is a vector of counts across the
response categories (multinomial likelihood), in the same format
[`sdt_cdp()`](https://popov-lab.github.io/bmm/dev/reference/sdt_cdp.md)
is fit to. See
[`sdt_cdp()`](https://popov-lab.github.io/bmm/dev/reference/sdt_cdp.md)
for the model and
[`sdt_rating()`](https://popov-lab.github.io/bmm/dev/reference/sdt_rating.md)
when no R/K split is available.

## Usage

``` r
dsdt_cdp(
  counts,
  stimulus,
  dfam,
  drec,
  thresholds,
  rcrit,
  kcrit = NULL,
  sigmar = 0,
  rho = 0,
  n_new = NULL,
  dist = "normal",
  log = FALSE
)

rsdt_cdp(
  n,
  n_trials,
  stimulus,
  dfam,
  drec,
  thresholds,
  rcrit,
  kcrit = NULL,
  sigmar = 0,
  rho = 0,
  n_new = NULL,
  dist = "normal"
)
```

## Arguments

- counts:

  Integer matrix with one row per observation and one column per
  response category in the canonical order `new`, `[guess]`, `know`,
  `remember` (each block ordered by confidence; see
  [`sdt_cdp()`](https://popov-lab.github.io/bmm/dev/reference/sdt_cdp.md)),
  or a vector for a single observation. `n_new + 2 * n_old` columns
  (R/K) or `n_new + 3 * n_old` (R/K/G).

- stimulus:

  Integer vector (0/1). Stimulus type: 0 = new/lure, 1 = old/target.

- dfam, drec:

  Numeric vectors. Familiarity and recollection sensitivities: the
  target means on the F and R axes, each in units of the corresponding
  lure SD (both lure SDs are 1). They are component means of a bivariate
  latent space, not the balanced \\d_a\\ that the other SDT models
  report as `d` – the old/new decision here is read off the aggregate
  strength S = F + R, so the model's discriminability is a derived
  quantity rather than either of these.

- thresholds:

  Numeric vector of `n_new + n_old - 1` ordered confidence thresholds on
  the aggregate F+R axis, or a matrix with one row per observation.

- rcrit:

  Numeric vector. Remember criterion on the recollection axis.

- kcrit:

  Know criterion on the familiarity axis. `NULL` (default) for the R/K
  model; finite values enable the R/K/G model.

- sigmar:

  Numeric vector. Log SD of the recollection target distribution (0 = SD
  1).

- rho:

  Numeric vector. F-R correlation on the unconstrained scale;
  `tanh(rho)` is the correlation. 0 (default) = independent processes
  (classic CDP).

- n_new:

  Integer number of "new" confidence levels. Defaults to half the number
  of confidence levels, rounded down (so 3 of 7); the number of "old"
  levels follows as `length(thresholds) + 1 - n_new`.

- dist:

  Noise distribution. Only `"normal"` is currently supported.

- log:

  Logical; if `TRUE` return the log-density (default `FALSE`).

- n:

  Integer. Number of observations to generate. `n_trials`, `stimulus`,
  `thresholds`, and the model parameters are recycled to this length.

- n_trials:

  Integer vector. Number of trials per observation.

## Value

`dsdt_cdp` returns the (log-)density (multinomial probability).
`rsdt_cdp` returns an integer matrix with one row per observation and
the canonical response count columns (`new1`, ..., `remember<K>`) that
[`sdt_cdp()`](https://popov-lab.github.io/bmm/dev/reference/sdt_cdp.md)
expects – [`cbind()`](https://rdrr.io/r/base/cbind.html) it to a design
data frame for a ready-to-fit data set.

## References

Wixted, J. T., & Mickes, L. (2010). A continuous dual-process model of
remember/know judgments. *Psychological Review*, *117*(4), 1025–1054.
[doi:10.1037/a0020874](https://doi.org/10.1037/a0020874)

## Examples

``` r
# CDP density (R/K, 1 new + 3 old levels: 3 thresholds, 7 count columns)
dsdt_cdp(
  counts = c(40, 5, 12, 30, 10, 18, 60), stimulus = 1,
  dfam = 0.8, drec = 1.0,
  thresholds = c(-0.5, 0.3, 1.0),
  rcrit = 0.7, n_new = 1
)
#> [1] 1.870555e-27
# Generate CDP count data (R/K, 3 new + 3 old levels) for 10 subjects
dat <- expand.grid(id = 1:10, stimulus = c(0L, 1L))
dat <- cbind(dat, rsdt_cdp(nrow(dat), 100, dat$stimulus,
                           dfam = 0.8, drec = 1.0,
                           thresholds = c(-1.1, -0.5, 0, 0.6, 1.3),
                           rcrit = 0.5, n_new = 3))
head(dat)
#>   id stimulus new1 new2 new3 know4 know5 know6 remember4 remember5 remember6
#> 1  1        0   19   17   18    14     6     2         7        10         7
#> 2  2        0   19   18   10    10    13     1         7         8        14
#> 3  3        0   21   15    9    16     3     5         7        11        13
#> 4  4        0   22   10   12    13     6     6         8         9        14
#> 5  5        0   15   13   13     4    13     5         6        10        21
#> 6  6        0   17   13   12    11     9     5         7        12        14
```
