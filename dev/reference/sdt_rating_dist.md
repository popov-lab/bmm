# Distribution functions for Confidence Rating SDT

Density and random generation for confidence rating signal detection
theory models. The response is a vector of counts across K ordered
rating categories (multinomial likelihood).

## Usage

``` r
dsdt_rating(
  counts,
  stimulus,
  d,
  thresholds,
  sdratio = 1,
  dist = c("normal", "gumbel_min", "gumbel_max", "logistic"),
  log = FALSE
)

rsdt_rating(
  n,
  n_trials,
  stimulus,
  d,
  thresholds,
  sdratio = 1,
  dist = c("normal", "gumbel_min", "gumbel_max", "logistic")
)
```

## Arguments

- counts:

  Integer matrix with one row per observation and one column per rating
  category, ordered from "definitely noise" (1) to "definitely signal"
  (K), or a vector for a single observation.

- stimulus:

  Integer vector (0/1). Stimulus type: 0 = noise, 1 = signal.

- d:

  Numeric vector. Sensitivity: \\d'\\ when `sdratio` is 1, and otherwise
  the balanced index \\d_a\\ (see
  [`sdt_yn()`](https://popov-lab.github.io/bmm/dev/reference/sdt_yn.md)).
  The separation between the distributions in noise units is
  `d * sqrt((1 + sdratio^2) / 2)`.

- thresholds:

  Numeric vector of length K-1 with the ordered decision thresholds, or
  an n-by-(K-1) matrix with one row per observation. The thresholds are
  positions on the noise-standardized axis, measured from the mean of
  the noise distribution (the signal mean sits at
  `d * sqrt((1 + sdratio^2) / 2)`), and are not rescaled by `sdratio`.

- sdratio:

  Numeric vector. Ratio of signal to noise standard deviations (default
  1, i.e., equal variance). Must be positive. Note that this is the
  natural scale: the `sdratio` parameter of
  [`sdt_rating()`](https://popov-lab.github.io/bmm/dev/reference/sdt_rating.md)
  is sampled on the log scale, so it corresponds to `log(sdratio)` here.

- dist:

  Character. Noise distribution: "normal" (default), "logistic",
  "gumbel_min", or "gumbel_max".

- log:

  Logical. If `TRUE`, returns log-density (default `FALSE`).

- n:

  Integer. Number of observations to generate. `n_trials`, `stimulus`,
  `thresholds`, and the model parameters are recycled to this length.

- n_trials:

  Integer vector. Number of trials per observation.

## Value

`dsdt_rating` returns the (log-)density (multinomial probability).
`rsdt_rating` returns an integer matrix with one row per observation and
one rating-count column per category (`r1` ... `rK`).

## References

Green, D. M., & Swets, J. A. (1966). *Signal detection theory and
psychophysics*. Wiley.

Selker, R., van den Bergh, D., Criss, A. H., & Wagenmakers, E.-J.
(2019). Parsimonious estimation of signal detection models from
confidence ratings. *Behavior Research Methods*, *51*(5), 1953–1967.
[doi:10.3758/s13428-019-01231-3](https://doi.org/10.3758/s13428-019-01231-3)

## Examples

``` r
# Density for a single observation (K=4)
dsdt_rating(counts = c(5, 15, 25, 55), stimulus = 1,
            d = 1.5, thresholds = c(0.25, 0.75, 1.25))
#> [1] 4.384333e-05
# Generate rating data (K=4) for 10 subjects and both stimulus types
dat <- expand.grid(id = 1:10, stimulus = c(0L, 1L))
dat <- cbind(dat, rsdt_rating(nrow(dat), 100, dat$stimulus,
                              d = 1.5, thresholds = c(0.25, 0.75, 1.25)))
head(dat)
#>   id stimulus r1 r2 r3 r4
#> 1  1        0 57 15 15 13
#> 2  2        0 67 18  7  8
#> 3  3        0 67 17  8  8
#> 4  4        0 69 13  9  9
#> 5  5        0 59 20 11 10
#> 6  6        0 58 14 20  8
```
