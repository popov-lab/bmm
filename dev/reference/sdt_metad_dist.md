# Distribution functions for meta-d' SDT

Density and random generation for the meta-d' model (Maniscalco & Lau,
2012). Confidence thresholds are placed using the metacognitive
sensitivity `metad`, then rescaled so the total "old"/"new" response
rates match what type-1 `d` predicts. These are the simulation
counterparts of the `metad` version of
[`sdt_rating()`](https://popov-lab.github.io/bmm/dev/reference/sdt_rating.md).
The old/new boundary is the middle threshold, so the number of rating
categories must be even (an odd number of `thresholds`).

## Usage

``` r
dsdt_metad(
  counts,
  stimulus,
  d,
  thresholds,
  metad,
  sdratio = 1,
  dist = c("normal", "gumbel_min", "gumbel_max", "logistic"),
  log = FALSE
)

rsdt_metad(
  n,
  n_trials,
  stimulus,
  d,
  thresholds,
  metad,
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
  on the noise-standardized axis and are not rescaled by `sdratio`.

- metad:

  Numeric vector. Metacognitive sensitivity (type-2 sensitivity), on the
  same scale as `d` (\\d'\\, or \\d_a\\ when `sdratio` is not 1), so the
  M-ratio `metad / d` is unaffected by `sdratio`. `metad = d`
  corresponds to ideal metacognition (recovers rating SDT).

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

`dsdt_metad` returns the (log-)density (multinomial probability).
`rsdt_metad` returns an integer matrix with one row per observation and
one rating-count column per category (`r1` ... `rK`).

## References

Maniscalco, B., & Lau, H. (2012). A signal detection theoretic approach
for estimating metacognitive sensitivity from confidence ratings.
*Consciousness and Cognition*, *21*(1), 422–430.
[doi:10.1016/j.concog.2011.09.021](https://doi.org/10.1016/j.concog.2011.09.021)

## Examples

``` r
# Density for a single observation (K=4) with imperfect metacognition
dsdt_metad(counts = c(5, 15, 25, 55), stimulus = 1,
           d = 1.5, thresholds = c(-0.5, 0.0, 0.5), metad = 1.0)
#> [1] 6.662953e-05
# Generate meta-d' rating data (K=4) for 10 subjects and both stimulus types
dat <- expand.grid(id = 1:10, stimulus = c(0L, 1L))
dat <- cbind(dat, rsdt_metad(nrow(dat), 100, dat$stimulus, d = 1.5,
                             thresholds = c(-0.5, 0, 0.5), metad = 1.0))
head(dat)
#>   id stimulus r1 r2 r3 r4
#> 1  1        0 53 31  8  8
#> 2  2        0 62 13 11 14
#> 3  3        0 57 22  8 13
#> 4  4        0 56 19 11 14
#> 5  5        0 50 23 16 11
#> 6  6        0 62 15 18  5
```
