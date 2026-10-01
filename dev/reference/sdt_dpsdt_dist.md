# Distribution functions for dual-process SDT (DPSDT)

Density and random generation for the dual-process signal detection
model (Yonelinas, 1994), with recall-to-reject of new items (`Rn`; as in
Yonelinas, 2024). Extends rating SDT with recollection probabilities
`Ro` (old items recollected as old) and `Rn` (new items recall-rejected)
that add mass to the most-confident rating category. These are the
simulation counterparts of the `dpsdt` version of
[`sdt_rating()`](https://popov-lab.github.io/bmm/dev/reference/sdt_rating.md);
here `Ro`/`Rn` are supplied directly as probabilities in `[0, 1]`.

## Usage

``` r
dsdt_dpsdt(
  counts,
  stimulus,
  d,
  thresholds,
  Ro,
  Rn = 0,
  sdratio = 1,
  dist = c("normal", "gumbel_min", "gumbel_max", "logistic"),
  log = FALSE
)

rsdt_dpsdt(
  n,
  n_trials,
  stimulus,
  d,
  thresholds,
  Ro,
  Rn = 0,
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

- Ro:

  Numeric vector in `[0, 1]`. Recollection probability for old (signal)
  items.

- Rn:

  Numeric vector in `[0, 1]`. Recall-to-reject probability for new
  (noise) items. Defaults to 0, the classic one-sided model, as the
  `dpsdt` version fixes it off unless `Rn` is in the formula.

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

`dsdt_dpsdt` returns the (log-)density (multinomial probability).
`rsdt_dpsdt` returns an integer matrix with one row per observation and
one rating-count column per category (`r1` ... `rK`).

## References

Yonelinas, A. P. (1994). Receiver-operating characteristics in
recognition memory: Evidence for a dual-process model. *Journal of
Experimental Psychology: Learning, Memory, and Cognition*, *20*(6),
1341–1354.
[doi:10.1037/0278-7393.20.6.1341](https://doi.org/10.1037/0278-7393.20.6.1341)

Yonelinas, A. P. (2024). The role of recollection and familiarity in
visual working memory: A mixture of threshold and signal detection
processes. *Psychological Review*, *131*(2), 321–348.
[doi:10.1037/rev0000432](https://doi.org/10.1037/rev0000432)

## Examples

``` r
# Density for a single observation (K=4) with recollection of old items;
# Rn defaults to 0 (no recall-to-reject), the one-sided model
dsdt_dpsdt(counts = c(2, 8, 20, 70), stimulus = 1,
           d = 1.5, thresholds = c(-0.5, 0.0, 0.5), Ro = 0.3)
#> [1] 3.023007e-05
# Generate DPSDT rating data (K=4) for 10 subjects and both stimulus types
dat <- expand.grid(id = 1:10, stimulus = c(0L, 1L))
dat <- cbind(dat, rsdt_dpsdt(nrow(dat), 100, dat$stimulus, d = 1.5,
                             thresholds = c(-0.5, 0, 0.5),
                             Ro = 0.3, Rn = 0.1))
head(dat)
#>   id stimulus r1 r2 r3 r4
#> 1  1        0 59 15 12 14
#> 2  2        0 61 12 12 15
#> 3  3        0 76  9  8  7
#> 4  4        0 66 14 10 10
#> 5  5        0 64 18 12  6
#> 6  6        0 61 24  7  8
```
