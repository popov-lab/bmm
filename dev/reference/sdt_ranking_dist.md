# Distribution functions for Ranking SDT

Density and random generation for ranking signal detection theory
(Meyer-Grant et al., 2026). Models rank ordering of m items by perceived
strength. Only `d` is estimated (no criterion or stimulus column).
Supports Gumbel-min (closed form) and Gaussian UV-SDT (numerical
integration).

## Usage

``` r
dsdt_ranking(
  counts,
  m,
  d,
  sdratio = 1,
  dist = c("gumbel_min", "normal"),
  log = FALSE
)

rsdt_ranking(n, n_trials, m, d, sdratio = 1, dist = c("gumbel_min", "normal"))
```

## Arguments

- counts:

  Integer matrix with one row per observation and one rank-count column
  per rank position (1 = most likely target), or a vector for a single
  observation. Columns beyond a row's set size `m` must be 0.

- m:

  Integer vector. Number of ranked items per observation. Must be at
  least 2 and no larger than the number of count columns.

- d:

  Numeric vector. Sensitivity: the distance between the target and lure
  distributions. It is \\d'\\ when `sdratio` is 1 and, for
  `dist = "gumbel_min"`, the \\g'\\ of Meyer-Grant et al. (2026). With
  another `sdratio` it is the balanced index \\d_a\\ that
  [`sdt_yn()`](https://popov-lab.github.io/bmm/dev/reference/sdt_yn.md)
  reports (in root-mean-square SD units).

- sdratio:

  Numeric vector. Ratio of signal to noise standard deviations (default
  1, i.e., equal variance). Must be positive. Only used when
  `dist = "normal"`.

- dist:

  Character. The distribution assumed for the latent evidence:
  "gumbel_min" (default), the smallest extreme value distribution with
  cumulative distribution function \\1 - \exp(-\exp(x))\\, evaluated in
  closed form; or "normal", Gaussian UV-SDT by numerical integration.

- log:

  Logical. If `TRUE`, returns log-density (default `FALSE`).

- n:

  Integer. Number of observations to generate. `n_trials`, `m`, `d`, and
  `sdratio` are recycled to this length.

- n_trials:

  Integer vector. Number of ranking trials per observation.

## Value

`dsdt_ranking` returns the (log-)density (multinomial probability).
`rsdt_ranking` returns an integer matrix with one row per observation
and one rank-count column per rank position (`rank1` ... `rank max(m)`);
rows with a smaller set size have structural zeros in the surplus
columns, matching the wide format
[`sdt_ranking()`](https://popov-lab.github.io/bmm/dev/reference/sdt_ranking.md)
expects.

## Parameter scales

These functions take `sdratio` as the ratio itself, while
[`sdt_ranking()`](https://popov-lab.github.io/bmm/dev/reference/sdt_ranking.md)
estimates its logarithm (0 = equal variance): a fitted `sdratio` of
0.375 is a ratio of `exp(0.375) = 1.455`, and passing `0.375` here
instead asks for a signal distribution 2.7 times narrower than the
noise. That is a legal value and raises no error, so exponentiate first.
`d` carries across unchanged.

## References

Meyer-Grant, C. G., Kellen, D., Harding, S. M., & Singmann, H. (2026).
Extreme-value signal detection theory for recognition memory: The
parametric road not taken. *Psychological Review*. Advance online
publication.
[doi:10.1037/rev0000615](https://doi.org/10.1037/rev0000615)

## Examples

``` r
# Gumbel-min ranking density
dsdt_ranking(counts = c(40, 30, 20, 10), m = 4, d = 1.0)
#> [1] 4.92936e-06
# Generate ranking data (m=4, Gumbel-min) for 10 subjects
dat <- data.frame(id = 1:10, set_size = 4L)
dat <- cbind(dat, rsdt_ranking(10, 100, m = 4, d = 1.0))
head(dat)
#>   id set_size rank1 rank2 rank3 rank4
#> 1  1        4    47    34     8    11
#> 2  2        4    53    23    16     8
#> 3  3        4    59    19    11    11
#> 4  4        4    52    16    17    15
#> 5  5        4    64    15    12     9
#> 6  6        4    57    14    19    10
```
