# Utility functions for Signal Detection Theory

Compute sensitivity and criterion from hit and false alarm rates for
different SDT distribution families. A single (hit, false alarm) pair
cannot identify the signal-to-noise SD ratio, so both quantities are the
equal-variance values; see
[`sdt_yn()`](https://venpopov.com/bmm/dev/reference/sdt_yn.md) for the
unequal-variance model.

## Usage

``` r
sdt_d(
  hit_rate,
  fa_rate,
  dist = c("normal", "gumbel_min", "gumbel_max", "logistic")
)

sdt_criterion(
  hit_rate,
  fa_rate,
  dist = c("normal", "gumbel_min", "gumbel_max", "logistic")
)
```

## Arguments

- hit_rate:

  Numeric. Proportion of hits (P("old" \| signal)).

- fa_rate:

  Numeric. Proportion of false alarms (P("old" \| noise)).

- dist:

  Character. The distribution assumed for the latent evidence, given
  here by its cumulative distribution function:

  - "normal" (default): Gaussian, \\\Phi(x)\\

  - "gumbel_min": smallest extreme value, \\1 - \exp(-\exp(x))\\ (the
    complementary log-log distribution)

  - "gumbel_max": largest extreme value, \\\exp(-\exp(-x))\\ (the
    log-log distribution, as in
    [`evd::pgumbel`](https://rdrr.io/pkg/evd/man/gumbel.html))

  - "logistic": \\1 / (1 + \exp(-x))\\

## Value

`sdt_d` returns the distance between the signal and noise distributions
on the latent evidence axis, obtained by inverting the decision rule
"respond old when the evidence exceeds the criterion": \\Q(1 - FA) -
Q(1 - H)\\, where \\Q\\ is the quantile function of `dist`. For
`dist = "normal"` this reduces to the familiar \\d' = \Phi^{-1}(H) -
\Phi^{-1}(FA)\\. Because one operating point implies equal variance,
this matches the `d` parameter of
[`sdt_yn()`](https://venpopov.com/bmm/dev/reference/sdt_yn.md) whenever
`sdratio` is at its default.

`sdt_criterion` returns the criterion (response bias) on the centred,
noise-standardized evidence axis used by
[`sdt_yn()`](https://venpopov.com/bmm/dev/reference/sdt_yn.md), where
the noise and signal distributions sit at -d'/2 and +d'/2: \\(Q(1 -
FA) + Q(1 - H)) / 2\\. For `dist = "normal"` this reduces to the
familiar \\-(\Phi^{-1}(H) + \Phi^{-1}(FA)) / 2\\.

## References

Green, D. M., & Swets, J. A. (1966). *Signal detection theory and
psychophysics*. Wiley.

## See also

[`sdt_yn()`](https://venpopov.com/bmm/dev/reference/sdt_yn.md), whose
`d` and `criterion` parameters these two functions compute in closed
form from observed rates: `sdt_d()` returns the same quantity as the `d`
parameter and `sdt_criterion()` the same quantity as `criterion`, on the
same axis, whenever `sdratio` is at its default.

## Examples

``` r
# Compute d from hit and false alarm rates (Gaussian SDT)
sdt_d(hit_rate = 0.8, fa_rate = 0.2, dist = "normal")
#> [1] 1.683242

# The extreme-value analogue
sdt_d(hit_rate = 0.75, fa_rate = 0.25, dist = "gumbel_min")
#> [1] 1.572534
# Compute criterion from hit and false alarm rates
sdt_criterion(hit_rate = 0.8, fa_rate = 0.2, dist = "normal")
#> [1] 0
```
