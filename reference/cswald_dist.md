# Distribution functions for the censored shifted Wald model (`cswald`)

These functions provide the density, distribution, quantile, and random
generation functions for the censored shifted Wald model: `cswald`. The
random generation (`rcswald`) and distribution functions (`pcswald`,
`qcswald`) use
[`rtdists::rdiffusion()`](https://rdrr.io/pkg/rtdists/man/Diffusion.html),
[`rtdists::pdiffusion()`](https://rdrr.io/pkg/rtdists/man/Diffusion.html),
and
[`rtdists::qdiffusion()`](https://rdrr.io/pkg/rtdists/man/Diffusion.html)
internally for the `"crisk"` version, which is theoretically consistent
since the censored shifted Wald model is an approximation to the Wiener
diffusion model for tasks with high accuracy (few errors).

## Usage

``` r
dcswald(
  rt,
  response,
  drift,
  bound,
  ndt,
  zr = 0.5,
  s = 1,
  version = c("simple", "crisk"),
  log = TRUE
)

rcswald(n, drift, bound, ndt, zr = 0.5, s = 1)

pcswald(
  q,
  response,
  drift,
  bound,
  ndt,
  zr = 0.5,
  s = 1,
  version = "simple",
  lower.tail = TRUE,
  log.p = FALSE
)

qcswald(
  p,
  response,
  drift,
  bound,
  ndt,
  zr = 0.5,
  s = 1,
  version = "simple",
  lower.tail = TRUE,
  log.p = FALSE
)
```

## Arguments

- rt:

  A vector of response times in seconds for which the likelihood should
  be evaluated

- response:

  A vector of responses coded numerically: 0 = lower response, 1 = upper
  response

- drift:

  The drift rate

- bound:

  The boundary separation

- ndt:

  The non-decision time

- zr:

  The relative starting point (proportion of boundary separation).
  Default is `0.5` (unbiased). Values must be between 0 and 1.

- s:

  The diffusion constant - the standard deviation of the noise in the
  evidence accumulation process. Default is `s = 1`

- version:

  A character string specifying the version of the `cswald` for which
  the likelihood should be returned. Available versions are "simple" and
  "crisk", the default is "simple."

- log:

  A single logical value indicating if log-likelihoods should be
  returned, the default is `TRUE`

- n:

  The number of random samples that should be generated

- q:

  A vector of quantiles (response times) at which to evaluate the CDF

- lower.tail:

  Logical; if `TRUE` (default), probabilities are P(RT \<= q), otherwise
  P(RT \> q)

- log.p:

  Logical; if `TRUE`, probabilities are returned on the log scale.
  Default is `FALSE`

- p:

  A vector of probabilities for which to compute quantiles

## Value

- `dcswald()` returns a numeric vector of (log-)likelihoods.

- `rcswald()` returns a data.frame with columns `rt` (response times)
  and `response` (1 = upper, 0 = lower).

- `pcswald()` returns a numeric vector of (log-)probabilities.

- `qcswald()` returns a numeric vector of quantiles (response times).

## Details

**Cumulative Distribution Function (`pcswald`)**

For the `"simple"` version, the CDF is only defined for `response = 1`
(correct responses), as errors are treated as censored observations. The
CDF returns the probability that a correct response occurs by time `q`.
For `response = 0`, `NA` is returned with a warning.

For the `"crisk"` version (competing risks), the CDF computes the
defective cumulative distribution P(RT \<= q, response = r), which is
the probability of responding with the specified response by time `q`.
This uses

[`rtdists::pdiffusion()`](https://rdrr.io/pkg/rtdists/man/Diffusion.html)
internally for accurate computation.

**Quantile Function (`qcswald`)**

The quantile function returns the response time `q` such that P(RT \<=
q) = p. Similar to the CDF, for the `"simple"` version this is only
defined for `response = 1`. For the `"crisk"` version, this uses
[`rtdists::qdiffusion()`](https://rdrr.io/pkg/rtdists/man/Diffusion.html)
internally.

## See also

[`rtdists::rdiffusion()`](https://rdrr.io/pkg/rtdists/man/Diffusion.html),
[`rtdists::pdiffusion()`](https://rdrr.io/pkg/rtdists/man/Diffusion.html),
[`rtdists::qdiffusion()`](https://rdrr.io/pkg/rtdists/man/Diffusion.html)
for the underlying functions

## Examples

``` r
dat <- rcswald(n = 1000, drift = 2, bound = 1, ndt = 0.3)
head(dat)
#>          rt response
#> 1 0.3586012        1
#> 2 0.3997922        1
#> 3 0.7422927        1
#> 4 0.8045497        1
#> 5 0.5244073        1
#> 6 0.5207823        1
hist(dat$rt)
```
