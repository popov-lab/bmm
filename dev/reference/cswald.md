# Censored-Shifted Wald Model

Censored-Shifted Wald Model

## Usage

``` r
cswald(rt, response, links = NULL, version = c("simple", "crisk"), ...)
```

## Arguments

- rt:

  The name of the variable in the dataset containing the response times.
  Response times should be coded in seconds (not milliseconds).

- response:

  The name of the variable in the dataset containing the
  response/decision. Responses should be coded as 0 (lower boundary) or
  1 (upper boundary). Alternatively, character values "lower" and
  "upper" or logical values (FALSE/TRUE) are accepted and will be
  converted automatically.

- links:

  A named list of link functions for the model parameters. Available
  parameters depend on the version: "simple" has `drift`, `bound`,
  `ndt`, and `s`; "crisk" additionally has `zr`. Default links are "log"
  for most parameters and "logit" for `zr`. For positive parameters,
  "softplus" is available as an alternative to "log" that grows linearly
  for large values and avoids the numerical blow-up of
  [`exp()`](https://rdrr.io/r/base/Log.html).

- version:

  A character string specifying which version of the cswald model to
  use. Options are:

  - `"simple"` (default): The standard censored shifted Wald model,
    which treats error responses as censored correct responses. Best
    suited for tasks with few errors (\<20%). **Note:** The `bound`
    parameter in the simple version represents the distance from the
    starting point to the correct boundary, which is half the total
    boundary separation in the diffusion model (assuming an unbiased
    starting point). To convert to the full boundary separation (as in
    DDM or crisk), multiply by 2.

  - `"crisk"`: The competing risks version, which models both response
    types as arising from racing accumulators toward opposite
    boundaries. Better suited for tasks with substantial error rates.
    The `bound` parameter represents the total boundary separation,
    consistent with the diffusion model parameterization.

  For more details, see Miller et al. (2017).

- ...:

  Additional arguments passed internally (for testing purposes).

## Value

An object of class `bmmodel`

## Details

- **Domain:** Decision Making / Response times

- **Task:** Choice Reaction Time tasks (with few errors)

- **Name:** Censored-Shifted Wald Model

- **Citation:**

  - Miller, R., Scherbaum, S., Heck, D. W., Goschke, T., & Enge, S.
    (2017). On the Relation Between the (Censored) Shifted Wald and the
    Wiener Distribution as Measurement Models for Choice Response Times.
    Applied Psychological Measurement, 42(2), 116-135.
    https://doi.org/10.1177/0146621617710465

- **Version:** simple

- **Requirements:**

  - Reaction times should be passed in seconds

&nbsp;

- The response variable should be passed numerically: 0 = lower
  response, 1 = upper response

&nbsp;

- **Parameters:**

  - `drift`: drift rate

  - `bound`: boundary (distance from starting point to correct boundary)

  - `ndt`: non-decision time

  - `s`: diffusion constant

- **Fixed parameters:**

  - `mu` = 0

  - `s` = 0

- **Default parameter links:**

  - drift = log; bound = log; ndt = log; s = log

- **Default priors:**

  - `drift`:

    - `main`: normal(0,1)

    - `effects`: normal(0,0.3)

  - `bound`:

    - `main`: normal(0,0.3)

    - `effects`: normal(0,0.3)

  - `ndt`:

    - `main`: normal(-2,0.3)

    - `effects`: normal(0,0.3)

  - `s`:

    - `main`: normal(0,0.3)

    - `effects`: normal(0,0.2)

## See also

[`dcswald()`](https://venpopov.com/bmm/dev/reference/cswald_dist.md) and
[`rcswald()`](https://venpopov.com/bmm/dev/reference/cswald_dist.md) for
the density and random generation functions.

## Examples

``` r
if (FALSE) { # isTRUE(Sys.getenv("BMM_EXAMPLES"))
# generate simulated data from the diffusion model
dat <- rcswald(n = 500, drift = 2, bound = 1.5, ndt = 0.3, zr = 0.5, s = 1)

# specify the model
model <- cswald(rt = "rt", response = "response", version = "simple")

# specify the formula
formula <- bmf(
  drift ~ 1,
  bound ~ 1,
  ndt ~ 1
)

# fit the model
fit <- bmm(
  formula = formula,
  data = dat,
  model = model,
  cores = 4,
  backend = "cmdstanr"
)
}
```
