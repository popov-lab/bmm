# Two-parameter mixture model by Zhang and Luck (2008).

Two-parameter mixture model by Zhang and Luck (2008).

## Usage

``` r
mixture2p(resp_error, ...)
```

## Arguments

- resp_error:

  The name of the variable in the provided dataset containing the
  response error. The response Error should code the response relative
  to the to-be-recalled target in radians. You can transform the
  response error in degrees to radian using the `deg2rad` function.

- ...:

  used internally for testing, ignore it

## Value

An object of class `bmmodel`

## Details

- **Domain:** Visual working memory

- **Task:** Continuous reproduction

- **Name:** Two-parameter mixture model by Zhang and Luck (2008).

- **Citation:**

  - Zhang, W., & Luck, S. J. (2008). Discrete fixed-resolution
    representations in visual working memory. Nature, 453(7192), 233-235

- **Requirements:**

  - The response vairable should be in radians and represent the angular
    error relative to the target

- **Parameters:**

  - `mu1`: Location parameter of the von Mises distribution for memory
    responses (in radians). Fixed internally to 0 by default.

  - `kappa`: Concentration parameter of the von Mises distribution

  - `thetat`: Mixture weight for target responses

- **Fixed parameters:**

  - `mu1` = 0

  - `mu2` = 0

  - `kappa2` = -100

- **Default parameter links:**

  - mu1 = tan_half; kappa = log; thetat = identity

- **Default priors:**

  - `mu1`:

    - `main`: student_t(1, 0, 1)

  - `kappa`:

    - `main`: normal(2, 1)

    - `effects`: normal(0, 1)

  - `thetat`:

    - `main`: logistic(0, 1)

## Examples

``` r
if (FALSE) { # isTRUE(Sys.getenv("BMM_EXAMPLES"))
# generate artificial data
dat <- data.frame(y = rmixture2p(n = 2000))

# define formula
ff <- bmmformula(kappa ~ 1, thetat ~ 1)

model <- mixture2p(resp_error = "y")

# fit the model
fit <- bmm(
  formula = ff,
  data = dat,
  model = model,
  cores = 4,
  iter = 500,
  backend = "cmdstanr"
)
}
```
