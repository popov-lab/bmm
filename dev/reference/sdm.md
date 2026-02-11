# Signal Discrimination Model (SDM) by Oberauer (2023)

Signal Discrimination Model (SDM) by Oberauer (2023)

## Usage

``` r
sdm(resp_error, version = "simple", ...)

sdmSimple(resp_error, version = "simple", ...)
```

## Arguments

- resp_error:

  The name of the variable in the dataset containing the response error.
  The response error should code the response relative to the
  to-be-recalled target in radians. You can transform the response error
  in degrees to radians using the `deg2rad` function.

- version:

  Character. The version of the model to use. Currently only "simple" is
  supported.

- ...:

  used internally for testing, ignore it

## Value

An object of class `bmmodel`

## Details

see [the online
article](https://venpopov.com/bmm/articles/bmm_sdm_simple.html) for a
detailed description of the model and how to use it. \* **Domain:**
Visual working memory

- **Task:** Continuous reproduction

- **Name:** Signal Discrimination Model (SDM) by Oberauer (2023)

- **Citation:**

  - Oberauer, K. (2023). Measurement models for visual working memory -
    A factorial model comparison. Psychological Review, 130(3), 841-852

- **Version:** simple

- **Requirements:**

  - The response variable should be in radians and represent the angular
    error relative to the target

- **Parameters:**

  - `mu`: Location parameter of the SDM distribution (in radians; by
    default fixed internally to 0)

  - `c`: Memory strength parameter of the SDM distribution

  - `kappa`: Precision parameter of the SDM distribution

- **Fixed parameters:**

  - `mu` = 0

- **Default parameter links:**

  - mu = tan_half; c = log; kappa = log

- **Default priors:**

  - `mu`:

    - `main`: student_t(1, 0, 1)

  - `kappa`:

    - `main`: student_t(5, 1.75, 0.75)

    - `effects`: normal(0, 1)

  - `c`:

    - `main`: student_t(5, 2, 0.75)

    - `effects`: normal(0, 1)

## Examples

``` r
if (FALSE) { # isTRUE(Sys.getenv("BMM_EXAMPLES"))
# simulate data from the model
dat <- data.frame(y = rsdm(n = 1000, c = 4, kappa = 3))

# specify formula
ff <- bmf(
  c ~ 1,
  kappa ~ 1
)

# specify the model
fit <- bmm(
  formula = ff,
  data = dat,
  model = sdm(resp_error = "y"),
  cores = 4,
  backend = "cmdstanr"
)
}
```
