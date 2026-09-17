# Three-parameter mixture model by Bays et al (2009).

Three-parameter mixture model by Bays et al (2009).

## Usage

``` r
mixture3p(resp_error, nt_features, set_size, regex = FALSE, ...)
```

## Arguments

- resp_error:

  The name of the variable in the dataset containing the response error.
  The response error should code the response relative to the
  to-be-recalled target in radians. You can transform the response error
  in degrees to radians using the `deg2rad` function.

- nt_features:

  A character vector with the names of the non-target feature values.
  The non_target feature values should be in radians and centered
  relative to the target. Alternatively, if regex=TRUE, a regular
  expression can be used to match the non-target feature columns in the
  dataset.

- set_size:

  Name of the column containing the set size variable (if set_size
  varies) or a numeric value for the set_size, if the set_size is fixed.

- regex:

  Logical. If TRUE, the `nt_features` argument is interpreted as a
  regular expression to match the non-target feature columns in the
  dataset.

- ...:

  used internally for testing, ignore it

## Value

An object of class `bmmodel`

## Details

- **Domain:** Visual working memory

- **Task:** Continuous reproduction

- **Name:** Three-parameter mixture model by Bays et al (2009).

- **Citation:**

  - Bays, P. M., Catalao, R. F. G., & Husain, M. (2009). The precision
    of visual working memory is set by allocation of a shared resource.
    Journal of Vision, 9(10), 1-11

- **Requirements:**

  - The response vairable should be in radians and represent the angular
    error relative to the target

&nbsp;

- The non-target features should be in radians and be centered relative
  to the target

&nbsp;

- **Parameters:**

  - `mu1`: Location parameter of the von Mises distribution for memory
    responses (in radians). Fixed internally to 0 by default.

  - `kappa`: Concentration parameter of the von Mises distribution

  - `thetat`: Mixture weight for target responses

  - `thetant`: Mixture weight for non-target responses

- **Fixed parameters:**

  - `mu1` = 0

  - `mu2` = 0

  - `kappa2` = -100

- **Default parameter links:**

  - mu1 = tan_half; kappa = log; thetat = softmax; thetant = softmax

- **Default priors:**

  - `mu1`:

    - `main`: student_t(1, 0, 1)

  - `kappa`:

    - `main`: normal(2, 1)

    - `effects`: normal(0, 1)

  - `thetat`:

    - `main`: logistic(0, 1)

  - `thetant`:

    - `main`: logistic(0, 1)

## Examples

``` r
if (FALSE) { # isTRUE(Sys.getenv("BMM_EXAMPLES"))
# generate artificial data from the Bays et al (2009) 3-parameter mixture model
dat <- data.frame(
  y = rmixture3p(n = 2000, mu = c(0, 1, -1.5, 2)),
  nt1_loc = 1,
  nt2_loc = -1.5,
  nt3_loc = 2
)

# define formula
ff <- bmmformula(
  kappa ~ 1,
  thetat ~ 1,
  thetant ~ 1
)

# specify the 3-parameter model with explicit column names for non-target features
model1 <- mixture3p(resp_error = "y", nt_features = paste0("nt", 1:3, "_loc"), set_size = 4)

# fit the model
fit <- bmm(
  formula = ff,
  data = dat,
  model = model1,
  cores = 4,
  iter = 500,
  backend = "cmdstanr"
)

# alternatively specify the 3-parameter model with a regular expression to match non-target features
# this is equivalent to the previous call, but more concise
model2 <- mixture3p(resp_error = "y", nt_features = "nt.*_loc", set_size = 4, regex = TRUE)

# fit the model
fit <- bmm(
  formula = ff,
  data = dat,
  model = model2,
  cores = 4,
  iter = 500,
  backend = "cmdstanr"
)
}
```
