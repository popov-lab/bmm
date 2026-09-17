# Interference measurement model by Oberauer and Lin (2017).

Three versions of the Interference measurement model by Oberauer and Lin
(2017). - the full, bsc, and abc. `IMMfull()`, `IMMbsc()`, and
`IMMabc()` are deprecated and will be removed in the future. Please use
`imm(version = 'full')`, `imm(version = 'bsc')`, or
`imm(version = 'abc')` instead.

## Usage

``` r
imm(
  resp_error,
  nt_features,
  nt_distances,
  set_size,
  regex = FALSE,
  version = "full",
  ...
)

IMMfull(resp_error, nt_features, nt_distances, set_size, regex = FALSE, ...)

IMMbsc(resp_error, nt_features, nt_distances, set_size, regex = FALSE, ...)

IMMabc(resp_error, nt_features, set_size, regex = FALSE, ...)
```

## Arguments

- resp_error:

  The name of the variable in the provided dataset containing the
  response error. The response Error should code the response relative
  to the to-be-recalled target in radians. You can transform the
  response error in degrees to radian using the `deg2rad` function.

- nt_features:

  A character vector with the names of the non-target variables. The
  non_target variables should be in radians and be centered relative to
  the target. Alternatively, if regex=TRUE, a regular expression can be
  used to match the non-target feature columns in the dataset.

- nt_distances:

  A vector of names of the columns containing the distances of
  non-target items to the target item. Alternatively, if regex=TRUE, a
  regular expression can be used to match the non-target distances
  columns in the dataset. Only necessary for the `bsc` and `full`
  versions.

- set_size:

  Name of the column containing the set size variable (if set_size
  varies) or a numeric value for the set_size, if the set_size is fixed.

- regex:

  Logical. If TRUE, the `nt_features` and `nt_distances` arguments are
  interpreted as a regular expression to match the non-target feature
  columns in the dataset.

- version:

  Character. The version of the IMM model to use. Can be one of `full`,
  `bsc`, or `abc`. The default is `full`.

- ...:

  used internally for testing, ignore it

## Value

An object of class `bmmodel`

## Details

- **Domain:** Visual working memory

- **Task:** Continuous reproduction

- **Name:** Interference measurement model by Oberauer and Lin (2017).

- **Citation:**

  - Oberauer, K., & Lin, H.Y. (2017). An interference model of visual
    working memory. Psychological Review, 124(1), 21-59

### Version: `full`

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

  - `a`: General activation of memory items

  - `c`: Context activation

  - `s`: Spatial similarity gradient

- **Fixed parameters:**

  - `mu1` = 0

  - `mu2` = 0

  - `kappa2` = -100

- **Default parameter links:**

  - mu1 = tan_half; kappa = log; a = log; c = log; s = log

- **Default priors:**

  - `mu1`:

    - `main`: student_t(1, 0, 1)

  - `kappa`:

    - `main`: normal(2, 1)

    - `effects`: normal(0, 1)

  - `a`:

    - `main`: normal(0, 1)

    - `effects`: normal(0, 1)

  - `c`:

    - `main`: normal(0, 1)

    - `effects`: normal(0, 1)

  - `s`:

    - `main`: normal(0, 1)

    - `effects`: normal(0, 1)

### Version: `bsc`

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

  - `c`: Context activation

  - `s`: Spatial similarity gradient

- **Fixed parameters:**

  - `mu1` = 0

  - `mu2` = 0

  - `kappa2` = -100

- **Default parameter links:**

  - mu1 = tan_half; kappa = log; c = log; s = log

- **Default priors:**

  - `mu1`:

    - `main`: student_t(1, 0, 1)

  - `kappa`:

    - `main`: normal(2, 1)

    - `effects`: normal(0, 1)

  - `c`:

    - `main`: normal(0, 1)

    - `effects`: normal(0, 1)

  - `s`:

    - `main`: normal(0, 1)

    - `effects`: normal(0, 1)

### Version: `abc`

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

  - `a`: General activation of memory items

  - `c`: Context activation

- **Fixed parameters:**

  - `mu1` = 0

  - `mu2` = 0

  - `kappa2` = -100

- **Default parameter links:**

  - mu1 = tan_half; kappa = log; a = log; c = log

- **Default priors:**

  - `mu1`:

    - `main`: student_t(1, 0, 1)

  - `kappa`:

    - `main`: normal(2, 1)

    - `effects`: normal(0, 1)

  - `a`:

    - `main`: normal(0, 1)

    - `effects`: normal(0, 1)

  - `c`:

    - `main`: normal(0, 1)

    - `effects`: normal(0, 1)

Additionally, all imm models have an internal parameter that is fixed to
0 to allow the model to be identifiable. This parameter is not estimated
and is not included in the model formula. The parameter is:

- b = "Background activation (internally fixed to 0)"

## Examples

``` r
if (FALSE) { # isTRUE(Sys.getenv("BMM_EXAMPLES"))
# load data
data <- oberauer_lin_2017

# define formula
ff <- bmmformula(
  kappa ~ 0 + set_size,
  c ~ 0 + set_size,
  a ~ 0 + set_size,
  s ~ 0 + set_size
)

# specify the full IMM model with explicit column names for non-target features and distances
# by default this fits the full version of the model
model1 <- imm(
  resp_error = "dev_rad",
  nt_features = paste0("col_nt", 1:7),
  nt_distances = paste0("dist_nt", 1:7),
  set_size = "set_size"
)

# fit the model
fit <- bmm(
  formula = ff,
  data = data,
  model = model1,
  cores = 4,
  backend = "cmdstanr"
)

# alternatively specify the IMM model with a regular expression to match non-target features
# this is equivalent to the previous call, but more concise
model2 <- imm(
  resp_error = "dev_rad",
  nt_features = "col_nt",
  nt_distances = "dist_nt",
  set_size = "set_size",
  regex = TRUE
)

# fit the model
fit <- bmm(
  formula = ff,
  data = data,
  model = model2,
  cores = 4,
  backend = "cmdstanr"
)

# you can also specify the `bsc` or `abc` versions of the model to fit a reduced version
model3 <- imm(
  resp_error = "dev_rad",
  nt_features = "col_nt",
  set_size = "set_size",
  regex = TRUE,
  version = "abc"
)
fit <- bmm(
  formula = ff,
  data = data,
  model = model3,
  cores = 4,
  backend = "cmdstanr"
)
}
```
