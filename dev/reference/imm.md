# Interference Measurement Model

Three versions of the Interference Measurement Model - the full, bsc,
and abc. `IMMfull()`, `IMMbsc()`, and `IMMabc()` are deprecated and will
be removed in the future. Please use `imm(version = 'full')`,
`imm(version = 'bsc')`, or `imm(version = 'abc')` instead.

## Usage

``` r
imm(
  resp_error,
  nt_features,
  nt_distances,
  set_size,
  regex = FALSE,
  version = c("full", "bsc", "abc"),
  variable_precision = FALSE,
  vp_nodes = 81L,
  links = NULL,
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

- variable_precision:

  Logical; if `TRUE`, the precision of memory varies from trial to trial
  (van den Berg et al., 2012). This adds the parameter `tau`, and
  `kappa` then describes the average trial: it is the concentration
  whose Fisher information `J(kappa)` is the mean of `J`; see
  [`mixture2p()`](https://venpopov.com/bmm/dev/reference/mixture2p.md)
  for the parameterisation. The capacity models of Zhang and Luck (2008)
  do not apply to the IMM, which has no discrete slots – its bottleneck
  is the ratio of context to general activation – so variable precision
  is the only extension offered here.

- vp_nodes:

  Number of quadrature nodes used when `variable_precision = TRUE`; must
  be an odd number of at least 41, default 81. See
  [`mixture2p()`](https://venpopov.com/bmm/dev/reference/mixture2p.md)
  for the accuracy this buys and the shapes it refuses.

- links:

  A named list of link functions for the model parameters.

- ...:

  used internally for testing, ignore it

## Value

An object of class `bmmodel`

## Details

- **Domain:** Working memory

- **Task:** Continuous reproduction

- **Name:** Interference Measurement Model

- **Citation:**

  - Oberauer, K., & Lin, H.-Y. (2017). An interference model of visual
    working memory. Psychological Review, 124(1), 21-59.
    https://doi.org/10.1037/rev0000044

### Version: `full`

- **Requirements:**

  - The response variable should be in radians and represent the angular
    error relative to the target

  - The non-target features should be in radians and be centered
    relative to the target

- **Parameters:**

  - `mu`: Location of the von Mises distribution of memory responses (in
    radians). Fixed internally to 0 by default.

  - `kappa`: Concentration of the von Mises distribution of memory
    responses

  - `a`: General activation of memory items

  - `c`: Context activation

  - `s`: Spatial similarity gradient

  - `b`: Background activation, the reference against which the
    activations of the target and the non-targets are weighed. Fixed
    internally to 0 on its log scale, i.e. to an activation of 1;
    freeing it leaves the remaining activations unidentified, since only
    their ratio to the background enters the likelihood.

- **Fixed parameters:**

  - `mu` = 0

  - `b` = 0

- **Default parameter links:**

  - mu = tan_half; kappa = log; a = log; c = log; s = log; b = log

- **Default priors:**

  - `mu`:

    - `main`: normal(0, 0.5)

    - `effects`: normal(0, 0.25)

    - `sd`: exponential(4)

  - `kappa`:

    - `main`: normal(2, 1)

    - `effects`: normal(0, 1)

    - `sd`: exponential(1)

  - `a`:

    - `main`: normal(0, 1)

    - `effects`: normal(0, 1)

    - `sd`: exponential(1)

  - `c`:

    - `main`: normal(0, 1)

    - `effects`: normal(0, 1)

    - `sd`: exponential(1)

  - `s`:

    - `main`: normal(0, 1)

    - `effects`: normal(0, 1)

    - `sd`: exponential(1)

  - `b`:

    - `main`: normal(0, 1)

    - `effects`: normal(0, 1)

    - `sd`: exponential(1)

### Version: `bsc`

- **Requirements:**

  - The response variable should be in radians and represent the angular
    error relative to the target

  - The non-target features should be in radians and be centered
    relative to the target

- **Parameters:**

  - `mu`: Location of the von Mises distribution of memory responses (in
    radians). Fixed internally to 0 by default.

  - `kappa`: Concentration of the von Mises distribution of memory
    responses

  - `c`: Context activation

  - `s`: Spatial similarity gradient

  - `b`: Background activation, the reference against which the
    activations of the target and the non-targets are weighed. Fixed
    internally to 0 on its log scale, i.e. to an activation of 1;
    freeing it leaves the remaining activations unidentified, since only
    their ratio to the background enters the likelihood.

- **Fixed parameters:**

  - `mu` = 0

  - `b` = 0

- **Default parameter links:**

  - mu = tan_half; kappa = log; c = log; s = log; b = log

- **Default priors:**

  - `mu`:

    - `main`: normal(0, 0.5)

    - `effects`: normal(0, 0.25)

    - `sd`: exponential(4)

  - `kappa`:

    - `main`: normal(2, 1)

    - `effects`: normal(0, 1)

    - `sd`: exponential(1)

  - `c`:

    - `main`: normal(0, 1)

    - `effects`: normal(0, 1)

    - `sd`: exponential(1)

  - `s`:

    - `main`: normal(0, 1)

    - `effects`: normal(0, 1)

    - `sd`: exponential(1)

  - `b`:

    - `main`: normal(0, 1)

    - `effects`: normal(0, 1)

    - `sd`: exponential(1)

### Version: `abc`

- **Requirements:**

  - The response variable should be in radians and represent the angular
    error relative to the target

  - The non-target features should be in radians and be centered
    relative to the target

- **Parameters:**

  - `mu`: Location of the von Mises distribution of memory responses (in
    radians). Fixed internally to 0 by default.

  - `kappa`: Concentration of the von Mises distribution of memory
    responses

  - `a`: General activation of memory items

  - `c`: Context activation

  - `b`: Background activation, the reference against which the
    activations of the target and the non-targets are weighed. Fixed
    internally to 0 on its log scale, i.e. to an activation of 1;
    freeing it leaves the remaining activations unidentified, since only
    their ratio to the background enters the likelihood.

- **Fixed parameters:**

  - `mu` = 0

  - `b` = 0

- **Default parameter links:**

  - mu = tan_half; kappa = log; a = log; c = log; b = log

- **Default priors:**

  - `mu`:

    - `main`: normal(0, 0.5)

    - `effects`: normal(0, 0.25)

    - `sd`: exponential(4)

  - `kappa`:

    - `main`: normal(2, 1)

    - `effects`: normal(0, 1)

    - `sd`: exponential(1)

  - `a`:

    - `main`: normal(0, 1)

    - `effects`: normal(0, 1)

    - `sd`: exponential(1)

  - `c`:

    - `main`: normal(0, 1)

    - `effects`: normal(0, 1)

    - `sd`: exponential(1)

  - `b`:

    - `main`: normal(0, 1)

    - `effects`: normal(0, 1)

    - `sd`: exponential(1)

## Priors from bmm 1.4

Since bmm 1.5.0, the parameters of this model are distributional
parameters (`dpar`) of a custom family, and the location `mu1` is called
`mu`. Calls written for bmm 1.4 keep working until bmm 1.6.0 at the
earliest: bmm translates their priors, formulas, `links` and
`emmeans(nlpar = )`, and warns with the code to use instead.

    # an intercept prior
    set_prior("normal(2, 1)", class = "b", coef = "Intercept", nlpar = "kappa") # 1.4
    set_prior("normal(2, 1)", class = "Intercept", dpar = "kappa")              # 1.5

    # a prior on all coefficients
    set_prior("normal(0, 0.5)", class = "b", nlpar = "kappa")                   # 1.4
    set_prior("normal(0, 0.5)", class = "Intercept", dpar = "kappa") +          # 1.5
      set_prior("normal(0, 0.5)", class = "b", dpar = "kappa")

    # a prior on the random-effect standard deviations
    set_prior("exponential(2)", class = "sd", nlpar = "kappa")                  # 1.4
    set_prior("exponential(2)", class = "sd", dpar = "kappa")                   # 1.5

    # the location, and emmeans()
    bmf(mu1 ~ 1 + cond, kappa ~ 1 + cond)                                      # 1.4
    bmf(mu ~ 1 + cond, kappa ~ 1 + cond)                                       # 1.5
    emmeans(fit, ~ cond, nlpar = "kappa")                                      # 1.4
    emmeans(fit, ~ cond, dpar = "kappa")                                       # 1.5

In bmm 1.4 a `class = "b"` prior also covered the intercept, so it now
becomes two priors with the same distribution, unless the intercept has
a prior of its own. One prior changes its meaning. When a formula has an
intercept and predictors, such as `kappa ~ 1 + cond`, brms centres the
predictors of a distributional parameter, so an intercept prior applies
at the average over the conditions, weighted by their number of trials,
not at the reference condition. The slopes keep their meaning, and so do
intercept-only formulas and formulas with one coefficient per condition
(`kappa ~ 0 + cond`). If your intercept prior was meant for the
reference condition, give it for the average instead.

Two things are not translated. The draws of the location are named after
`mu` (`b_Intercept` instead of `b_mu1_Intercept`), and calls that go
straight to brms need `dpar =`:
`posterior_linpred(fit, dpar = "kappa")`, not `nlpar = "kappa"`. Fits
saved with bmm 1.4 work as before with
[`summary()`](https://rdrr.io/r/base/summary.html),
[`prior_info()`](https://venpopov.com/bmm/dev/reference/prior_info.md),
`emmeans()`,
[`native_parameters()`](https://venpopov.com/bmm/dev/reference/native_parameters.md)
and
[`conditional_effects()`](https://venpopov.com/bmm/dev/reference/conditional_effects.bmmfit.md).
[`update()`](https://rdrr.io/r/stats/update.html) samples them again
under the custom family, translates their priors and says so in a
message.

## References

Oberauer, K., & Lin, H.-Y. (2017). An interference model of visual
working memory. Psychological Review, 124(1), 21-59.

van den Berg, R., Shin, H., Chou, W.-C., George, R., & Ma, W. J. (2012).
Variability in encoding precision accounts for visual short-term memory
limitations. PNAS, 109(22), 8780-8785.

Zhang, W., & Luck, S. J. (2008). Discrete fixed-resolution
representations in visual working memory. Nature, 453(7192), 233-235.

## Examples

``` r
if (FALSE) { # identical(Sys.getenv("BMM_EXAMPLES"), "true")
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
