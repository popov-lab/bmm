# Three-Parameter Mixture Model

Three-Parameter Mixture Model

## Usage

``` r
mixture3p(
  resp_error,
  nt_features,
  set_size,
  regex = FALSE,
  version = c("simple", "slot", "slot_averaging"),
  variable_precision = FALSE,
  vp_nodes = 81L,
  links = NULL,
  ...
)
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
  In the `"slot"` and `"slot_averaging"` versions a single set size
  cannot test the capacity rule: `K` is then estimated from the
  responses at that set size, and only below it.

- regex:

  Logical. If TRUE, the `nt_features` argument is interpreted as a
  regular expression to match the non-target feature columns in the
  dataset.

- version:

  A character string specifying which version of the model to use.
  Options are:

  - `"simple"` (default): the standard three-parameter mixture model.
    The weights of the target (`thetat`), the non-targets (`thetant`)
    and guessing are normalised across the components a trial actually
    has, with guessing as the reference.

  - `"slot"`: the probability that a response comes from memory follows
    the fixed-resolution slot model of Zhang and Luck (2008),
    `min(1, K / set_size)`, and `ptarget` is the probability that a
    remembered response reports the target rather than a non-target.

  - `"slot_averaging"`: as `"slot"`, but `K` slots are distributed over
    the items and precision grows with the number of slots an item
    holds.

  Storage and swapping are separate in the capacity versions: the
  capacity rule says how often a response comes from memory and
  `ptarget` says which stored item is reported. In `"simple"` both
  weights are log odds against guessing, `thetat` for the target and
  `thetant` for all non-targets together, so neither is a probability on
  its own;
  [`native_parameters()`](https://venpopov.com/bmm/dev/reference/native_parameters.md)
  returns the probabilities.

- variable_precision:

  Logical; if `TRUE`, the precision of memory varies from trial to trial
  (van den Berg et al., 2012). This adds the parameter `tau`, and
  `kappa` then describes the average trial: it is the concentration
  whose Fisher information `J(kappa)` is the mean of `J`. See
  [`mixture2p()`](https://venpopov.com/bmm/dev/reference/mixture2p.md)
  for details of the parameterisation. Under
  `version = "slot_averaging"`, the precision of an item held in `j`
  slots is the sum of `j` independent slot draws,
  `gamma(shape = j * J(kappa) / tau, scale = tau)`, so its variance
  grows with `j`.

- vp_nodes:

  Number of quadrature nodes used when `variable_precision = TRUE`; must
  be an odd number of at least 41, and defaults to 81. See
  [`mixture2p()`](https://venpopov.com/bmm/dev/reference/mixture2p.md)
  for the accuracy this gives.

- links:

  A named list of link functions for the model parameters. The links of
  `thetat` and `thetant` cannot be changed, because the softmax over
  them is part of the likelihood.

- ...:

  used internally for testing, ignore it

## Value

An object of class `bmmodel`

## Details

- **Domain:** Working memory

- **Task:** Continuous reproduction

- **Name:** Three-Parameter Mixture Model

- **Citation:**

  - Bays, P. M., Catalao, R. F. G., & Husain, M. (2009). The precision
    of visual working memory is set by allocation of a shared resource.
    Journal of Vision, 9(10), Article 7. https://doi.org/10.1167/9.10.7

  - Zhang, W., & Luck, S. J. (2008). Discrete fixed-resolution
    representations in visual working memory. Nature, 453(7192),
    233-235. https://doi.org/10.1038/nature06860

  - van den Berg, R., Shin, H., Chou, W.-C., George, R., & Ma, W. J.
    (2012). Variability in encoding precision accounts for visual
    short-term memory limitations. Proceedings of the National Academy
    of Sciences, 109(22), 8780-8785.
    https://doi.org/10.1073/pnas.1117465109

- **Requirements:**

  - The response variable should be in radians and represent the angular
    error relative to the target

  - The non-target features should be in radians and be centered
    relative to the target

  - In the 'slot' and 'slot_averaging' versions, the set size should
    vary for the capacity rule to be testable

### Version: `simple`

- **Parameters:**

  - `mu`: Location of the von Mises distribution of memory responses (in
    radians). Fixed internally to 0 by default.

  - `kappa`: Concentration of the von Mises distribution of memory
    responses

  - `thetat`: Mixture weight for target responses. The weights of the
    target, the non-targets and guessing are normalised across the
    components present on a trial, with guessing as the reference, so
    thetat maps to the probability of a target response only jointly
    with thetant.

  - `thetant`: Mixture weight for non-target responses, summed over the
    non-targets, which share it equally. Like thetat, it is weighed
    against guessing.

- **Fixed parameters:**

  - `mu` = 0

- **Default parameter links:**

  - mu = tan_half; kappa = log; thetat = softmax; thetant = softmax

- **Default priors:**

  - `mu`:

    - `main`: normal(0, 0.5)

    - `effects`: normal(0, 0.25)

    - `sd`: exponential(4)

  - `kappa`:

    - `main`: normal(2, 1)

    - `effects`: normal(0, 1)

    - `sd`: exponential(1)

  - `thetat`:

    - `main`: logistic(0, 1)

    - `effects`: normal(0, 2.5)

    - `sd`: exponential(1)

  - `thetant`:

    - `main`: logistic(0, 1)

    - `effects`: normal(0, 2.5)

    - `sd`: exponential(1)

### Version: `slot`

- **Parameters:**

  - `mu`: Location of the von Mises distribution of memory responses (in
    radians). Fixed internally to 0 by default.

  - `kappa`: Concentration of the von Mises distribution of memory
    responses

  - `K`: Capacity: the number of items that can be held in memory. An
    item is either held at full precision or not held at all, so the
    probability that a response comes from memory is min(1, K /
    set_size).

  - `ptarget`: Probability that a response reports the target rather
    than a non-target, given that it comes from memory.

- **Fixed parameters:**

  - `mu` = 0

- **Default parameter links:**

  - mu = tan_half; kappa = log; K = log; ptarget = logit

- **Default priors:**

  - `mu`:

    - `main`: normal(0, 0.5)

    - `effects`: normal(0, 0.25)

    - `sd`: exponential(4)

  - `kappa`:

    - `main`: normal(2, 1)

    - `effects`: normal(0, 1)

    - `sd`: exponential(1)

  - `K`:

    - `main`: normal(1, 0.5)

    - `effects`: normal(0, 0.3)

    - `sd`: exponential(2)

  - `ptarget`:

    - `main`: logistic(1, 1)

    - `effects`: normal(0, 0.5)

    - `sd`: exponential(1)

### Version: `slot_averaging`

- **Parameters:**

  - `mu`: Location of the von Mises distribution of memory responses (in
    radians). Fixed internally to 0 by default.

  - `kappa`: Concentration of the von Mises distribution for an item
    held in a single slot. An item holding several slots is more
    precise, because averaging independent samples adds their Fisher
    information.

  - `K`: Capacity: the number of memory slots distributed over the items
    in the array. An item receives floor(K / set_size) or one more slot,
    and an item that receives none is guessed.

  - `ptarget`: Probability that a response reports the target rather
    than a non-target, given that it comes from memory.

- **Fixed parameters:**

  - `mu` = 0

- **Default parameter links:**

  - mu = tan_half; kappa = log; K = log; ptarget = logit

- **Default priors:**

  - `mu`:

    - `main`: normal(0, 0.5)

    - `effects`: normal(0, 0.25)

    - `sd`: exponential(4)

  - `kappa`:

    - `main`: normal(2, 1)

    - `effects`: normal(0, 1)

    - `sd`: exponential(1)

  - `K`:

    - `main`: normal(1, 0.5)

    - `effects`: normal(0, 0.3)

    - `sd`: exponential(2)

  - `ptarget`:

    - `main`: logistic(1, 1)

    - `effects`: normal(0, 0.5)

    - `sd`: exponential(1)

### Variable precision

With `variable_precision = TRUE` every version gains `tau`, and `kappa`
becomes the precision of the average trial. Shown for `"simple"`:

- **Parameters:**

  - `mu`: Location of the von Mises distribution of memory responses (in
    radians). Fixed internally to 0 by default.

  - `kappa`: Concentration of the von Mises distribution of memory
    responses. With variable precision, kappa describes the average
    trial: it is the concentration whose Fisher information J(kappa) is
    the mean of J across trials.

  - `tau`: Scale of the trial-to-trial variability in memory precision.
    The Fisher information J of a memory representation is drawn from
    gamma(mean = J(kappa), scale = tau), so tau -\> 0 recovers constant
    precision and larger values spread precision more widely across
    trials.

  - `thetat`: Mixture weight for target responses. The weights of the
    target, the non-targets and guessing are normalised across the
    components present on a trial, with guessing as the reference, so
    thetat maps to the probability of a target response only jointly
    with thetant.

  - `thetant`: Mixture weight for non-target responses, summed over the
    non-targets, which share it equally. Like thetat, it is weighed
    against guessing.

- **Default parameter links:**

  - mu = tan_half; kappa = log; tau = log; thetat = softmax; thetant =
    softmax

- **Default priors:**

  - `mu`:

    - `main`: normal(0, 0.5)

    - `effects`: normal(0, 0.25)

    - `sd`: exponential(4)

  - `kappa`:

    - `main`: normal(2, 1)

    - `effects`: normal(0, 1)

    - `sd`: exponential(1)

  - `thetat`:

    - `main`: logistic(0, 1)

    - `effects`: normal(0, 2.5)

    - `sd`: exponential(1)

  - `thetant`:

    - `main`: logistic(0, 1)

    - `effects`: normal(0, 2.5)

    - `sd`: exponential(1)

  - `tau`:

    - `main`: normal(0, 1)

    - `effects`: normal(0, 0.5)

    - `sd`: exponential(2)

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

Bays, P. M., Catalao, R. F. G., & Husain, M. (2009). The precision of
visual working memory is set by allocation of a shared resource. Journal
of Vision, 9(10), Article 7.

Zhang, W., & Luck, S. J. (2008). Discrete fixed-resolution
representations in visual working memory. Nature, 453(7192), 233-235.

van den Berg, R., Shin, H., Chou, W.-C., George, R., & Ma, W. J. (2012).
Variability in encoding precision accounts for visual short-term memory
limitations. PNAS, 109(22), 8780-8785.

## Examples

``` r
if (FALSE) { # identical(Sys.getenv("BMM_EXAMPLES"), "true")
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
