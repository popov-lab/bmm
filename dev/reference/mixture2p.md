# Two-Parameter Mixture Model

Two-Parameter Mixture Model

## Usage

``` r
mixture2p(
  resp_error,
  set_size = NULL,
  version = c("simple", "slot", "slot_averaging"),
  variable_precision = FALSE,
  vp_nodes = 81L,
  links = NULL,
  ...
)
```

## Arguments

- resp_error:

  The name of the variable in the provided dataset containing the
  response error. The response Error should code the response relative
  to the to-be-recalled target in radians. You can transform the
  response error in degrees to radian using the `deg2rad` function.

- set_size:

  Name of the variable in the dataset containing the set size, or a
  single number if the set size is constant. Required by the `"slot"`
  and `"slot_averaging"` versions and ignored by `"simple"`. With a
  single set size the capacity rule cannot be tested: `K` is then
  estimated from the responses at that set size, and only below it.

- version:

  A character string specifying which version of the model to use.
  Options are:

  - `"simple"` (default): the standard two-parameter mixture model, in
    which the probability of a memory response (`thetat`) and the
    precision of memory (`kappa`) are free parameters.

  - `"slot"`: the fixed-resolution slot model of Zhang and Luck (2008).
    An item is either held at full precision or not held at all, so the
    probability of a memory response is `min(1, K / set_size)` and
    `kappa` does not vary with set size. Note that `K` is not identified
    above `max(set_size)`, where every item is remembered whatever `K`
    is, and that the cap introduces a kink in the gradient at
    `K = set_size`.

  - `"slot_averaging"`: the slots-plus-averaging model of Zhang and Luck
    (2008). `K` slots are distributed over the items, an item receives
    `floor(K / set_size)` or one more, and precision grows with the
    number of slots because averaging independent samples adds their
    Fisher information.

- variable_precision:

  Logical; if `TRUE`, the precision of memory varies from trial to trial
  rather than being constant (van den Berg et al., 2012; Fougnie et al.,
  2012). This adds the parameter `tau`, and the likelihood integrates
  the Fisher information `J` of a memory representation over
  `gamma(mean = J(kappa), scale = tau)`. `kappa` then describes the
  average trial: it is the concentration whose Fisher information
  `J(kappa)` is the mean of `J`, and `tau -> 0` recovers the
  constant-precision model. Under `version = "slot_averaging"`, the
  precision of an item held in `j` slots is the sum of `j` independent
  slot draws, `gamma(shape = j * J(kappa) / tau, scale = tau)`, so its
  variance grows with `j`. The integral is evaluated by quadrature,
  which costs about ten times as much per gradient at the default
  `vp_nodes`; expect a larger factor in wall time, since the posterior
  geometry also needs more leapfrog steps. `tau` is the parameter that
  needs the most data, and is often only weakly identified in small
  designs.

- vp_nodes:

  Number of quadrature nodes used when `variable_precision = TRUE`; must
  be an odd number of at least 41. With the default of 81, the log
  likelihood holds to about 1e-5 as long as the implied gamma shape
  `J(kappa) / tau` stays above `60 / vp_nodes` (0.74 at the default),
  and to about 2e-6 above a shape of 1. Below that threshold the model
  rejects rather than returning a biased likelihood, and asks for more
  nodes; 121 nodes reach a shape of about 0.5, and 241 about 0.25.

- links:

  A named list of link functions for the model parameters, e.g.
  `list(kappa = "softplus")`. Defaults are `"log"` for `kappa`, `K` and
  `tau`, and `"logit"` for `thetat`.

- ...:

  used internally for testing, ignore it

## Value

An object of class `bmmodel`

## Details

- **Domain:** Working memory

- **Task:** Continuous reproduction

- **Name:** Two-Parameter Mixture Model

- **Citation:**

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

  - The 'slot' and 'slot_averaging' versions additionally require a
    set_size variable, which should vary for the capacity rule to be
    testable

### Version: `simple`

- **Parameters:**

  - `mu`: Location of the von Mises distribution of memory responses (in
    radians). Fixed internally to 0 by default.

  - `kappa`: Concentration of the von Mises distribution of memory
    responses

  - `thetat`: Mixture weight for target responses, i.e. the probability
    that a response comes from memory rather than from guessing.

- **Fixed parameters:**

  - `mu` = 0

- **Default parameter links:**

  - mu = tan_half; kappa = log; thetat = logit

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

### Version: `slot`

- **Parameters:**

  - `mu`: Location of the von Mises distribution of memory responses (in
    radians). Fixed internally to 0 by default.

  - `kappa`: Concentration of the von Mises distribution of memory
    responses

  - `K`: Capacity: the number of items that can be held in memory. An
    item is either held at full precision or not held at all, so the
    probability of a memory response is min(1, K / set_size).

- **Fixed parameters:**

  - `mu` = 0

- **Default parameter links:**

  - mu = tan_half; kappa = log; K = log

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

- **Fixed parameters:**

  - `mu` = 0

- **Default parameter links:**

  - mu = tan_half; kappa = log; K = log

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

  - `thetat`: Mixture weight for target responses, i.e. the probability
    that a response comes from memory rather than from guessing.

- **Default parameter links:**

  - mu = tan_half; kappa = log; tau = log; thetat = logit

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

Zhang, W., & Luck, S. J. (2008). Discrete fixed-resolution
representations in visual working memory. Nature, 453(7192), 233-235.

van den Berg, R., Shin, H., Chou, W.-C., George, R., & Ma, W. J. (2012).
Variability in encoding precision accounts for visual short-term memory
limitations. PNAS, 109(22), 8780-8785.

Fougnie, D., Suchow, J. W., & Alvarez, G. A. (2012). Variability in the
quality of visual working memory. Nature Communications, 3(1), 1229.

## Examples

``` r
if (FALSE) { # identical(Sys.getenv("BMM_EXAMPLES"), "true")
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
