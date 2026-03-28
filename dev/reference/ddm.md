# Diffusion Decision Model

Diffusion Decision Model

## Usage

``` r
ddm(rt, response, links = NULL, ...)
```

## Arguments

- rt:

  Name of the reaction time variable coding reaction time in seconds in
  the data.

- response:

  Name of the response variable coding the response numerically (0 =
  lower response / incorrect, 1 = upper response / correct)

- links:

  A list of links for the parameters.

- ...:

  used internally for testing, ignore it

## Value

An object of class `bmmodel`

## Details

- **Domain:** Decision Making / Response times

- **Task:** Two-Alternative Force Choice RT

- **Name:** Diffusion Decision Model

- **Citation:**

  - Ratcliff, R. (1978). A theory of memory retrieval. Psychological
    Review, 85(2), 59-108. https://doi.org/10/fjwm2f;

- **Requirements:**

  - The response time should be in seconds and represent the time
    between onset of the target stimulus until the response execution

&nbsp;

- The response should be coded numerically: 0 = lower response, 1 =
  upper response

&nbsp;

- **Parameters:**

  - `drift`: Drift rate = Average rate of evidence accumulation of the
    decision processes

  - `bound`: Boundary separation = Distance between the decision
    boundaries that need to be reached

  - `ndt`: Non-decision time = Additional time required beyond the
    evidence accumulation process

  - `zr`: Relative starting point = Starting point between the decision
    thresholds relative to the upper bound.

- **Fixed parameters:**

  - `zr` = 0

  - `mu` = 0

- **Default parameter links:**

  - drift = identity; bound = log; ndt = log; zr = logit

- **Default priors:**

  - `drift`:

    - `main`: cauchy(0,1)

    - `effects`: normal(0,0.5)

  - `bound`:

    - `main`: normal(0,0.5)

    - `effects`: normal(0,0.5)

  - `ndt`:

    - `main`: normal(-1.5,0.5)

    - `effects`: normal(0,0.3)

  - `zr`:

    - `main`: normal(0,0.5)

    - `effects`: normal(0,0.3)

## Default behavior

By default, `zr` is fixed at 0. If you want to estimate `zr`, add a
formula for `zr` in your
[`bmf()`](https://venpopov.com/bmm/dev/reference/bmmformula.md) call.

## Examples

``` r
if (FALSE) { # isTRUE(Sys.getenv("BMM_EXAMPLES"))
# Simulate data for one subject using rddm
set.seed(123)
n_trials <- 500

# Simulate DDM data with fixed parameters
sim_data <- rddm(
  n = n_trials,
  drift = 1.5,    # drift rate
  bound = 1.2,    # boundary separation
  ndt = 0.3,      # non-decision time
  zr = 0.5        # relative starting point
)

# Prepare data frame
dat <- data.frame(
  rt = sim_data$rt,
  response = sim_data$response
)

# Define formula (intercept-only model)
ff <- bmmformula(
  drift ~ 1,
  bound ~ 1,
  ndt ~ 1
)

# Specify the DDM model
model <- ddm(rt = "rt", response = "response")

# Fit the model
fit <- bmm(
  formula = ff,
  data = dat,
  model = model,
  cores = 4,
  iter = 1000,
  backend = "cmdstanr"
)

# Check parameter recovery
summary(fit)
}
```
