# EZ-Diffusion Model

EZ-Diffusion Model

## Usage

``` r
ezdm(mean_rt, var_rt, n_upper, n_trials, links = NULL, version = "3par", ...)
```

## Arguments

- mean_rt:

  The names of the variable or variables (for 4par version) coding the
  mean reaction time in seconds in the data.

- var_rt:

  The names of the variable or variables (for 4par version) coding the
  variance of the reaction time in seconds in the data

- n_upper:

  The name of the variable coding the number of responses that hit the
  upper response threshold (typically the number of correct responses)
  in the data.

- n_trials:

  The name of the variable coding the number of trials that was used to
  calculated the aggregated statistics.

- links:

  A list of links for the parameters.

- version:

  A character label for the version of the model. There is a
  three-parameter version (version = "3par") of the `ezdm` that fixes
  the relative starting point `zr` to 0.5, and a four parameter version
  (version = "4par"), that allows to freely estimate the starting point.

- ...:

  used internally for testing, ignore it

## Value

An object of class `bmmodel`

## Details

- **Domain:** Processing Speed, Decision Making

- **Task:** Choice Reaction Time tasks

- **Name:** EZ-Diffusion Model

- **Citation:**

  - Wagenmakers, E.-J., Van Der Maas, H. L. J., & Grasman, R. P. P. P.
    (2007). An EZ-diffusion model for response time and accuracy.
    Psychonomic Bulletin & Review, 14(1), 3-22.
    https://doi.org/10/fk447c

&nbsp;

- Chávez De la Peña, A. F., & Vandekerckhove, J. (2025). An EZ Bayesian
  hierarchical drift diffusion model for response time and accuracy.
  Psychonomic Bulletin & Review.
  https://doi.org/10.3758/s13423-025-02729-y

&nbsp;

- **Version:** 4par

- **Requirements:**

  Provide aggregated statistics for each subject and condition that
  model parameters should vary over:

&nbsp;

- Mean reaction times (mean_rt) in seconds

- Variance of reaction times (var_rt) in seconds

- Number of responses to the upper decision threshold (n_upper)

- Total number of trials used to calculate aggregated statistics
  (n_trials)

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

  - `s`: The diffusion constant, that is the standard deviation of the
    Gaussian noise during sampling

- **Fixed parameters:**

  - `s` = 0

  - `mu` = 0

- **Default parameter links:**

  - drift = identity; bound = log; ndt = log; zr = logit; s = log

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

  - `s`:

    - `main`: normal(0,1)

    - `effects`: normal(0,0.3)

## Examples

``` r
if (FALSE) { # \dontrun{
# Minimal parameter recovery example with 3-parameter EZDM

# Simulate data from known parameters
set.seed(123)
sim_data <- rezdm(
  n = 10,
  n_trials = 100,
  drift = 2,
  bound = 1.5,
  ndt = 0.3,
  version = "3par"
)

# Add subject ID
sim_data$id <- 1:10

# Specify model
model <- ezdm(
  mean_rt = "mean_rt",
  var_rt = "var_rt",
  n_upper = "n_upper",
  n_trials = "n_trials",
  version = "3par"
)

# Specify formula with random effects
formula <- bmf(
  drift ~ 1 + (1 | id),
  bound ~ 1 + (1 | id),
  ndt ~ 1
)

# Fit model (using cmdstanr backend)
fit <- bmm(
  formula = formula,
  data = sim_data,
  model = model,
  backend = "cmdstanr",
  cores = 4,
  chains = 4,
  iter = 2000,
  warmup = 1000
)

# Check parameter recovery
summary(fit)

# Extract population-level effects
# True values: drift = 2, bound = 1.5, ndt = 0.3 (on log scale for drift/bound)
exp(brms::fixef(fit))
} # }
```
