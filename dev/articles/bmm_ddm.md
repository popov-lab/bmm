# The Diffusion Decision Model (ddm)

## 1 Introduction

The Diffusion Decision Model (DDM) is one of the most widely used
cognitive models for analyzing response time (RT) and accuracy data in
two-alternative forced-choice (2AFC) tasks. Originally developed by
Ratcliff (1978), the DDM has become a cornerstone of cognitive
psychology and neuroscience research on decision-making processes.

The DDM is ideal when:

- You have RT data from a two-choice decision task with both correct and
  error responses
- You want to decompose decision-making into separate cognitive
  processes
- You need to account for speed-accuracy tradeoffs
- You want to model trial-to-trial variability in decision processes
- You’re interested in linking cognitive parameters to neural or
  behavioral measures

The `bmm` package provides a single DDM specification with a default
fixed starting point. If you include a formula for `zr` in your
[`bmf()`](https://venpopov.com/bmm/dev/reference/bmmformula.md) call,
the starting point becomes a freely estimated parameter. Trial-to-Trial
variability parameters of the DDM are currently not supported, due to
feasibility for Bayesian hierarchical estimation.

## 2 Theoretical Background

### 2.1 The Diffusion Process

The DDM models decision-making as a noisy accumulation of evidence over
time toward one of two decision boundaries. The evidence accumulation
process is described by the Wiener diffusion process:

\\ dX(t) = v \cdot dt + s \cdot dW(t) \\

where:

- \\X(t)\\ is the accumulated evidence at time \\t\\
- \\v\\ is the drift rate (average rate of evidence accumulation)
- \\s\\ is the within-trial noise (standard deviation of the diffusion
  process, typically fixed at 0.1 or 1)
- \\dW(t)\\ is a Wiener process (Brownian motion)

The process starts at position \\z \cdot a\\ (where \\z\\ is the
starting point parameter and \\a\\ is the boundary separation) and
continues until it reaches one of two absorbing boundaries at 0
(lower/incorrect) or \\a\\ (upper/correct).

### 2.2 Model Parameters

#### 2.2.1 Core Parameters

1.  **Drift rate (\\v\\)**: The average rate at which evidence
    accumulates toward the correct boundary. Higher drift rates lead to
    faster and more accurate decisions.
    - Typical range: -5 to 5 (task-dependent)
    - Interpretation: Quality of evidence, task difficulty, or attention
2.  **Boundary separation (\\a\\)**: The amount of evidence required
    before making a decision. Larger boundaries lead to slower but more
    accurate responses.
    - Typical range: 0.5 to 3
    - Interpretation: Response caution, speed-accuracy tradeoff
3.  **Non-decision time (\\t_0\\ or ndt)**: Time for processes unrelated
    to the decision (e.g., stimulus encoding, motor response).
    - Typical range: 0.1 to 0.5 seconds
    - Interpretation: Perceptual and motor latencies

#### 2.2.2 Starting Point Bias

4.  **Starting point bias (\\z_r\\)**: Relative starting point of
    evidence accumulation (0 = lower boundary, 1 = upper boundary).

- Default: fixed at 0 (which implies \\z_r = 0.5\\ on the response
  scale)
- Interpretation: A priori response bias

#### 2.2.3 Variability Parameters (7-parameter DDM)

The full DDM, as described by Ratcliff and McKoon (2008) and implemented
by Henrich et al. (2024) in STAN, includes three additional parameters
that capture trial-to-trial variability in the core parameters:

5.  **Drift rate variability (\\s_v\\ or sdrift)**: Between-trial
    variability in drift rate, accounting for fluctuations in evidence
    quality across trials (typical range: 0 to 2).

6.  **Non-decision time variability (\\s\_{t_0}\\ or sndt)**:
    Between-trial variability in non-decision time, accounting for
    variable encoding/motor times (typical range: 0 to 0.2).

7.  **Starting point variability (\\s_z\\ or szr)**: Between-trial
    variability in starting point, accounting for variable response bias
    across trials (typical range: 0 to 0.5).

**Why the 7-parameter version is not implemented in `bmm`**:

While these variability parameters are theoretically important and can
improve model fit, Henrich et al. (2024) demonstrated that the
7-parameter likelihood has prohibitively long computation times,
especially for hierarchical Bayesian estimation. Key limitations
include:

- **Single-subject feasibility only**: Even with 100-300 trials, model
  fitting can up to 20 minutes when using within-chain parallelization.
  With 1500+ trials typical for robust DDM estimation, runtimes likely
  extend to hours or even days.

- **Hierarchical models infeasible**: The most common use case in
  cognitive research (multiple subjects, conditions, and trials) becomes
  computationally intractable without high-performance computing
  resources. A typical hierarchical study (20+ subjects, 2-3 conditions,
  200+ trials per condition) would likely require weeks of computation
  time.

- **Limited practical benefit**: For most research questions, fixed- and
  free-starting-point models provide adequate fits and interpretable
  parameters, making the substantial computational cost of the
  7-parameter version difficult to justify.

For researchers specifically interested in trial-to-trial variability
parameters, we recommend considering simulation-based inference methods
(e.g., BayesFlow) or fitting separate single-subject models rather than
hierarchical approaches.

### 2.3 Predictions

The DDM predicts complete RT distributions for both correct and error
responses, including:

- Mean RT and accuracy
- RT variability (standard deviation)
- Shape of RT distributions (skewness, leading edge)
- Speed-accuracy tradeoffs
- Effects of experimental manipulations on different cognitive processes

## 3 Implementation in bmm

### 3.1 Computational Efficiency Considerations

The DDM likelihood requires numerical solution of the Fokker-Planck
equation, which can be computationally demanding for large datasets.
Runtime scales approximately linearly with the number of trials. For
typical datasets (100-300 trials), models fit quickly (\< 20 seconds on
modern hardware with 4 chains × 2000 iterations). Larger datasets (1500+
trials) may take several minutes for single-subject models and several
hours for hierarchical models with many subjects.

**Practical recommendations**: Start with the default model (fixed `zr`)
to test your model specification. Use the default `cores = 4` to
parallelize chains efficiently. For hierarchical models, consider
starting with random intercepts only before adding random slopes. Test
on a subset of data first to verify model specification and estimate
total runtime.

### 3.2 Example Dataset: Ratcliff & Rouder (1998)

The `rr98` dataset from the `rtdists` package provides a classic example
of choice response time data. The dataset contains responses from a
brightness discrimination task where participants judged whether pixel
arrays were “high” or “low” brightness.

``` r

library(bmm)
library(ggplot2)
library(dplyr)

# Load the dataset from rtdists
data(rr98, package = "rtdists")

# Remove outliers as in the original analysis
rr98 <- rr98[!rr98$outlier, ]

# Examine structure
str(rr98[, c("id", "rt", "response_num", "instruction", "strength")])

# Check accuracy by condition
rr98 %>%
  mutate(correct = (response_num == 1 & strength <= 16) | 
                   (response_num == 2 & strength > 16)) %>%
  group_by(id, instruction) %>%
  summarise(accuracy = mean(correct), .groups = "drop")
```

The dataset includes:

- **rt**: Response time in seconds
- **response_num**: Binary response (1 = “dark”, 2 = “light”)
- **instruction**: Speed vs accuracy emphasis
- **strength**: Task difficulty (0-32, brightness level)
- **id**: Participant identifier

This dataset is ideal for DDM analysis because:

1.  **Sufficient error trials**: Unlike tasks with near-ceiling accuracy
2.  **Multiple conditions**: Speed vs accuracy emphasis
3.  **Clear manipulation**: Brightness strength affects evidence quality
4.  **Rich structure**: Individual differences and trial-level
    variability

### 3.3 Preparing Data for DDM

The DDM requires:

1.  A numeric RT variable (positive values)
2.  A binary response variable coded as 0/1 or “lower”/“upper”

``` r

# Prepare data for DDM analysis
# Focus on one participant and instruction condition for demonstration
demo_data <- rr98 %>%
  dplyr::filter(id == "jf", instruction == "accuracy") %>%
  mutate(
    # Recode response: 1 = "upper", 0 = "lower"
    response = ifelse(response_num == 1, 0, 1)
  )

# Check response coding
table(demo_data$response)

# Visualize RT distributions by response
ggplot(demo_data, aes(x = rt, fill = factor(response))) +
  facet_wrap(~ strength) +
  geom_histogram(alpha = 0.6, bins = 30, position = "identity") +
  labs(x = "Response Time (s)", y = "Count", fill = "Response") +
  theme_minimal()
```

## 4 Fitting DDM Models in bmm

### 4.1 DDM with Fixed Starting Point (Default)

The simplest DDM with drift, boundary, and non-decision time:

``` r

# Define the model
model_fixed_zr <- ddm(rt = "rt", 
                      response = "response")

# Define the formula - let drift vary by strength
formula_fixed_zr <- bmf(
  drift ~ 0 + strength,
  bound ~ 1,
  ndt ~ 1
)

# Fit the model 
fit_fixed_zr <- bmm(
  formula = formula_fixed_zr,
  data = demo_data,
  model = model_fixed_zr,
  backend = "cmdstanr",
  cores = 4,
  chains = 4,
  iter = 2000,
  warmup = 1000,
  file = "assets/bmmfit_ddm_fixed_zr_vignette"
)

# Examine convergence
summary(fit_fixed_zr)

# Check diagnostics
plot(fit_fixed_zr)
```

**Interpreting the fixed-`zr` model**:

- **drift**: Higher values indicate faster evidence accumulation (easier
  trials)
- **bound**: Higher values indicate more cautious responding
- **ndt**: Time spent on non-decision processes

### 4.2 DDM with Free Starting Point Bias

Adding starting point bias allows modeling a priori response
preferences:

``` r

# Define the model with free `zr`
model_free_zr <- ddm(rt = "rt", 
                     response = "response")

# Define the formula
formula_free_zr <- bmf(
  drift ~ 0 + strength,
  bound ~ 1,
  ndt ~ 1,
  zr ~ 1  # starting point bias
)

# Fit the model
fit_free_zr <- bmm(
  formula = formula_free_zr,
  data = demo_data,
  model = model_free_zr,
  backend = "cmdstanr",
  cores = 4,
  file = "assets/bmmfit_ddm_free_zr_vignette"
)

# Compare models
loo_fixed_zr <- loo(fit_fixed_zr)
loo_free_zr <- loo(fit_free_zr)
loo_compare(loo_fixed_zr, loo_free_zr)
```

**When to use starting point bias**:

- When you suspect a priori response preferences
- When conditions systematically favor one response
- When testing theoretical predictions about response bias

## 5 Hierarchical DDM Models

One of the key advantages of fitting DDMs in a Bayesian framework with
`bmm` is the ability to easily specify hierarchical models that account
for individual differences.

### 5.1 Multi-Subject Analysis

``` r

# Use multiple subjects from rr98
multi_subject_data <- rr98 %>%
  filter(instruction == "accuracy") %>%
  mutate(response = ifelse(response_num == 1, 0, 1)) %>%
  select(id, rt, response, strength)

# Hierarchical model
model_hier <- ddm(rt = "rt", 
                  response = "response")

# Formula with random effects
formula_hier <- bmf(
  drift ~ 0 + strength + (0 + strength | id),
  bound ~ 1 + (1 | id),
  ndt ~ 1 + (1 | id)
)

# Fit hierarchical model
fit_hier <- bmm(
  formula = formula_hier,
  data = multi_subject_data,
  model = model_hier,
  backend = "cmdstanr",
  cores = 4
)

# Extract individual-level parameters
# This gives you person-specific drift, bound, and ndt estimates
```

**Benefits of hierarchical modeling**:

1.  **Partial pooling**: Borrows strength across subjects for better
    estimates
2.  **Individual differences**: Quantifies between-subject variability
3.  **Handles unbalanced data**: Works even with different trial counts
    per subject
4.  **Shrinkage**: Extreme individual estimates are pulled toward group
    mean

### 5.2 Condition Effects

Test how experimental manipulations affect different cognitive
processes:

``` r

# Compare speed vs accuracy emphasis
full_data <- rr98 %>%
  mutate(response = ifelse(response_num == 1, 0, 1)) %>%
  select(id, rt, response, instruction, strength)

# Model with instruction effects
formula_instruction <- bmf(
  drift ~ 0 + strength + (0 + strength | id),
  bound ~ instruction + (instruction | id),
  ndt ~ 1 + (1 | id)
)

fit_instruction <- bmm(
  formula = formula_instruction,
  data = full_data,
  model = model_hier,
  backend = "cmdstanr",
  cores = 4
)
```

**Theoretical predictions**:

- **Speed emphasis**: Lower boundary (faster, less accurate)
- **Accuracy emphasis**: Higher boundary (slower, more accurate)
- **Drift**: Should primarily reflect stimulus strength, not instruction

## 6 Model Comparison and Selection

### 6.1 Comparing Fixed vs Free `zr`

Use information criteria to compare model fit:

``` r

# Compute LOO for each model
loo_fixed_zr <- loo(fit_fixed_zr)
loo_free_zr <- loo(fit_free_zr)

# Compare models
loo_compare(loo_fixed_zr, loo_free_zr)

# Check for influential observations
plot(loo_free_zr)
```

**General findings** (from simulation studies):

- Fixed-`zr` model often sufficient for most applications
- Free-`zr` model improves fit when response bias is present
- Model selection should balance fit quality with interpretability

### 6.2 Posterior Predictive Checks

Assess model adequacy by comparing predictions to data:

``` r

# Generate posterior predictions
pp_fixed_zr <- posterior_predict(fit_fixed_zr)

# Compare predicted vs observed RT distributions
library(bayesplot)
ppc_dens_overlay(demo_data$rt, pp_fixed_zr[1:50, ])

# Quantile-probability plots (QP plots)
# These are the gold standard for evaluating DDM fits
# Shows predicted vs observed RT quantiles for correct/error responses
```

## 7 Hypothesis Testing

### 7.1 Testing Parameter Differences

Use Bayesian hypothesis testing to compare parameters across conditions:

``` r

# Test if drift differs between conditions
hypothesis(fit_instruction, "driftinstructionspeed < 0")

# Test if boundary differs between conditions  
hypothesis(fit_instruction, "boundinstructionspeed < boundinstructionaccuracy")

# Joint hypothesis about multiple parameters
hypothesis(fit_instruction, c(
  "driftstrength > 0",
  "boundinstructionaccuracy > boundinstructionspeed"
))
```

### 7.2 Bayes Factors

Compare specific hypotheses using bridge sampling:

``` r

# Fit constrained model (drift doesn't vary by strength)
formula_null <- bmf(
  drift ~ 1 + (1 | id),
  bound ~ instruction + (instruction | id),
  ndt ~ 1 + (1 | id)
)

fit_null <- bmm(
  formula = formula_null,
  data = full_data,
  model = model_hier,
  backend = "cmdstanr"
)

# Compute Bayes factor (requires bridgesampling package)
# bf <- bayes_factor(fit_instruction, fit_null)
```

## 8 Advanced Topics

### 8.1 Custom Priors

Specify informative priors based on prior research or pilot data:

``` r

# Define custom priors
custom_priors <- c(
  prior(normal(2, 1), nlpar = "drift", coef = "Intercept"),
  prior(normal(1, 0.5), nlpar = "bound", coef = "Intercept"),
  prior(normal(-1, 0.3), nlpar = "ndt", coef = "Intercept")  # on log scale
)

# Fit with custom priors
fit_custom <- bmm(
  formula = formula_fixed_zr,
  data = demo_data,
  model = model_fixed_zr,
  prior = custom_priors,
  backend = "cmdstanr"
)
```

### 8.2 Fixing Parameters

Fix parameters to specific values for hypothesis testing:

``` r

# Fix starting point to 0.5 (unbiased) while estimating others
formula_fixed <- bmf(
  drift ~ strength,
  bound ~ 1,
  ndt ~ 1,
  zr = 0.5  # fixed parameter
)

fit_fixed <- bmm(
  formula = formula_fixed,
  data = demo_data,
  model = model_fixed_zr,
  backend = "cmdstanr"
)
```

### 8.3 Link Functions

By default, `bmm` uses appropriate link functions:

- **drift**: identity (can be positive or negative)
- **bound**: log (must be positive)
- **ndt**: log (must be positive)
- **zr**: logit (must be between 0 and 1)
- **sdrift, sndt**: log (must be positive)
- **szr**: logit (must be between 0 and 1)

You can customize links if needed:

``` r

# Use log link for drift (forces positive values)
model_custom_links <- ddm(
  rt = "rt",
  response = "response",
  links = list(drift = "log")
)
```

## 9 Comparison with Other Models

### 9.1 DDM vs Censored Shifted Wald (cswald)

**DDM advantages**:

- Models both correct and error responses
- Accounts for response bias and variability
- Richer theoretical interpretation
- Standard in the field

**cswald advantages**:

- Much faster to fit (minutes vs hours)
- Simpler interpretation with fewer parameters
- Adequate for high-accuracy tasks
- More stable with small samples

**Recommendation**: Use cswald for initial exploration or when errors
are rare (\<10%); use DDM for complete analysis when computational
resources allow.

### 9.2 DDM vs EZ-DDM (ezdm)

The EZ-DDM provides closed-form parameter estimates from summary
statistics (mean RT, accuracy, RT variance). It’s extremely fast but:

- Only provides point estimates (no uncertainty quantification)
- Cannot model individual differences hierarchically
- Cannot include covariates
- Assumes specific parameter values (e.g., \\z = 0.5\\, \\s_z = 0\\)

**Recommendation**: Use ezdm for quick initial estimates or when you
need to process many datasets rapidly; use full DDM for formal
inference.

## 10 Practical Guidelines

### 10.1 Sample Size Recommendations

Based on simulation studies and practical experience:

- **Minimum**: 50 trials per condition (fixed `zr`)
- **Adequate**: 100-200 trials per condition (free `zr`)
- **Recommended**: 200+ trials per condition for stable estimates
- **Hierarchical models**: At least 5-10 subjects per group

### 10.2 Convergence Diagnostics

Always check:

1.  **Rhat \< 1.01** for all parameters (indicates chain convergence)
2.  **ESS \> 400** (effective sample size, indicates sufficient
    sampling)
3.  **No divergent transitions** (indicates numerical issues)
4.  **Trace plots** look like “hairy caterpillars”
5.  **Posterior predictive checks** show good fit

### 10.3 Common Issues and Solutions

**Problem**: Divergent transitions - **Solution**: Increase
`adapt_delta` (e.g., 0.95 or 0.99)

**Problem**: Slow mixing / low ESS - **Solution**: Increase iterations,
use more informative priors, check reparameterization

**Problem**: Boundary estimates at extremes - **Solution**: Check data
quality, consider more informative priors

**Problem**: Non-decision time longer than minimum RT - **Solution**:
Check for fast outliers, consider data trimming

**Problem**: Starting point estimates at 0 or 1 - **Solution**: May
indicate model misspecification or insufficient data

## 11 References

Henrich, Felix, Raphael Hartmann, Valentin Pratz, Andreas Voss, and Karl
Christoph Klauer. 2024. “The Seven-Parameter Diffusion Model: An
Implementation in Stan for Bayesian Analyses.” *Behavior Research
Methods* 56 (4): 3102–16. <https://doi.org/10.3758/s13428-023-02179-1>.

Ratcliff, Roger. 1978. “A Theory of Memory Retrieval.” *Psychological
Review* 85 (2): 59–108. <https://doi.org/10.1037/0033-295X.85.2.59>.

Ratcliff, Roger, and Gail McKoon. 2008. “The Diffusion Decision Model:
Theory and Data for Two-Choice Decision Tasks.” *Neural Computation* 20
(4): 873–922. <https://doi.org/10.1162/neco.2008.12-06-420>.
