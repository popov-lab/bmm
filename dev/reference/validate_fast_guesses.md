# Test if fast contaminants show random guessing behavior

Uses Bayesian Beta-Binomial conjugate analysis to test whether fast
flagged contaminants show random guessing (~50% accuracy for 2AFC). The
test computes the posterior distribution for the proportion of "upper"
responses and uses a Savage-Dickey Bayes Factor to quantify evidence for
or against the guessing hypothesis.

## Usage

``` r
validate_fast_guesses(
  contam_flag,
  rt_data,
  response,
  threshold_type = c("quantile", "absolute"),
  rt_threshold = 0.25,
  prior_alpha = 1,
  prior_beta = 1,
  guess_prob = 0.5,
  credible_mass = 0.95
)
```

## Arguments

- contam_flag:

  Logical vector indicating which trials were flagged as contaminants

- rt_data:

  Numeric vector of reaction times (in seconds)

- response:

  Response data in any format accepted by `.convert_response_to_upper()`
  (numeric 0/1, logical, character, factor)

- threshold_type:

  Character. How to interpret `rt_threshold`:

  - `"quantile"` (default): Use rt_threshold as quantile (0-1)

  - `"absolute"`: Use rt_threshold as absolute RT in seconds

- rt_threshold:

  Numeric. Threshold for defining "fast" trials. Interpretation depends
  on `threshold_type`:

  - If `threshold_type = "quantile"` (default): Quantile of RT
    distribution (e.g., 0.25 = 25th percentile). Default 0.25.

  - If `threshold_type = "absolute"`: Absolute RT value in seconds
    (e.g., 0.25 = 250ms).

- prior_alpha, prior_beta:

  Numeric. Parameters for Beta prior distribution. Default 1,1 gives
  uniform prior. Values \> 1 express prior belief about response
  proportions.

- guess_prob:

  Numeric. Null hypothesis value for guessing probability. Default 0.5
  (equal probability of upper/lower responses).

- credible_mass:

  Numeric. Probability mass for Highest Density Interval. Default 0.95
  for 95% HDI. Common alternatives: 0.90, 0.99.

## Value

List with components:

- `method`: "bayesian"

- `prop_upper`: Observed proportion of upper responses

- `hdi_lower`, `hdi_upper`: 95% Highest Density Interval bounds

- `bf_01`: Bayes Factor for H0 (guessing) vs H1 (non-random)

- `guess_in_hdi`: Logical, whether `guess_prob` is in HDI

- `bf_evidence`: Character, evidence category on Jeffreys scale

- `posterior_alpha`, `posterior_beta`: Posterior Beta parameters

- `n_tested`: Number of fast flagged trials tested

- `rt_threshold`: Actual RT threshold value used (in seconds)

- `threshold_type`: Type of threshold used ("quantile" or "absolute")

- `credible_mass`: Credible mass used for HDI computation

- `mean_rt_tested`: Mean RT of tested trials

## Details

The function performs a Bayesian test using the Beta-Binomial conjugate
prior-posterior relationship. With a Beta(alpha, beta) prior and
observing n_upper "upper" responses out of n_tested trials, the
posterior is:

Beta(alpha + n_upper, beta + n_lower)

The Savage-Dickey Bayes Factor compares the posterior and prior
densities at the null hypothesis value (default 0.5):

BF_01 = posterior_density(guess_prob) / prior_density(guess_prob)

Evidence categories follow Jeffreys (1961) scale:

- BF \> 10: Strong evidence for guessing

- BF \> 3: Moderate evidence for guessing

- BF \> 1: Anecdotal evidence for guessing

- BF \< 1/3: Moderate evidence against guessing

- BF \< 1/10: Strong evidence against guessing

**Note**: This function can be used as a standalone validation step
after obtaining contamination probabilities from
[`flag_contaminant_rts()`](https://venpopov.com/bmm/dev/reference/flag_contaminant_rts.md).

## References

Jeffreys, H. (1961). Theory of Probability (3rd ed.). Oxford University
Press.

## See also

[`flag_contaminant_rts()`](https://venpopov.com/bmm/dev/reference/flag_contaminant_rts.md)
for obtaining contamination probabilities

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate data with random guessing on fast trials
set.seed(123)
n <- 100
rt <- c(runif(20, 0.15, 0.30), rgamma(80, 5, 10))
response <- c(rbinom(20, 1, 0.5), rbinom(80, 1, 0.7))
contam_flag <- rt < 0.35

# Test using quantile threshold (default, adaptive)
result1 <- validate_fast_guesses(
  contam_flag = contam_flag,
  rt_data = rt,
  response = response,
  threshold_type = "quantile",
  rt_threshold = 0.30 # 30th percentile
)

# Test using absolute threshold (fixed RT)
result2 <- validate_fast_guesses(
  contam_flag = contam_flag,
  rt_data = rt,
  response = response,
  threshold_type = "absolute",
  rt_threshold = 0.30 # 300ms
)

print(result1$bf_01) # Bayes Factor
print(result1$bf_evidence) # Evidence category
print(result1$guess_in_hdi) # Is 0.5 in 95% HDI?
print(result1$threshold_type) # "quantile"
} # }
```
