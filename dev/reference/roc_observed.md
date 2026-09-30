# Observed ROC points from a fitted SDT model

Computes empirical (observed) ROC points from the response-count data
used to fit a
[`sdt_rating()`](https://popov-lab.github.io/bmm/dev/reference/sdt_rating.md)
or [`sdt_yn()`](https://popov-lab.github.io/bmm/dev/reference/sdt_yn.md)
model. Rating models pool counts across observations within each
stimulus type (and optional condition) to produce cumulative
hit/false-alarm rates at each confidence threshold; binary models
produce one operating point per criterion level (the empirical
counterpart of the model-implied points from
[`roc_sdt()`](https://popov-lab.github.io/bmm/dev/reference/roc_sdt.md)).

## Usage

``` r
roc_observed(fit, conditions = NULL)
```

## Arguments

- fit:

  A `bmmfit` object from a rating or binary SDT model.

- conditions:

  Optional character vector of column names to condition on. For rating
  models, `NULL` pools all observations into one ROC. For binary models,
  `NULL` auto-detects the criterion-varying predictor(s).

## Value

A data frame of class `"bmm_sdt_roc_observed"` with columns `FA`, `Hit`,
and any condition columns. Rating models additionally include the (0,0)
and (1,1) endpoints. The attribute `model_type` is `"rating"` or
`"binary"`; rating results also carry `n_ratings`.

## See also

[`roc_sdt()`](https://popov-lab.github.io/bmm/dev/reference/roc_sdt.md),
[`plot.bmm_sdt_roc()`](https://popov-lab.github.io/bmm/dev/reference/plot.bmm_sdt_roc.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Three base-rate conditions shift the criterion, which identifies sdratio
dat <- expand.grid(id = 1:20, stimulus = c(0L, 1L),
                   condition = c("liberal", "neutral", "strict"))
dat$n_trials <- 100L
criteria <- c(liberal = -0.5, neutral = 0, strict = 0.5)
dat$n_old <- rsdt_yn(nrow(dat), dat$n_trials, dat$stimulus, d = 1.5,
                     criterion = criteria[as.character(dat$condition)],
                     sdratio = 1.3)

fit <- bmm(
  formula = bmf(d ~ 1, criterion ~ 0 + condition, sdratio ~ 1),
  data = dat,
  model = sdt_yn(response = "n_old", stimulus = "stimulus",
                 n_trials = "n_trials"),
  cores = 4,
  backend = "cmdstanr"
)

# One point per base-rate condition, from the response counts
obs <- roc_observed(fit)
obs
plot(roc_sdt(fit), observed = obs)
} # }
```
