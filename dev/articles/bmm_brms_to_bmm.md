# From brms to bmm

## 1 Who this article is for

You know `brms`. You have written `brm(y ~ x + (1 | id))`, read its
summary, and perhaps written a non-linear formula or two. Then you open
a `bmm` article and the formula has no response on its left-hand side,
the model is an object of its own, and the estimates come with links you
did not choose. This article translates between the two packages.

We fit one model twice: once in plain `brms`, once with `bmm`. The model
is the two-parameter mixture model for continuous reproduction tasks
(Zhang and Luck 2008), which `brms` can express with its built-in
mixture family. Both fits have the same likelihood and the same priors,
so they give the same estimates up to the noise of the sampler. With
that established, the rest of the article explains each difference in
the code: where the response went, why every parameter gets its own
formula, what `0 + setsize` and `1 + setsize` each give you, which scale
the estimates are on, and why `bmm` sets informative default priors. A
later section does the same comparison for signal detection:
[trial-level `sdt_yn()` against a probit generalized linear mixed
model](#signal-detection-as-a-probit-glmm).

## 2 The same model, written twice

### 2.1 The data

`zhang_luck_2008` ships with `bmm`. Eight participants saw one, two,
three or six colored squares, and after a delay they reproduced the
color of one of them on a color wheel. The response error is the
distance between the reported and the true color, in radians:

``` r

library(bmm)
library(brms)

dat <- zhang_luck_2008
head(dat[, c("subID", "setsize", "response_error")])
```

``` fansi
#> # A tibble: 6 × 3
#>   subID setsize response_error
#>   <dbl> <fct>            <dbl>
#> 1     1 6              -1.85  
#> 2     1 1              -0.0698
#> 3     1 2              -0.0349
#> 4     1 2               0.0349
#> 5     1 6              -1.54  
#> 6     1 3               0.733
```

`setsize` is a factor. With a numeric column, `0 + setsize` below would
fit a single slope instead of one value per set size.

The mixture model assumes that a response comes either from memory or
from a guess. Memory responses scatter around the true color with a von
Mises distribution, the normal distribution of the circle, whose
precision is `kappa`. Guesses fall anywhere on the circle with equal
probability. The second parameter is the probability that a response
comes from memory. We let both parameters differ between set sizes and
between participants.

### 2.2 The model in brms

In `brms`, this is a mixture of two von Mises distributions. Each
component has a location (`mu1`, `mu2`) and a precision (`kappa1`,
`kappa2`), and `theta1` is the weight of the first component. Turning
this into the mixture model takes three constraints. Both locations are
fixed at 0, because the response error is already relative to the
target. The precision of the second component is fixed at a value so
small that the distribution is flat, which makes it the guessing
distribution. Constant priors do the fixing:

``` r

brms_fit <- brm(
  formula = bf(
    response_error ~ 1,
    kappa1 ~ 0 + setsize + (1 | subID),
    theta1 ~ 0 + setsize + (1 | subID),
    mu2 ~ 1,
    kappa2 ~ 1
  ),
  family = mixture(von_mises, von_mises, order = "none"),
  prior = prior(constant(0), class = Intercept, dpar = mu1) +
    prior(constant(0), class = Intercept, dpar = mu2) +
    prior(constant(-100), class = Intercept, dpar = kappa2) +
    prior(normal(2, 1), class = b, dpar = kappa1) +
    prior(logistic(0, 1), class = b, dpar = theta1) +
    prior(exponential(1), class = sd, dpar = kappa1) +
    prior(exponential(1), class = sd, dpar = theta1),
  data = dat,
  cores = 4,
  control = list(adapt_delta = 0.95),
  seed = 681,
  backend = "cmdstanr",
  file = "fits/brmsfit_brms_to_bmm"
)
```

`kappa2` is fixed at -100 because `brms` estimates `kappa` on the log
scale, and `exp(-100)` is as close to zero as anyone needs. At the
default `adapt_delta` of 0.8, the sampler occasionally reports a
divergent transition for this model, so we raise it; this is a `brms`
argument and works the same in
[`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md). `file =` saves
the fit and loads it on later runs instead of refitting, so change the
name, or drop the argument, when you change the data or the model. The
last four priors are not the `brms` defaults. They are the priors `bmm`
uses for this model, and we pass them here so that both fits are the
same model. We come back to them in [the section on
priors](#why-the-default-priors-are-informative).

### 2.3 The model in bmm

The `bmm` version needs a model object, one formula per parameter, and
the data:

``` r

bmm_fit <- bmm(
  formula = bmf(
    kappa ~ 0 + setsize + (1 | subID),
    thetat ~ 0 + setsize + (1 | subID)
  ),
  data = dat,
  model = mixture2p(resp_error = "response_error"),
  cores = 4,
  control = list(adapt_delta = 0.95),
  seed = 681,
  backend = "cmdstanr",
  file = "fits/bmmfit_brms_to_bmm"
)
```

The constraints are gone from the code.
[`mixture2p()`](https://venpopov.com/bmm/dev/reference/mixture2p.md)
knows that the locations are 0, that the second component is the
guessing distribution, and which priors to use. What is left is what you
want to know: how `kappa` and `thetat`, the weight of the memory
responses (the `t` stands for target), depend on set size.

### 2.4 Same model, same estimates

`bmm` translates its formula into a `brmsformula` and passes it to
`brms`, which writes the Stan model. You can look at the formula `bmm`
generated:

``` r

bmm_fit$formula
#> response_error ~ 1 
#> mu ~ 1
#> kappa ~ 0 + setsize + (1 | subID)
#> thetat ~ 0 + setsize + (1 | subID)
```

`kappa` and `thetat` are parameters of a family that `bmm` defines
itself, and `mu`, the location of the memory responses, is fixed at 0.
`brms` does not know this family, so `bmm` adds the Stan code of its
likelihood: a mixture of a von Mises distribution centered on the target
and a uniform distribution for the guesses. The `brms` call above writes
the same likelihood with two von Mises components, the second one flat.
The two Stan programs differ, but the likelihood and the priors are the
same, so the posteriors agree up to the noise of the sampler. Both
columns are on the link scale, log `kappa` and logit `thetat`, which we
explain [below](#which-scale-the-estimates-are-on):

``` r

brms_est <- fixef(brms_fit)
bmm_est <- fixef(bmm_fit)
ci <- function(est, row) {
  sprintf(
    "%.2f [%.2f, %.2f]",
    est[row, "Estimate"], est[row, "Q2.5"], est[row, "Q97.5"]
  )
}
set_sizes <- levels(dat$setsize)
comparison <- data.frame(
  parameter = rep(c("kappa", "thetat"), each = length(set_sizes)),
  set_size = rep(set_sizes, 2),
  brms = c(
    sapply(paste0("kappa1_setsize", set_sizes), ci, est = brms_est),
    sapply(paste0("theta1_setsize", set_sizes), ci, est = brms_est)
  ),
  bmm = c(
    sapply(paste0("kappa_setsize", set_sizes), ci, est = bmm_est),
    sapply(paste0("thetat_setsize", set_sizes), ci, est = bmm_est)
  )
)
knitr::kable(comparison, row.names = FALSE)
```

| parameter | set_size | brms                   | bmm                    |
|:----------|:---------|:-----------------------|:-----------------------|
| kappa     | 1        | 2.81 \[2.69, 2.93\]    | 2.81 \[2.69, 2.93\]    |
| kappa     | 2        | 2.22 \[2.09, 2.35\]    | 2.22 \[2.09, 2.35\]    |
| kappa     | 3        | 1.99 \[1.83, 2.14\]    | 1.99 \[1.83, 2.14\]    |
| kappa     | 6        | 1.87 \[1.58, 2.15\]    | 1.87 \[1.58, 2.17\]    |
| thetat    | 1        | 5.03 \[4.06, 6.18\]    | 5.04 \[4.09, 6.26\]    |
| thetat    | 2        | 2.77 \[2.31, 3.24\]    | 2.77 \[2.27, 3.23\]    |
| thetat    | 3        | 1.53 \[1.13, 1.91\]    | 1.53 \[1.11, 1.91\]    |
| thetat    | 6        | -0.60 \[-1.03, -0.24\] | -0.60 \[-1.03, -0.26\] |

Why does `bmm` bring its own family when the `brms` mixture family fits
this model? Its family also covers versions of the model that a `brms`
mixture cannot write: the capacity versions (`version = "slot"` and
`"slot_averaging"`), in which the probability of a memory response
follows from a capacity `K` and the set size, and variable precision
(`variable_precision = TRUE`), in which the precision of memory varies
from trial to trial.
[`?mixture2p`](https://venpopov.com/bmm/dev/reference/mixture2p.md)
describes them. Besides the likelihood, `bmm` writes the `brms` code,
checks the data, and picks priors and starting values that suit the
model, as it does for models that `brms` has no family for at all, such
as [`sdm()`](https://venpopov.com/bmm/dev/reference/sdm.md) or
[`cswald()`](https://venpopov.com/bmm/dev/reference/cswald.md).

## 3 Where the response lives

In `brms`, the response is on the left-hand side of the main formula:
`response_error ~ 1`. In `bmm`, it is an argument of the model object:
`mixture2p(resp_error = "response_error")`. The formulas only ever have
parameters on their left-hand side.

The reason is that a measurement model usually needs more than one
column, and each column plays a different role. The three-parameter
mixture model needs the response error, the colors of the items that
were not cued, and the set size. A yes/no signal detection model fit to
counts needs the number of “old” responses, the number of trials, and
whether the item was old or new; with one row per trial, `n_trials`
drops out. In `brms`, these roles end up spread over the response,
addition terms such as `| trials()`, and predictors inside non-linear
formulas. In `bmm`, the model object names each column once, by its
role:

``` r

mixture3p(
  resp_error = "response_error", nt_features = "col_lure",
  set_size = "setsize", regex = TRUE
)
sdt_yn(response = "n_old", stimulus = "stimulus", n_trials = "n_trials")
```

Knowing the role of each column lets `bmm` check the data before
anything is compiled. For instance,
[`mixture2p()`](https://venpopov.com/bmm/dev/reference/mixture2p.md)
warns if the response error looks like degrees rather than radians. It
also lets `bmm` build the columns that the Stan model needs from the
ones you have, such as which of the non-target columns are filled at
each set size. The help page of each model, such as
[`?mixture2p`](https://venpopov.com/bmm/dev/reference/mixture2p.md),
lists the columns it expects.

## 4 Signal detection as a probit GLMM

A signal detection model can be fit without calling it that. With one
row per trial, `resp` coded 1 for an “old” response and `stimulus` coded
1 for a studied item and 0 for a new one, a yes/no recognition
experiment can be analyzed as a generalized linear mixed model (GLMM)
with a probit link:

``` r

brm(
  resp ~ 0 + Intercept + stimulus + (1 + stimulus | id) + (1 + stimulus | item),
  family = bernoulli("probit")
)
```

This is the equal-variance Gaussian signal detection model (DeCarlo
1998). In `bmm`, the same model is
[`sdt_yn()`](https://venpopov.com/bmm/dev/reference/sdt_yn.md) with one
row per trial: leave out `n_trials`, and `response` holds the 0 or 1 of
each trial.

### 4.1 How the parameters map

The probability of an “old” response is `pnorm(Intercept)` for a new
item and `pnorm(Intercept + stimulus)` for a studied one, where
`stimulus` stands for its coefficient. In signal detection terms, these
are `pnorm(-criterion)` and `pnorm(d - criterion)`: the evidence of a
new item is normal with mean 0 and standard deviation (SD) 1, that of a
studied item normal with mean `d`, and the participant says “old” when
the evidence lies above the criterion. The criterion is thus measured
from the mean of the new items. Matching the two pairs of probabilities
gives the intercept as minus the criterion and the `stimulus`
coefficient as `d`. So each part of the `brms` model has a `bmm` name:

| `bmm` ([`sdt_yn()`](https://venpopov.com/bmm/dev/reference/sdt_yn.md)) | `brms` probit GLMM |
|----|----|
| `d` | coefficient of `stimulus` |
| `criterion` | minus the intercept |
| `(1 \| id)` in the formula of `d` | `stimulus` slope by `id` |
| `(1 \| id)` in the formula of `criterion` | intercept by `id` |
| `(1 \| item)` in the formula of `d` | `stimulus` slope by `item` |
| `(1 \| item)` in the formula of `criterion` | intercept by `item` |
| `sdratio` (unequal variance) | none: the probit GLMM fixes both standard deviations at 1 |
| `bias` (criterion measured from the midpoint between the two distributions) | minus the intercept if `stimulus` is coded -0.5 / +0.5; with the 0 / 1 coding used here, `-(Intercept + stimulus / 2)` |
| `dist` (the distribution of the evidence) | partly: `bernoulli("logit")` is `dist = "logistic"` and `bernoulli("cloglog")` is `dist = "gumbel_max"`; no link gives `dist = "gumbel_min"` for a response coded 1 = “old” (see below the table) |

In the `dist` row, the link is that of the probability of an “old”
response, whereas
[`?sdt_yn`](https://venpopov.com/bmm/dev/reference/sdt_yn.md) gives each
`dist` by the cumulative distribution function of the evidence. The
function `1 - exp(-exp(x))`, the inverse of the complementary log-log
link, is listed there under `"gumbel_min"`; the probability of “old”
under `"gumbel_max"` is that same function of `d - criterion`. With the
response coded 1 for “new” instead, `bernoulli("cloglog")` gives
`dist = "gumbel_min"`: the probability of “new” is that function of
`criterion - d * stimulus`.

`stimulus` itself appears in no `bmm` formula: it is named in the model
object, `sdt_yn(stimulus = "stimulus")`, and `d` is its effect, so a
predictor in the formula of `d` is always one that interacts with
`stimulus`. A predictor of the `brms` intercept is a predictor of
`criterion`, with its sign flipped, and its interaction with `stimulus`
is a predictor of `d`. For a design factor `cond`,
`resp ~ 0 + Intercept + stimulus * cond` is
`bmf(d ~ 1 + cond, criterion ~ 1 + cond)`. The `d` coefficients are the
`stimulus` and `stimulus:cond` coefficients, and the `criterion`
coefficients are the intercept and `cond` coefficients, each multiplied
by -1.

The intercept by participant changes sign along with the intercept, so
its SD is the SD of the criterion. Beyond the parameters, `bmm` adds
helpers for the receiver operating characteristic (ROC) and the area
under the ROC curve (AUC), which the [signal detection
article](https://venpopov.com/bmm/articles/bmm_sdt.html) describes.

### 4.2 Both fits on one data set

To check the mapping, we simulate a recognition experiment in which 30
participants each see the same 60 items. Half of the participants
studied items 1 to 30, the other half items 31 to 60, so every item is
old for 15 participants and new for the other 15. `d` and `criterion`
vary between participants and between items. The data come from the
latent normal distributions of the model directly, not from a `bmm`
function, so that the check does not depend on `bmm`’s own code:

``` r

set.seed(4646)
truth <- c(
  d = 1.3, criterion = 0.6,
  sd_id_d = 0.4, sd_id_criterion = 0.25,
  sd_item_d = 0.3, sd_item_criterion = 0.2
)
n_id <- 30
n_item <- 60
rec <- expand.grid(id = seq_len(n_id), item = seq_len(n_item))
rec$stimulus <- as.integer((rec$item <= n_item / 2) == (rec$id %% 2 == 1))

id_d <- rnorm(n_id, 0, truth[["sd_id_d"]])
id_criterion <- rnorm(n_id, 0, truth[["sd_id_criterion"]])
item_d <- rnorm(n_item, 0, truth[["sd_item_d"]])
item_criterion <- rnorm(n_item, 0, truth[["sd_item_criterion"]])
trial_d <- truth[["d"]] + id_d[rec$id] + item_d[rec$item]
trial_criterion <- truth[["criterion"]] + id_criterion[rec$id] +
  item_criterion[rec$item]

# new items draw their evidence from normal(0, 1), old items from normal(d, 1)
evidence <- rnorm(nrow(rec), mean = rec$stimulus * trial_d)
rec$resp <- as.integer(evidence > trial_criterion)
rec$id <- factor(rec$id)
rec$item <- factor(rec$item)
head(rec)
#>   id item stimulus resp
#> 1  1    1        1    0
#> 2  2    1        0    1
#> 3  3    1        1    0
#> 4  4    1        0    0
#> 5  5    1        1    1
#> 6  6    1        0    1
```

The `brms` fit is the GLMM from above. As in the first example, we give
it the priors that `bmm` uses for
[`sdt_yn()`](https://venpopov.com/bmm/dev/reference/sdt_yn.md), so that
the priors on the means and the SDs agree.
[`?sdt_yn`](https://venpopov.com/bmm/dev/reference/sdt_yn.md) lists
them: normal(1, 1) on `d`, normal(0.5, 1.5) on `criterion` and
exponential(1) on every SD. On the intercept, the criterion’s prior has
its sign flipped:

``` r

glmm_fit <- brm(
  resp ~ 0 + Intercept + stimulus + (1 + stimulus | id) + (1 + stimulus | item),
  family = bernoulli("probit"),
  prior = prior(normal(-0.5, 1.5), class = b, coef = Intercept) +
    prior(normal(1, 1), class = b, coef = stimulus) +
    prior(exponential(1), class = sd),
  data = rec,
  chains = 4,
  cores = 4,
  seed = 646,
  backend = "cmdstanr",
  file = "fits/brmsfit_sdt_glmm_brms_to_bmm"
)
```

The `bmm` fit names the two columns in
[`sdt_yn()`](https://venpopov.com/bmm/dev/reference/sdt_yn.md) and gives
`d` and `criterion` a formula each:

``` r

trial_fit <- bmm(
  formula = bmf(
    d ~ 1 + (1 | id) + (1 | item),
    criterion ~ 1 + (1 | id) + (1 | item)
  ),
  data = rec,
  model = sdt_yn(response = "resp", stimulus = "stimulus"),
  chains = 4,
  cores = 4,
  seed = 646,
  backend = "cmdstanr",
  file = "fits/bmmfit_sdt_trial_brms_to_bmm"
)
```

`stimulus` appears in neither formula.
[`sdt_yn()`](https://venpopov.com/bmm/dev/reference/sdt_yn.md) reads it
from the column named in the model object and uses it to decide whether
`d` enters a trial, the job that the `stimulus` predictor does in
`brms`. The `1` in `d ~ 1` is the population `d`, the coefficient of
`stimulus` in the GLMM, and `(1 | id)` lets it vary between
participants, as the `stimulus` slope by `id` does.

The fits still differ in the correlations. `(1 + stimulus | id)` in
`brms` also estimates the correlation between each participant’s
intercept and slope, and the same for items. The `bmm` formula above
leaves these correlations out. Write `(1 + stimulus || id)` to drop them
in `brms`, or give the participant terms of `d` and `criterion` a shared
ID, `(1 | p | id)`, to add them in `bmm`. The `bmm` correlation is that
of `d` and `criterion`, which is minus the `brms` correlation of slope
and intercept, because the intercept is minus the criterion. The default
priors on a correlation also differ: `brms` uses LKJ(1) and `bmm` uses
LKJ(2):
[`make_stancode()`](https://paulbuerkner.com/brms/reference/stancode.html)
of the `brms` model shows LKJ(1), and
[`stancode()`](https://venpopov.com/bmm/dev/reference/stancode.bmmformula.md)
of a `bmm` model with `(1 | p | id)` shows LKJ(2).
`prior(lkj(2), class = cor)` gives the `brms` model the `bmm` prior.

The table puts the estimates of the two fits next to the values we
simulated from, each as the posterior mean and its 95% credible interval
(CrI). The `brms` intercept and its draws are multiplied by -1 to give
the criterion:

``` r

glmm_draws <- as_draws_df(glmm_fit)
trial_draws <- as_draws_df(trial_fit)
glmm_par <- list(
  d = glmm_draws$b_stimulus,
  criterion = -glmm_draws$b_Intercept,
  sd_id_d = glmm_draws$sd_id__stimulus,
  sd_id_criterion = glmm_draws$sd_id__Intercept,
  sd_item_d = glmm_draws$sd_item__stimulus,
  sd_item_criterion = glmm_draws$sd_item__Intercept
)
trial_par <- list(
  d = trial_draws$b_d_Intercept,
  criterion = trial_draws$b_criterion_Intercept,
  sd_id_d = trial_draws$sd_id__d_Intercept,
  sd_id_criterion = trial_draws$sd_id__criterion_Intercept,
  sd_item_d = trial_draws$sd_item__d_Intercept,
  sd_item_criterion = trial_draws$sd_item__criterion_Intercept
)
mean_cri <- function(x) {
  q <- quantile(x, c(0.025, 0.975))
  sprintf("%.2f [%.2f, %.2f]", mean(x), q[1], q[2])
}
knitr::kable(data.frame(
  parameter = c("d", "criterion", "SD of d by id", "SD of criterion by id",
                "SD of d by item", "SD of criterion by item"),
  `brms term` = c("stimulus", "-Intercept", "sd(stimulus) by id",
                  "sd(Intercept) by id", "sd(stimulus) by item",
                  "sd(Intercept) by item"),
  simulated = truth[names(glmm_par)],
  brms = vapply(glmm_par, mean_cri, character(1)),
  bmm = vapply(trial_par, mean_cri, character(1)),
  check.names = FALSE
), row.names = FALSE)
```

| parameter | brms term | simulated | brms | bmm |
|:---|:---|---:|:---|:---|
| d | stimulus | 1.30 | 1.24 \[1.04, 1.45\] | 1.24 \[1.05, 1.45\] |
| criterion | -Intercept | 0.60 | 0.59 \[0.45, 0.74\] | 0.59 \[0.46, 0.72\] |
| SD of d by id | sd(stimulus) by id | 0.40 | 0.37 \[0.16, 0.59\] | 0.33 \[0.16, 0.51\] |
| SD of criterion by id | sd(Intercept) by id | 0.25 | 0.25 \[0.09, 0.41\] | 0.23 \[0.10, 0.36\] |
| SD of d by item | sd(stimulus) by item | 0.30 | 0.27 \[0.03, 0.49\] | 0.25 \[0.04, 0.44\] |
| SD of criterion by item | sd(Intercept) by item | 0.20 | 0.18 \[0.03, 0.34\] | 0.16 \[0.02, 0.28\] |

For every row, each fit’s estimate lies inside the other fit’s 95% CrI,
and the two estimates differ by at most 0.05. Every simulated value lies
inside both intervals. The intervals of the SDs are wide, because 30
participants and 60 items are few groups from which to estimate a
spread.

### 4.3 The item effects on `d`

The item SD of `d` rests on fewer trials than the item SD of the
criterion. On a trial with a new item, `d` does not enter the
probability of an “old” response, `pnorm(-criterion)`, so an item’s
effect on `d` is estimated only from the trials on which the item was
old: 15 per item here, against 30 for its effect on the criterion. The
`brms` model says the same thing in its own terms: the `stimulus` slope
by item multiplies a predictor that is 0 on every trial with a new item.
In a design in which some items are never old, their effects on `d` stay
at the prior. With `bias` instead of `criterion`, `d` also moves the
criterion and enters every trial.

### 4.4 Aggregated counts

Item terms need trial rows, because counts summed over a participant’s
items no longer say which item each response came from. The same holds
for predictors that change from trial to trial. Without item terms or
trial-level predictors, you do not need trial rows.
[`sdt_yn()`](https://venpopov.com/bmm/dev/reference/sdt_yn.md) also
takes one row per participant and stimulus type, with the number of
“old” responses and the number of trials. This layout is faster, because
the time a fit takes grows with the number of rows: here, 60 rows
instead of 1800.

``` r

rec_counts <- aggregate(
  cbind(n_old = resp, n_trials = 1) ~ id + stimulus,
  data = rec, FUN = sum
)
agg_fit <- bmm(
  formula = bmf(d ~ 1 + (1 | id), criterion ~ 1 + (1 | id)),
  data = rec_counts,
  model = sdt_yn(response = "n_old", stimulus = "stimulus",
                 n_trials = "n_trials"),
  chains = 4,
  cores = 4,
  seed = 646,
  backend = "cmdstanr",
  file = "fits/bmmfit_sdt_agg_brms_to_bmm"
)
```

For this model, the 1800 trials and the 60 counts give the same
posterior, because the counts keep everything the trials say about the
parameters. To check this, we take every posterior draw of `agg_fit` and
compute two log-likelihoods: that of the 1800 trial rows, each written
as a count of 0 or 1 out of 1 trial, and that of the 60 counts. Their
difference, the gap, is printed as its range over all draws. If the gap
is the same in every draw, the two log-likelihoods are the same function
of the parameters up to a constant, and with the same priors they give
the same posterior. The second line computes the constant we expect:

``` r

rec_single <- transform(rec, n_old = resp, n_trials = 1)
loglik_gap <- rowSums(log_lik(agg_fit, newdata = rec_single)) -
  rowSums(log_lik(agg_fit))
range(loglik_gap)
#> [1] -889.9357 -889.9357
-sum(lchoose(rec_counts$n_trials, rec_counts$n_old))
#> [1] -889.9357
```

The gap is the same in every draw, up to rounding, and equals minus the
sum of the log binomial coefficients of the counts. That term does not
depend on any parameter, so it drops out of the posterior.

### 4.5 Running and comparing trial-level fits

Trial-level fits take longer because they evaluate the likelihood of
every trial. With `backend = "cmdstanr"`,
`bmm(..., threads = threading(4))` evaluates it on four threads in each
chain, which shortens the fit if your machine has cores to spare beyond
the chains.
[`threading()`](https://paulbuerkner.com/brms/reference/threading.html)
is the `brms` function that sets this up, and `threads = 4` is short for
it.

[`loo()`](https://mc-stan.org/loo/reference/loo.html) on a trial-level
fit leaves out one trial at a time; on the aggregated fit, it leaves out
one count, that is, all trials of a participant and stimulus type at
once. These are different predictions, so compare the
[`loo()`](https://mc-stan.org/loo/reference/loo.html) results only
between fits to the same rows. To compare a trial-level fit with an
aggregated one, use `kfold(fit, group = "id")`, which leaves out whole
participants in either layout.

## 5 One formula per parameter

Because the response lives in the model object, a `bmmformula` is a list
of regressions, one for each parameter of the model. Each one is a
regular `brms` formula, with population-level effects, group-level
effects in the `lme4` syntax, and non-linear terms if you need them. The
[bmmformula
article](https://venpopov.com/bmm/articles/bmm_bmmformula.html) covers
the syntax in full, including how to fix a parameter to a constant.

### 5.1 `0 + setsize` or `1 + setsize`

Most `bmm` articles write `0 + setsize`, and you may be used to
`1 + setsize` (or just `setsize`, which is the same thing). Both
describe the same regression of each parameter on set size. They differ
in what each coefficient means, and, as we show below, in which default
prior each coefficient gets.

With `0 + setsize`, the formula has no intercept, and there is one
coefficient per set size. Each coefficient is the parameter in that cell
of the design: `kappa_setsize1` is `kappa` at set size 1,
`kappa_setsize6` is `kappa` at set size 6, both on the log scale
explained [below](#which-scale-the-estimates-are-on). This is often
easier when you start out, because the summary shows the estimate for
each condition directly, and those are the numbers you want to plot or
report. If you are not sure which coding to use, start with this one.

With `1 + setsize`, the intercept is the parameter in the reference
condition, set size 1, and the other coefficients are differences from
it. `kappa_setsize6` is now how much log `kappa` at set size 6 differs
from log `kappa` at set size 1. This coding is just as suitable, and it
reads naturally when the difference between conditions is your question.
To get the parameter in a cell, you add the coefficients:

``` r

bmm_fit_ref <- bmm(
  formula = bmf(
    kappa ~ 1 + setsize + (1 | subID),
    thetat ~ 1 + setsize + (1 | subID)
  ),
  data = dat,
  model = mixture2p(resp_error = "response_error"),
  cores = 4,
  control = list(adapt_delta = 0.95),
  seed = 681,
  backend = "cmdstanr",
  file = "fits/bmmfit_brms_to_bmm_reference"
)
```

``` r

fixef(bmm_fit_ref)[c("kappa_Intercept", "kappa_setsize6"), ]
#>                  Estimate  Est.Error      Q2.5      Q97.5
#> kappa_Intercept  2.811732 0.06189117  2.689923  2.9343583
#> kappa_setsize6  -0.926667 0.15239452 -1.219430 -0.6323925

draws_ref <- as_draws_df(bmm_fit_ref)
draws_cell <- as_draws_df(bmm_fit)
kappa6_ref <- draws_ref$b_kappa_Intercept + draws_ref$b_kappa_setsize6
kappa6_cell <- draws_cell$b_kappa_setsize6
rbind(
  `1 + setsize` = quantile(kappa6_ref, c(0.025, 0.5, 0.975)),
  `0 + setsize` = quantile(kappa6_cell, c(0.025, 0.5, 0.975))
)
#>                 2.5%      50%    97.5%
#> 1 + setsize 1.593512 1.883484 2.175504
#> 0 + setsize 1.583832 1.875491 2.166188
```

Add the coefficients draw by draw, as above, and summarize afterwards.
The draws carry the uncertainty of each coefficient and how the
coefficients covary, which adding the two point estimates loses. The
reverse also works: a difference between set sizes is the difference of
two cell coefficients of the `0 + setsize` fit, again taken per draw.

Doing this by hand for every cell and parameter gets tedious.
[`native_parameters()`](https://venpopov.com/bmm/dev/reference/native_parameters.md)
does it for you. It evaluates every parameter in every cell of the
design, whatever the coding of the formula, and reports it on the scale
of the model (more on scales in the next section). Here are both fits
side by side:

``` r

np_cell <- native_parameters(bmm_fit, re_formula = NA, summary = TRUE)
np_ref <- native_parameters(bmm_fit_ref, re_formula = NA, summary = TRUE)
fmt <- function(np) sprintf("%.2f [%.2f, %.2f]", np$Estimate, np$Q2.5, np$Q97.5)
np_cell <- np_cell[np_cell$parameter %in% c("kappa", "thetat"), ]
np_ref <- np_ref[np_ref$parameter %in% c("kappa", "thetat"), ]
knitr::kable(data.frame(
  parameter = np_cell$parameter,
  set_size = np_cell$setsize,
  `0 + setsize` = fmt(np_cell),
  `1 + setsize` = fmt(np_ref[match(paste(np_cell$parameter, np_cell$setsize),
                                   paste(np_ref$parameter, np_ref$setsize)), ]),
  check.names = FALSE
), row.names = FALSE)
```

| parameter | set_size | 0 + setsize            | 1 + setsize            |
|:----------|:---------|:-----------------------|:-----------------------|
| kappa     | 1        | 16.64 \[14.67, 18.70\] | 16.67 \[14.73, 18.81\] |
| kappa     | 2        | 9.27 \[8.07, 10.53\]   | 9.27 \[8.12, 10.56\]   |
| kappa     | 3        | 7.31 \[6.24, 8.52\]    | 7.31 \[6.22, 8.49\]    |
| kappa     | 6        | 6.58 \[4.87, 8.72\]    | 6.66 \[4.92, 8.81\]    |
| thetat    | 1        | 0.99 \[0.98, 1.00\]    | 0.99 \[0.98, 1.00\]    |
| thetat    | 2        | 0.94 \[0.91, 0.96\]    | 0.94 \[0.91, 0.96\]    |
| thetat    | 3        | 0.82 \[0.75, 0.87\]    | 0.83 \[0.77, 0.88\]    |
| thetat    | 6        | 0.36 \[0.26, 0.44\]    | 0.36 \[0.28, 0.45\]    |

Unlike the coefficients in the first table, these values are on the
scale of the model: `kappa` is the precision itself, not its log, and
`thetat` is shown as the probability of a memory response, not its
logit. `re_formula = NA` sets the participant effects to zero, which
gives the parameters of a typical participant;
[`?native_parameters`](https://venpopov.com/bmm/dev/reference/native_parameters.md)
explains how this differs from the population average under a
non-identity link.

On the scale of the model, the two codings give practically the same
estimates. In every cell, the estimates of the two fits differ by less
than 2.5% of the `0 + setsize` estimate; the largest gap is thetat at
set size 6, 0.355 against 0.364. Pick the coding whose coefficients you
want to read.

What the coding does change is the default prior of each coefficient.
`bmm` sets its default priors per coefficient, and the value of a
parameter in one cell is a different quantity from the difference
between two cells: the first lies somewhere in the plausible range of
the parameter, the second is centered on zero. With `0 + setsize`, every
coefficient is a cell, so every one gets the prior for the parameter
itself, normal(2, 1) for `kappa` on the log scale and logistic(0, 1) for
`thetat`. With `1 + setsize`, the differences get the prior for effects,
normal(0, 1) for `kappa` and normal(0, 2.5) for `thetat`. The prior for
the parameter itself, normal(2, 1), does not go on the reference cell,
though. Like `brms`, `bmm` centers the predictors before fitting, so
that prior applies to the intercept of the centered predictors, which in
this balanced design is the average of log `kappa` over the four set
sizes. The intercept in the summary is still the reference cell: that
average minus a quarter of each difference. The prior of every cell
therefore comes from both priors. For `kappa` it is a normal with
standard deviation 1.09 at set size 1 and 1.30 at the other set sizes,
wider than the standard deviation of 1 that every cell gets with
`0 + setsize`. With this data set, the table above shows that the
difference barely matters. The table hides one exception, visible only
on the logit scale: `thetat` at set size 1. Near a probability of 1, the
data say little about the logit, and with `1 + setsize` the effects
prior, centered on zero, pulls the reference cell towards the other set
sizes. On the logit scale, the two fits estimate `thetat` at set size 1
as 5.04 and 4.85, and the set-size differences of `thetat` differ by up
to 0.25 logits. With fewer trials per cell, the priors carry more
weight, and the two codings can give different estimates. So check
[`default_prior()`](https://venpopov.com/bmm/dev/reference/default_prior.bmmformula.md)
for the coding you use, as shown in [the section on
priors](#why-the-default-priors-are-informative), and set your own
priors if the defaults do not match what you know.

The participant effects can be coded either way, too. `(1 | subID)`, as
here, gives every participant one shift for all set sizes.
`(0 + setsize | subID)` and `(1 + setsize | subID)` both let
participants differ per set size, coded as cell values or as differences
from the reference cell.

## 6 Which scale the estimates are on

The summary of a `bmm` fit lists the link of each parameter:

``` r

summary(bmm_fit)
```

``` fansi
  Model: mixture2p(resp_error = "response_error") 
  Links: mu = tan_half; kappa = log; thetat = logit 
Formula: mu = 0
         kappa ~ 0 + setsize + (1 | subID)
         thetat ~ 0 + setsize + (1 | subID) 
   Data: dat (Number of observations: 4000)
  Draws: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
         total post-warmup draws = 4000

Multilevel Hyperparameters:
~subID (Number of levels: 8) 
                     Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
sd(kappa_Intercept)      0.09      0.06     0.01     0.24 1.01     1009     1591
sd(thetat_Intercept)     0.40      0.17     0.17     0.82 1.00     1333     1804

Regression Coefficients:
                Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
kappa_setsize1      2.81      0.06     2.69     2.93 1.00     1940     1440
kappa_setsize2      2.22      0.07     2.09     2.35 1.00     2071     1819
kappa_setsize3      1.99      0.08     1.83     2.14 1.00     2359     2101
kappa_setsize6      1.87      0.15     1.58     2.17 1.00     3294     2287
thetat_setsize1     5.04      0.54     4.09     6.26 1.00     3003     2142
thetat_setsize2     2.77      0.24     2.27     3.23 1.00     1709     1758
thetat_setsize3     1.53      0.20     1.11     1.91 1.00     1378     1218
thetat_setsize6    -0.60      0.19    -1.03    -0.26 1.00     1308     1371

Constant Parameters:
                 Value
mu_Intercept      0.00

Draws were sampled using sample(hmc). For each parameter, Bulk_ESS
and Tail_ESS are effective sample size measures, and Rhat is the potential
scale reduction factor on split chains (at convergence, Rhat = 1).
```

`kappa` must be positive, and `thetat` turns into a probability. Neither
can be predicted by a linear model directly, so `bmm` estimates them on
an unbounded scale and transforms them back. `kappa` has a log link, so
the coefficients are log precisions. `thetat` has a logit link: the
probability that a response comes from memory is `plogis(thetat)`, so a
`thetat` of 0 means that memory and guessing are equally likely. Every
coefficient in the summary, every prior, and every value you fix a
parameter to lives on this link scale.

To report a parameter on its own scale, transform the draws and
summarize afterwards, not the other way around. The mean of
[`exp()`](https://rdrr.io/r/base/Log.html) of the draws is not
[`exp()`](https://rdrr.io/r/base/Log.html) of their mean. For set size
6:

``` r

c(
  `exp(mean)` = exp(mean(kappa6_cell)),
  `mean(exp)` = mean(exp(kappa6_cell))
)
#> exp(mean) mean(exp) 
#>  6.508365  6.578416
```

[`native_parameters()`](https://venpopov.com/bmm/dev/reference/native_parameters.md),
used above, follows this order for every parameter. For the mixture
models, [`k2sd()`](https://venpopov.com/bmm/dev/reference/k2sd.md)
converts `kappa` into the circular standard deviation in radians, which
is the measure of precision many papers report:

``` r

kappa_draws <- native_parameters(bmm_fit, re_formula = NA, pars = "kappa")
sd_deg <- k2sd(kappa_draws$value) * 180 / pi
tapply(sd_deg, kappa_draws$setsize, median)
#>        1        2        3        6 
#> 14.27779 19.38580 22.03948 23.42801
```

A coefficient on its own transforms to something else. With
`1 + setsize`, `exp(b_kappa_setsize6)` is the factor by which `kappa` at
set size 6 is larger or smaller than at set size 1, not `kappa` at set
size 6. The help page of each model lists the links of its parameters.

## 7 Why the default priors are informative

`brms` uses flat priors on regression coefficients unless you set them.
For the `brms` version of this model, the defaults are:

``` r

default_prior(
  bf(
    response_error ~ 1,
    kappa1 ~ 0 + setsize + (1 | subID),
    theta1 ~ 0 + setsize + (1 | subID),
    mu2 ~ 1,
    kappa2 ~ 1
  ),
  data = dat,
  family = mixture(von_mises, von_mises, order = "none")
) |>
  as.data.frame() |>
  transform(prior = ifelse(prior == "", "(flat)", prior)) |>
  subset(dpar %in% c("kappa1", "theta1") & coef == "" & group == "",
         select = c(dpar, class, prior)) |>
  knitr::kable(row.names = FALSE)
```

| dpar   | class | prior                |
|:-------|:------|:---------------------|
| kappa1 | b     | (flat)               |
| kappa1 | sd    | student_t(3, 0, 2.5) |
| theta1 | b     | (flat)               |
| theta1 | sd    | student_t(3, 0, 2.5) |

`class = b` covers every coefficient of `kappa1` and `theta1`, so each
of them has a flat prior. That suits a regression whose coefficients
could take any value. It does not suit a measurement model. Its
parameters stand for cognitive processes, and we know roughly which
values are plausible before seeing any data. `bmm` puts this knowledge
into its default priors:

``` r

default_prior(
  bmf(kappa ~ 0 + setsize + (1 | subID), thetat ~ 0 + setsize + (1 | subID)),
  data = dat,
  model = mixture2p(resp_error = "response_error")
) |>
  as.data.frame() |>
  subset(dpar %in% c("kappa", "thetat") & coef == "" & group == "",
         select = c(dpar, class, prior)) |>
  knitr::kable(row.names = FALSE)
```

| dpar   | class | prior          |
|:-------|:------|:---------------|
| kappa  | sd    | exponential(1) |
| kappa  | b     | normal(2, 1)   |
| thetat | sd    | exponential(1) |
| thetat | b     | logistic(0, 1) |

The normal(2, 1) on the log of `kappa` places 95% of its mass between a
`kappa` of 1.0 and 52.5, that is, between a circular standard deviation
of 71 and 8 degrees. Values in the thousands, which would mean a
circular standard deviation of 1.8 degrees or less, get next to no prior
mass. The logistic(0, 1) on `thetat` is flat on the probability scale:
before seeing the data, every probability of a memory response between 0
and 1 is equally likely. It still does work that the `brms` default does
not: a flat prior on the logit is improper and, on the probability
scale, piles its mass at 0 and 1.

There is a second reason, which has to do with the mixture. When `kappa`
comes close to 0, the memory distribution becomes as flat as the
guessing distribution, and the data cannot tell the two apart. In that
region, `thetat` can take any value without changing the fit. A flat
prior does nothing to keep the sampler out of that region; a prior that
holds `kappa` in a plausible range does.

What this means for you:

- The more trials per participant and condition, the less the estimates
  depend on the priors. With few trials, the priors pull the estimates
  towards their center.
- Look at the priors before you fit, with
  [`default_prior()`](https://venpopov.com/bmm/dev/reference/default_prior.bmmformula.md)
  as above, and at the help page of the model, which lists them.
- If you know more than the defaults, or something different, set your
  own priors with `prior =` in
  [`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md). They are
  `brms` priors on the `bmm` parameter names, with `nlpar` or `dpar` as
  [`default_prior()`](https://venpopov.com/bmm/dev/reference/default_prior.bmmformula.md)
  shows it, for example
  `prior(normal(1.5, 0.5), class = b, dpar = kappa)`, again on the log
  scale. With `0 + setsize`, this sets the prior of every cell of
  `kappa`. With `1 + setsize`, it sets the prior of the differences
  only. The intercept takes `class = Intercept`, as
  [`default_prior()`](https://venpopov.com/bmm/dev/reference/default_prior.bmmformula.md)
  shows for that coding, and that prior applies to the centered
  intercept, not to the reference cell. The [article on extracting model
  information](https://venpopov.com/bmm/articles/bmm_extract_info.html)
  shows the details.
- Bayes factors depend on the priors much more than estimates do. If you
  plan to compute them, choose the priors on purpose rather than taking
  the defaults.

## 8 What stays brms

A `bmm` fit is a `brmsfit` with extra information, so the `brms` toolbox
works on it:
[`fixef()`](https://rdrr.io/pkg/nlme/man/fixed.effects.html),
[`ranef()`](https://rdrr.io/pkg/nlme/man/random.effects.html),
[`as_draws_df()`](https://mc-stan.org/posterior/reference/draws_df.html),
[`pp_check()`](https://venpopov.com/bmm/dev/reference/pp_check.bmmfit.md),
[`conditional_effects()`](https://venpopov.com/bmm/dev/reference/conditional_effects.bmmfit.md),
[`hypothesis()`](https://paulbuerkner.com/brms/reference/hypothesis.brmsfit.html),
[`loo()`](https://mc-stan.org/loo/reference/loo.html) and the rest. As
shown above, `bmm_fit$formula` holds the `brmsformula` that `bmm`
generated, and `stancode(bmm_fit)` the Stan code.
`summary(bmm_fit, backend = "brms")` prints the `brms` summary instead
of the `bmm` one. It shows the formula as `brms` received it, with
`response_error ~ 1` and `mu ~ 1`, and it lists the fixed `mu` as a
coefficient named `Intercept`, with the value 0 and no diagnostics. If
you want to change a model beyond what `bmm` offers, these give you the
`brms` code to start from.

## 9 Where to go next

- [Get started](https://venpopov.com/bmm/articles/bmm.html) lists the
  models by task and walks through a first fit.
- The [bmmformula
  article](https://venpopov.com/bmm/articles/bmm_bmmformula.html) covers
  the formula syntax.
- The [mixture models
  article](https://venpopov.com/bmm/articles/bmm_mixture_models.html)
  fits the two- and three-parameter mixture models to a larger data set.
- The [signal detection
  article](https://venpopov.com/bmm/articles/bmm_sdt.html) covers the
  yes/no, rating, m-AFC and ranking models, with counts or trial rows.
  The section [Signal detection as a probit
  GLMM](#signal-detection-as-a-probit-glmm) above is the short version
  for [`sdt_yn()`](https://venpopov.com/bmm/dev/reference/sdt_yn.md).
- [Extracting model
  information](https://venpopov.com/bmm/articles/bmm_extract_info.html)
  shows how to get the default priors, the Stan code and the Stan data
  before you fit.

## References

DeCarlo, Lawrence T. 1998. “Signal Detection Theory and Generalized
Linear Models.” *Psychological Methods* 3 (2): 186–205.
<https://doi.org/10.1037/1082-989X.3.2.186>.

Zhang, Weiwei, and Steven J. Luck. 2008. “Discrete Fixed-Resolution
Representations in Visual Working Memory.” *Nature* 453 (7192): 233–35.
<https://doi.org/10.1038/nature06860>.
