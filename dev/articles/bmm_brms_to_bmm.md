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
mixture family. Both fits run the same Stan model, so they give the same
estimates. With that established, the rest of the article explains each
difference in the code: where the response went, why every parameter
gets its own formula, what `0 + setsize` and `1 + setsize` each give
you, which scale the estimates are on, and why `bmm` sets informative
default priors.

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
  backend = "cmdstanr",
  file = "assets/brmsfit_brms_to_bmm"
)
```

`kappa2` is fixed at -100 because `brms` estimates `kappa` on the log
scale, and `exp(-100)` is as close to zero as anyone needs. At the
default `adapt_delta` of 0.8, the sampler occasionally reports a
divergent transition for this model, so we raise it; this is a `brms`
argument and works the same in
[`bmm()`](https://popov-lab.github.io/bmm/dev/reference/bmm.md).
`file =` saves the fit and loads it on later runs instead of refitting,
so change the name, or drop the argument, when you change the data or
the model. The last four priors are not the `brms` defaults. They are
the priors `bmm` uses for this model, and we pass them here so that both
fits are the same model. We come back to them in [the section on
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
  backend = "cmdstanr",
  file = "assets/bmmfit_brms_to_bmm"
)
```

The constraints are gone from the code.
[`mixture2p()`](https://popov-lab.github.io/bmm/dev/reference/mixture2p.md)
knows that the locations are 0, that the second component is the
guessing distribution, and which priors to use. What is left is what you
want to know: how `kappa` and `thetat`, the weight of the memory
responses (the `t` stands for target), depend on set size.

### 2.4 Same model, same estimates

`bmm` translates its formula into a `brmsformula`, and `brms` writes the
Stan model from that. You can look at the formula `bmm` generated:

``` r

bmm_fit$formula
#> response_error ~ 1 
#> mu1 ~ 1
#> kappa ~ 0 + setsize + (1 | subID)
#> thetat ~ 0 + setsize + (1 | subID)
#> kappa2 ~ 1
#> mu2 ~ 1
#> kappa1 ~ kappa
#> theta1 ~ thetat
```

`kappa` and `thetat` are non-linear parameters that `bmm` passes on to
`kappa1` and `theta1`. `mu1`, `mu2` and `kappa2` are fixed with constant
priors, as in the `brms` call. Apart from the parameter names, the Stan
code of the two fits is the same program with the same priors, so the
posteriors agree up to the noise of the sampler. Both columns are on the
link scale, log `kappa` and logit `thetat`, which we explain
[below](#which-scale-the-estimates-are-on):

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
| kappa     | 1        | 2.81 \[2.69, 2.93\]    | 2.81 \[2.69, 2.92\]    |
| kappa     | 2        | 2.22 \[2.10, 2.35\]    | 2.23 \[2.10, 2.35\]    |
| kappa     | 3        | 1.98 \[1.82, 2.14\]    | 1.99 \[1.83, 2.14\]    |
| kappa     | 6        | 1.87 \[1.58, 2.16\]    | 1.87 \[1.57, 2.15\]    |
| thetat    | 1        | 5.04 \[4.13, 6.26\]    | 5.06 \[4.08, 6.32\]    |
| thetat    | 2        | 2.77 \[2.28, 3.23\]    | 2.77 \[2.32, 3.22\]    |
| thetat    | 3        | 1.53 \[1.11, 1.92\]    | 1.54 \[1.14, 1.90\]    |
| thetat    | 6        | -0.60 \[-1.02, -0.24\] | -0.59 \[-0.99, -0.27\] |

For this model, `bmm` adds nothing to the Stan program. It writes the
`brms` code, checks the data, and picks priors and starting values that
suit the model. For models that `brms` has no family for, such as
[`sdm()`](https://popov-lab.github.io/bmm/dev/reference/sdm.md) or
[`cswald()`](https://popov-lab.github.io/bmm/dev/reference/cswald.md),
`bmm` also supplies the Stan code of the likelihood.

## 3 Where the response lives

In `brms`, the response is on the left-hand side of the main formula:
`response_error ~ 1`. In `bmm`, it is an argument of the model object:
`mixture2p(resp_error = "response_error")`. The formulas only ever have
parameters on their left-hand side.

The reason is that a measurement model usually needs more than one
column, and each column plays a different role. The three-parameter
mixture model needs the response error, the colors of the items that
were not cued, and the set size. A yes/no signal detection model needs
the number of “old” responses, the number of trials, and whether the
item was old or new. In `brms`, these roles end up spread over the
response, addition terms such as `| trials()`, and predictors inside
non-linear formulas. In `bmm`, the model object names each column once,
by its role:

``` r

mixture3p(
  resp_error = "response_error", nt_features = "col_lure",
  set_size = "setsize", regex = TRUE
)
sdt_yn(response = "n_old", stimulus = "stimulus", n_trials = "n_trials")
```

Knowing the role of each column lets `bmm` check the data before
anything is compiled. For instance,
[`mixture2p()`](https://popov-lab.github.io/bmm/dev/reference/mixture2p.md)
warns if the response error looks like degrees rather than radians. It
also lets `bmm` build the columns that the Stan model needs from the
ones you have, such as which of the non-target columns are filled at
each set size. The help page of each model, such as
[`?mixture2p`](https://popov-lab.github.io/bmm/dev/reference/mixture2p.md),
lists the columns it expects.

## 4 One formula per parameter

Because the response lives in the model object, a `bmmformula` is a list
of regressions, one for each parameter of the model. Each one is a
regular `brms` formula, with population-level effects, group-level
effects in the `lme4` syntax, and non-linear terms if you need them. The
[bmmformula
article](https://popov-lab.github.io/bmm/articles/bmm_bmmformula.html)
covers the syntax in full, including how to fix a parameter to a
constant.

### 4.1 `0 + setsize` or `1 + setsize`

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
  backend = "cmdstanr",
  file = "assets/bmmfit_brms_to_bmm_reference"
)
```

``` r

fixef(bmm_fit_ref)[c("kappa_Intercept", "kappa_setsize6"), ]
#>                   Estimate  Est.Error      Q2.5      Q97.5
#> kappa_Intercept  2.8105353 0.06038921  2.694992  2.9271767
#> kappa_setsize6  -0.9254779 0.15041941 -1.222528 -0.6306745

draws_ref <- as_draws_df(bmm_fit_ref)
draws_cell <- as_draws_df(bmm_fit)
kappa6_ref <- draws_ref$b_kappa_Intercept + draws_ref$b_kappa_setsize6
kappa6_cell <- draws_cell$b_kappa_setsize6
rbind(
  `1 + setsize` = quantile(kappa6_ref, c(0.025, 0.5, 0.975)),
  `0 + setsize` = quantile(kappa6_cell, c(0.025, 0.5, 0.975))
)
#>                 2.5%      50%    97.5%
#> 1 + setsize 1.596731 1.883600 2.176051
#> 0 + setsize 1.570206 1.869862 2.153408
```

Add the coefficients draw by draw, as above, and summarize afterwards.
The draws carry the uncertainty of each coefficient and how the
coefficients covary, which adding the two point estimates loses. The
reverse also works: a difference between set sizes is the difference of
two cell coefficients of the `0 + setsize` fit, again taken per draw.

Doing this by hand for every cell and parameter gets tedious.
[`native_parameters()`](https://popov-lab.github.io/bmm/dev/reference/native_parameters.md)
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
| kappa     | 1        | 16.63 \[14.74, 18.60\] | 16.65 \[14.81, 18.67\] |
| kappa     | 2        | 9.28 \[8.15, 10.53\]   | 9.27 \[8.09, 10.51\]   |
| kappa     | 3        | 7.31 \[6.26, 8.49\]    | 7.31 \[6.26, 8.46\]    |
| kappa     | 6        | 6.54 \[4.81, 8.61\]    | 6.66 \[4.94, 8.81\]    |
| thetat    | 1        | 0.99 \[0.98, 1.00\]    | 0.99 \[0.98, 1.00\]    |
| thetat    | 2        | 0.94 \[0.91, 0.96\]    | 0.94 \[0.91, 0.96\]    |
| thetat    | 3        | 0.82 \[0.76, 0.87\]    | 0.83 \[0.76, 0.88\]    |
| thetat    | 6        | 0.36 \[0.27, 0.43\]    | 0.36 \[0.28, 0.45\]    |

Unlike the coefficients in the first table, these values are on the
scale of the model: `kappa` is the precision itself, not its log, and
`thetat` is shown as the probability of a memory response, not its
logit. `re_formula = NA` sets the participant effects to zero, which
gives the parameters of a typical participant;
[`?native_parameters`](https://popov-lab.github.io/bmm/dev/reference/native_parameters.md)
explains how this differs from the population average under a
non-identity link.

On the scale of the model, the two codings give practically the same
estimates. In every cell, the estimates of the two fits differ by less
than 1.8% of the `0 + setsize` estimate; the largest gap is kappa at set
size 6, 6.54 against 6.66. Pick the coding whose coefficients you want
to read.

What the coding does change is the default prior of each coefficient.
`bmm` sets its default priors per coefficient, and the value of a
parameter in one cell is a different quantity from the difference
between two cells: the first lies somewhere in the plausible range of
the parameter, the second is centered on zero. With `0 + setsize`, every
coefficient is a cell, so every one gets the prior for the parameter
itself, normal(2, 1) for `kappa` on the log scale and logistic(0, 1) for
`thetat`. With `1 + setsize`, only the intercept, the reference cell,
gets that prior. The differences get the prior for effects, normal(0, 1)
for `kappa` and normal(0, 2.5) for `thetat`. Any other cell is the
intercept plus a difference, so its prior is that of a sum of two
independent coefficients. For `kappa`, the two normal priors add up to
normal(2, 1.41), whose standard deviation is the square root of 1² + 1²:
wider than the normal(2, 1) that the same cell gets with `0 + setsize`.
With this data set, the table above shows that the difference barely
matters. The table hides one exception, visible only on the logit scale:
`thetat` at set size 1. Near a probability of 1, the data say little
about the logit, and with `1 + setsize` the effects prior, centered on
zero, pulls the reference cell towards the other set sizes. On the logit
scale, the two fits estimate `thetat` at set size 1 as 5.06 and 4.70,
and the set-size differences of `thetat` differ by up to 0.42 logits.
With fewer trials per cell, the priors carry more weight, and the two
codings can give different estimates. So check
[`default_prior()`](https://popov-lab.github.io/bmm/dev/reference/default_prior.bmmformula.md)
for the coding you use, as shown in [the section on
priors](#why-the-default-priors-are-informative), and set your own
priors if the defaults do not match what you know.

The participant effects can be coded either way, too. `(1 | subID)`, as
here, gives every participant one shift for all set sizes.
`(0 + setsize | subID)` and `(1 + setsize | subID)` both let
participants differ per set size, coded as cell values or as differences
from the reference cell.

## 5 Which scale the estimates are on

The summary of a `bmm` fit lists the link of each parameter:

``` r

summary(bmm_fit)
```

``` fansi
  Model: mixture2p(resp_error = "response_error") 
  Links: mu1 = tan_half; kappa = log; thetat = logit 
Formula: mu1 = 0
         kappa ~ 0 + setsize + (1 | subID)
         thetat ~ 0 + setsize + (1 | subID) 
   Data: dat (Number of observations: 4000)
  Draws: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
         total post-warmup draws = 4000

Multilevel Hyperparameters:
~subID (Number of levels: 8) 
                     Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
sd(kappa_Intercept)      0.09      0.06     0.01     0.23 1.01      675     1420
sd(thetat_Intercept)     0.39      0.16     0.17     0.77 1.00     1236     1794

Regression Coefficients:
                Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
kappa_setsize1      2.81      0.06     2.69     2.92 1.00     2722     2078
kappa_setsize2      2.23      0.07     2.10     2.35 1.00     2361     2233
kappa_setsize3      1.99      0.08     1.83     2.14 1.00     3120     2769
kappa_setsize6      1.87      0.15     1.57     2.15 1.00     2948     2364
thetat_setsize1     5.06      0.56     4.08     6.32 1.00     3247     2143
thetat_setsize2     2.77      0.23     2.32     3.22 1.00     2053     2104
thetat_setsize3     1.54      0.19     1.14     1.90 1.00     1648     1481
thetat_setsize6    -0.59      0.19    -0.99    -0.27 1.00     1602     1516

Constant Parameters:
                  Value
mu1_Intercept      0.00

Draws were sampled using sample(hmc). For each parameter, Bulk_ESS
and Tail_ESS are effective sample size measures, and Rhat is the potential
scale reduction factor on split chains (at convergence, Rhat = 1).
```

`kappa` must be positive, and `thetat` turns into a probability. Neither
can be predicted by a linear model directly, so `bmm` estimates them on
an unbounded scale and transforms them back. `kappa` has a log link, so
the coefficients are log precisions. `thetat` has a logit link. `brms`
fixes the weight of the second, guessing component (`theta2`) at 0 on
this scale, so a `thetat` of 0 means that memory and guessing are
equally likely, and the probability that a response comes from memory is
`plogis(thetat)`. Every coefficient in the summary, every prior, and
every value you fix a parameter to lives on this link scale.

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
#>  6.475063  6.543754
```

[`native_parameters()`](https://popov-lab.github.io/bmm/dev/reference/native_parameters.md),
used above, follows this order for every parameter. For the mixture
models,
[`k2sd()`](https://popov-lab.github.io/bmm/dev/reference/k2sd.md)
converts `kappa` into the circular standard deviation in radians, which
is the measure of precision many papers report:

``` r

kappa_draws <- native_parameters(bmm_fit, re_formula = NA, pars = "kappa")
sd_deg <- k2sd(kappa_draws$value) * 180 / pi
tapply(sd_deg, kappa_draws$setsize, median)
#>        1        2        3        6 
#> 14.27673 19.38314 22.02748 23.50076
```

A coefficient on its own transforms to something else. With
`1 + setsize`, `exp(b_kappa_setsize6)` is the factor by which `kappa` at
set size 6 is larger or smaller than at set size 1, not `kappa` at set
size 6. The help page of each model lists the links of its parameters.

## 6 Why the default priors are informative

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
  subset(nlpar %in% c("kappa", "thetat") & coef == "" & group == "",
         select = c(nlpar, class, prior)) |>
  knitr::kable(row.names = FALSE)
```

| nlpar  | class | prior          |
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
  [`default_prior()`](https://popov-lab.github.io/bmm/dev/reference/default_prior.bmmformula.md)
  as above, and at the help page of the model, which lists them.
- If you know more than the defaults, or something different, set your
  own priors with `prior =` in
  [`bmm()`](https://popov-lab.github.io/bmm/dev/reference/bmm.md). They
  are `brms` priors on the `bmm` parameter names, with `nlpar` or `dpar`
  as
  [`default_prior()`](https://popov-lab.github.io/bmm/dev/reference/default_prior.bmmformula.md)
  shows it, for example
  `prior(normal(1.5, 0.5), class = b, nlpar = kappa)`, again on the log
  scale. With `0 + setsize`, this sets the prior of every cell of
  `kappa`. With `1 + setsize`, it sets the prior of the differences
  only; the intercept takes `coef = "Intercept"`, as
  [`default_prior()`](https://popov-lab.github.io/bmm/dev/reference/default_prior.bmmformula.md)
  shows for that coding. The [article on extracting model
  information](https://popov-lab.github.io/bmm/articles/bmm_extract_info.html)
  shows the details.
- Bayes factors depend on the priors much more than estimates do. If you
  plan to compute them, choose the priors on purpose rather than taking
  the defaults.

## 7 What stays brms

A `bmm` fit is a `brmsfit` with extra information, so the `brms` toolbox
works on it:
[`fixef()`](https://rdrr.io/pkg/nlme/man/fixed.effects.html),
[`ranef()`](https://rdrr.io/pkg/nlme/man/random.effects.html),
[`as_draws_df()`](https://mc-stan.org/posterior/reference/draws_df.html),
[`pp_check()`](https://popov-lab.github.io/bmm/dev/reference/pp_check.bmmfit.md),
[`conditional_effects()`](https://popov-lab.github.io/bmm/dev/reference/conditional_effects.bmmfit.md),
[`hypothesis()`](https://paulbuerkner.com/brms/reference/hypothesis.brmsfit.html),
[`loo()`](https://mc-stan.org/loo/reference/loo.html) and the rest. As
shown above, `bmm_fit$formula` holds the `brmsformula` that `bmm`
generated, and `stancode(bmm_fit)` the Stan code.
`summary(bmm_fit, backend = "brms")` prints the `brms` summary instead
of the `bmm` one. It lists the link of `theta1` as `identity`, because
`brms` applies the softmax inside the mixture family. If you want to
change a model beyond what `bmm` offers, these give you the `brms` code
to start from.

## 8 Where to go next

- [Get started](https://popov-lab.github.io/bmm/articles/bmm.html) lists
  the models by task and walks through a first fit.
- The [bmmformula
  article](https://popov-lab.github.io/bmm/articles/bmm_bmmformula.html)
  covers the formula syntax.
- The [mixture models
  article](https://popov-lab.github.io/bmm/articles/bmm_mixture_models.html)
  fits the two- and three-parameter mixture models to a larger data set.
- [Extracting model
  information](https://popov-lab.github.io/bmm/articles/bmm_extract_info.html)
  shows how to get the default priors, the Stan code and the Stan data
  before you fit.

## References

Zhang, Weiwei, and Steven J. Luck. 2008. “Discrete Fixed-Resolution
Representations in Visual Working Memory.” *Nature* 453 (7192): 233–35.
<https://doi.org/10.1038/nature06860>.
