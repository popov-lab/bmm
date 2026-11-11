# Simulate data from a bmm model

`simulate_data()` draws data from a model given the values of its
parameters and a design, in the layout
[`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md) expects for
that model: the design's columns with the response columns added under
the names the model was specified with. The result can be passed to
[`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md) as it is, which
is what a parameter recovery study needs.

## Usage

``` r
simulate_data(model, pars, design = NULL, scale = c("link", "native"), ...)

# Default S3 method
simulate_data(model, pars, design = NULL, scale = c("link", "native"), ...)
```

## Arguments

- model:

  A `bmmodel` object, such as `mixture2p(resp_error = "y")`.

- pars:

  A named list of parameter values, named as
  [`parameter_info()`](https://venpopov.com/bmm/dev/reference/parameter_info.md)
  lists the model's parameters, or a data frame of them with one row or
  one row per row of `design`, such as the attribute `simulate_pars` of
  an earlier result. Each value is a single number or a vector with one
  element per row of `design`, element `i` going to row `i`. Parameters
  the model fixes may be left out, and then take their fixed value; a
  value given for one of them replaces it. Of alternative parameters,
  such as `criterion` and `bias` of
  [`sdt_yn()`](https://venpopov.com/bmm/dev/reference/sdt_yn.md), give
  one.

- design:

  A data frame with one row per observation of the model's data: one
  trial for trial-level models, one cell for the models of response
  counts. It holds the columns the model reads (see the section Design)
  and any others, such as condition labels or a subject identifier,
  which are returned unchanged. `NULL` for the models that read no
  column, and then `n_trials` in `...` gives the number of observations.

- scale:

  The scale of the values in `pars`. `"link"` (the default) is the scale
  [`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md) estimates the
  parameters on and [`summary()`](https://rdrr.io/r/base/summary.html)
  reports, for example `kappa` on the log scale. `"native"` is the scale
  after the inverse of the model's link, the scale of
  [`native_parameters()`](https://venpopov.com/bmm/dev/reference/native_parameters.md):
  `kappa` as a precision, `thetat` of
  [`mixture2p()`](https://venpopov.com/bmm/dev/reference/mixture2p.md)
  as the probability of recall. The weights of a softmax group, such as
  `thetat` and `thetant` of
  [`mixture3p()`](https://venpopov.com/bmm/dev/reference/mixture3p.md),
  are given as the probabilities themselves. They are the weights of a
  display with non-targets; where the set size is 1, the non-target
  weight goes to the target and to guessing in proportion, as in the
  fit, so `thetat = 0.6` with `thetant = 0.2` gives a probability of
  recall of 0.75 there.

- ...:

  `n_trials`, the number of observations when `design` is `NULL`;
  `formula`, the
  [`bmmformula()`](https://venpopov.com/bmm/dev/reference/bmmformula.md)
  with the activation functions, for the custom version of
  [`m3()`](https://venpopov.com/bmm/dev/reference/m3.md), whose
  parameters exist only in that formula.

## Value

A data frame: `design` with the response columns added. It carries the
attribute `simulate_pars`, a data frame with the value of every
parameter (fixed ones included) for every row on the link scale, the
values a recovery study compares the estimates with, and
`simulate_call`, the model's class and the scale. For
[`ezdm()`](https://venpopov.com/bmm/dev/reference/ezdm.md), rows whose
simulated mean response time is not positive are dropped with a warning,
because [`bmm()`](https://venpopov.com/bmm/dev/reference/bmm.md) rejects
them; their row numbers are in the attribute `dropped_rows`. A mean
response time below the non-decision time is a valid draw and is kept.

## Details

The values are converted from the link scale with the inverse of each
parameter's link, by the code
[`native_parameters()`](https://venpopov.com/bmm/dev/reference/native_parameters.md)
uses, so that a softmax group and the set-size gate of the non-target
weights of
[`mixture3p()`](https://venpopov.com/bmm/dev/reference/mixture3p.md) (no
non-target weight where no non-target was shown) are applied as in the
fit. Some parameters are defined on a transformed scale behind an
identity link, such as `sdratio` of
[`sdt_rating()`](https://venpopov.com/bmm/dev/reference/sdt_rating.md),
the log of the ratio of the standard deviations. For those, `"native"`
is the value itself; the transformation to the argument of the random
number generator (here `exp(sdratio)`) is applied on both scales.

The data are drawn with the model's `r*()` function, such as
[`rmixture2p()`](https://venpopov.com/bmm/dev/reference/mixture2p_dist.md)
or
[`rsdt_rating()`](https://venpopov.com/bmm/dev/reference/sdt_rating_dist.md).
They are not independent of the model's density functions: a recovery
study that has to rule out an error shared by the two needs a generator
written separately.

## Design

The columns each model reads from `design`, under the names the model
was specified with:

|  |  |
|----|----|
| Model | Columns |
| [`sdm()`](https://venpopov.com/bmm/dev/reference/sdm.md), [`ddm()`](https://venpopov.com/bmm/dev/reference/ddm.md), [`cswald()`](https://venpopov.com/bmm/dev/reference/cswald.md) | none |
| [`mixture2p()`](https://venpopov.com/bmm/dev/reference/mixture2p.md) | `set_size` (if it names a column) for the versions `"slot"` and `"slot_averaging"`; none for `"simple"` |
| [`mixture3p()`](https://venpopov.com/bmm/dev/reference/mixture3p.md) | `set_size` (if it names a column), `nt_features`, of which a row's first `set_size - 1` are read |
| [`imm()`](https://venpopov.com/bmm/dev/reference/imm.md) | as [`mixture3p()`](https://venpopov.com/bmm/dev/reference/mixture3p.md), plus `nt_distances` for the versions `"full"` and `"bsc"` |
| [`m3()`](https://venpopov.com/bmm/dev/reference/m3.md) | `n_trials`, the number of trials of each row; `num_options` when it names columns; for the custom version, the columns its activation functions read besides the parameters, as numbers (code a factor as 0/1 indicators) |
| [`mpt()`](https://venpopov.com/bmm/dev/reference/mpt.md) | `n_trials` (none with `response`, one row per trial); `tree_id` and `covariates` when the model has them |
| [`ezdm()`](https://venpopov.com/bmm/dev/reference/ezdm.md) | `n_trials` |
| [`sdt_yn()`](https://venpopov.com/bmm/dev/reference/sdt_yn.md) | `stimulus`, `n_trials` |
| [`sdt_mafc()`](https://venpopov.com/bmm/dev/reference/sdt_mafc.md) | `n_trials`; `m` when it names a column |
| [`sdt_ranking()`](https://venpopov.com/bmm/dev/reference/sdt_ranking.md) | `n_trials`; `m` when it names a column |
| [`sdt_rating()`](https://venpopov.com/bmm/dev/reference/sdt_rating.md), [`sdt_cdp()`](https://venpopov.com/bmm/dev/reference/sdt_cdp.md) | `stimulus`, `n_trials` |
| [`psychometric_function()`](https://venpopov.com/bmm/dev/reference/psychometric_function.md) | `intensity`; `n_trials` when it names a column (else each row is one 0/1 trial); `guess` and `lapse` when they name columns |

The simplex groups of
[`mpt()`](https://venpopov.com/bmm/dev/reference/mpt.md) are given as
probabilities that sum to 1 on both scales.

## Examples

``` r
# circular responses: two-parameter mixture model, 100 trials
data <- bmm::simulate_data(mixture2p(resp_error = "y"),
  pars = list(kappa = log(8), thetat = 0.5), n_trials = 100
)
head(data)
#>             y
#> 1  0.55736215
#> 2 -0.50792284
#> 3  0.36035066
#> 4  0.27127972
#> 5 -0.06787571
#> 6 -0.35931991
attr(data, "simulate_pars")[1, ]
#>   mu    kappa thetat
#> 1  0 2.079442    0.5

# the same on the native scale: kappa as a precision, thetat as P(memory)
data <- bmm::simulate_data(mixture2p(resp_error = "y"),
  pars = list(kappa = 8, thetat = 0.62), scale = "native", n_trials = 100
)

# response times and decisions: one value of the drift rate per trial
design <- data.frame(cond = rep(c("easy", "hard"), each = 50))
data <- bmm::simulate_data(ddm(rt = "rt", response = "response"),
  pars = list(drift = ifelse(design$cond == "easy", 2, 0.5), bound = 0, ndt = -1.5),
  design = design
)
head(data)
#>   cond        rt response
#> 1 easy 0.3394879        1
#> 2 easy 0.2731111        1
#> 3 easy 0.3965125        1
#> 4 easy 0.4204438        1
#> 5 easy 0.2694648        0
#> 6 easy 0.4030829        1

# counts: yes/no recognition, 50 trials per stimulus type
design <- data.frame(stimulus = c(0, 1), n_trials = 50)
bmm::simulate_data(sdt_yn(response = "n_old", stimulus = "stimulus", n_trials = "n_trials"),
  pars = list(d = 1.5, criterion = 0.5), design = design
)
#>   stimulus n_trials n_old
#> 1        0       50    16
#> 2        1       50    39

# categorical responses: memory measurement model, 100 trials per row
model <- m3(
  resp_cats = c("corr", "other", "npl"), num_options = c(1, 4, 5),
  choice_rule = "softmax", version = "ss"
)
bmm::simulate_data(model, pars = list(c = 2, a = 0.5), design = data.frame(n_trials = c(100, 100)))
#>   n_trials corr other npl
#> 1      100   48    29  23
#> 2      100   55    27  18
```
