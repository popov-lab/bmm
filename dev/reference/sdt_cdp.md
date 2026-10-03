# Continuous Dual-Process Signal Detection Theory Model

Continuous Dual-Process Signal Detection Theory Model

## Usage

``` r
sdt_cdp(
  response = "",
  stimulus,
  n_new,
  n_old,
  dist = "normal",
  threshold_type = c("parsimonious", "equidistant", "log_distance"),
  links = NULL,
  ...
)
```

## Arguments

- response:

  An optional common prefix for the response count columns (default
  `""`, i.e. the bare canonical names `new1`, `know2`, ... that
  [`aggregate_sdt_cdp_data()`](https://popov-lab.github.io/bmm/dev/reference/aggregate_sdt_cdp_data.md)
  and
  [`rsdt_cdp()`](https://popov-lab.github.io/bmm/dev/reference/sdt_cdp_dist.md)
  produce).

- stimulus:

  The name of the variable coding the stimulus type. Must be coded as 0
  (new/lure) and 1 (old/target).

- n_new, n_old:

  Integer numbers of "new" and "old" confidence levels. Together with
  `response` they determine the response count columns and the threshold
  parameters.

- dist:

  The noise distribution. Only `"normal"` is currently supported (the
  CDP model is inherently Gaussian).

- threshold_type:

  Character. Threshold parameterization on the strength axis:
  `"parsimonious"` (default) and `"equidistant"` use a single `spacing`
  parameter; `"log_distance"` (Paulewicz & Blaut, 2022) estimates the
  `n_new + n_old - 2` distances between adjacent thresholds freely, each
  as a `deltaN` parameter on the log scale (the log width of the
  interval between thresholds `N` and `N + 1`, as in
  [`sdt_rating()`](https://popov-lab.github.io/bmm/dev/reference/sdt_rating.md)).
  `criterion` stays on the old/new boundary, threshold `n_new`.

- links:

  A named list of link functions for the parameters, e.g.
  `links = list(drec = "log")`. Only `dfam`, `drec`, `criterion` and
  `rcrit` can be set, to `"log"`, `"softplus"`, `"logit"`, `"probit"` or
  `"identity"`, the links the model can invert inside its formula.
  `sigmar`, `rho`, `kcrit` and the threshold parameters keep their
  identity links, because their fixed values and the model's own
  transformations ([`exp()`](https://rdrr.io/r/base/Log.html),
  [`tanh()`](https://rdrr.io/r/base/Hyperbolic.html)) assume it. A
  positive-range link (`"log"`, `"softplus"`) on `criterion` or `rcrit`
  states that the criterion is positive: the sampler then starts next to
  0 on the native scale and cannot reach negative values.

- ...:

  used internally for testing, ignore it

## Value

An object of class `bmmodel`

## Details

- **Domain:** Recognition Memory

- **Task:** Old/New Recognition with Remember/Know Judgments

- **Name:** Continuous Dual-Process Signal Detection Theory (CDP)

- **Citation:**

  - Wixted, J. T., & Mickes, L. (2010). A continuous dual-process model
    of remember/know judgments. Psychological Review, 117(4), 1025-1054.
    https://doi.org/10.1037/a0020874

- **Requirements:**

  Provide aggregated data with one row per cell (unique combination of
  predictors and stimulus class) and one integer count column per
  response category:

  - `new1` ... `new<n_new>` for 'new' judgments

  - `know<k>` and `remember<k>` (and optionally `guess<k>`) for 'old'
    judgments, with k on the unified confidence scale n_new+1 ...
    n_new+n_old

  - an optional common column prefix is set via `response`

  - stimulus: 0 = new/lure, 1 = old/target

Use aggregate_sdt_cdp_data() to build these columns from long-format
(trial-level) data

- **Parameters:**

  - `dfam`: Familiarity sensitivity: target mean on the familiarity
    axis, in lure-SD units (d' on that axis, where targets and lures
    share SD 1)

  - `drec`: Recollection sensitivity: target mean on the recollection
    axis, in lure-SD units. The target SD is exp(sigmar), so once sigmar
    is estimated drec is not the d_a the other SDT models report as d

  - `criterion`: Response bias: old/new boundary on the strength (F+R)
    axis

  - `spacing`: Threshold spacing: controls distance between adjacent
    thresholds (exp(spacing) ensures positive spacing)

  - `rcrit`: Remember criterion: threshold on the recollection axis

  - `sigmar`: Log SD of the recollection target distribution, so
    exp(sigmar) is the SD and 0 means SD = 1

  - `rho`: Familiarity-recollection correlation on an unconstrained
    scale; tanh(rho) is the correlation and 0 means independent
    processes

  - `kcrit`: Know criterion: threshold on the familiarity axis that
    splits Know from Guess. Active only when the data include 'guess'
    counts

- **Fixed parameters:**

  - `sigmar` = 0

  - `rho` = 0

  - `kcrit` = -100

- **Default parameter links:**

  - dfam = identity; drec = identity; criterion = identity; spacing =
    identity; rcrit = identity; sigmar = identity; rho = identity; kcrit
    = identity

- **Default priors:**

  - `dfam`:

    - `main`: normal(1, 1)

    - `effects`: normal(0, 0.5)

    - `sd`: exponential(1)

  - `drec`:

    - `main`: normal(1, 1)

    - `effects`: normal(0, 0.5)

    - `sd`: exponential(1)

  - `criterion`:

    - `main`: normal(0, 1.5)

    - `effects`: normal(0, 0.5)

    - `sd`: exponential(2)

  - `spacing`:

    - `main`: normal(0, 0.5)

    - `effects`: normal(0, 0.3)

    - `sd`: exponential(2)

  - `rcrit`:

    - `main`: normal(0, 1)

    - `effects`: normal(0, 0.5)

    - `sd`: exponential(2)

  - `sigmar`:

    - `main`: normal(0, 0.5)

    - `effects`: normal(0, 0.3)

    - `sd`: exponential(2)

  - `rho`:

    - `main`: normal(0, 0.5)

    - `effects`: normal(0, 0.3)

    - `sd`: exponential(2)

  - `kcrit`:

    - `main`: normal(0, 1)

    - `effects`: normal(0, 0.5)

    - `sd`: exponential(2)

The continuous dual-process model (Wixted & Mickes, 2010) assumes two
continuous memory signals per item, Familiarity (F) and Recollection
(R), with correlation `tanh(rho)`. Old/new confidence is read off the
aggregate strength S = F + R; the Remember/Know judgment splits "old"
responses on R (against `rcrit`); an optional Know/Guess split uses F
(against `kcrit`).

**Response format.** The model is fit to aggregated counts, like
[`sdt_rating()`](https://popov-lab.github.io/bmm/dev/reference/sdt_rating.md):
one row per cell (e.g. participant x stimulus class x condition) and one
integer count column per response category. The columns follow a fixed
naming scheme – `new1 ... new<n_new>` for "new" judgments, and `know<k>`
/ `remember<k>` (optionally `guess<k>`) for "old" judgments, where `k`
runs over the unified confidence scale `n_new + 1 ... n_new + n_old` (1
= most confident "new", K = most confident "old") – so the constructor
only needs the numbers of confidence levels, not a long vector of column
names. An optional common prefix is set via `response` (e.g.
`response = "cdp_"` for columns `cdp_new1`, `cdp_know2`, ...). Use
[`aggregate_sdt_cdp_data()`](https://popov-lab.github.io/bmm/dev/reference/aggregate_sdt_cdp_data.md)
to build these columns from long-format (trial-level or count) data;
[`rsdt_cdp()`](https://popov-lab.github.io/bmm/dev/reference/sdt_cdp_dist.md)
generates them directly. The numbers of "new" and "old" confidence
levels need not be equal (e.g. the 1-new / 5-old scale of Rotello et
al., 2005). Whether the Know/Guess split is modelled is **driven by the
data**: include `guess<k>` columns to fit the three-way
Remember/Know/Guess model.

**Variants via the formula.** By default `sigmar`, `rho`, and `kcrit`
are fixed (equal recollection variance, independent F and R, no
Know/Guess split), recovering the classic independent-process CDP. Free
them through the formula:

- `sigmar ~ 1` estimates unequal recollection variance.

- `rho ~ 1` estimates the within-item F-R correlation; `rho` also
  accepts predictors (e.g. `rho ~ condition`), bounded to (-1, 1) via an
  internal `tanh`. This structural, within-item correlation is distinct
  from a between-subject correlation of random effects, which is
  available separately through brms syntax such as `(1 |p| id)` on any
  parameter.

- `kcrit ~ 1` estimates the Know/Guess criterion. The data must then
  contain the `guess<k>` count columns, and data with those columns must
  free `kcrit`: with `kcrit` at its fixed value every Guess response has
  probability 0, so
  [`bmm()`](https://popov-lab.github.io/bmm/dev/reference/bmm.md)
  refuses that combination. The fixed value -100 is a sentinel that
  switches the split off;
  [`dsdt_cdp()`](https://popov-lab.github.io/bmm/dev/reference/sdt_cdp_dist.md)
  and
  [`rsdt_cdp()`](https://popov-lab.github.io/bmm/dev/reference/sdt_cdp_dist.md)
  use `kcrit = NULL` for the same purpose.

In the three-way model a Guess response in confidence bin k is possible
only while `rcrit + kcrit` exceeds the lower threshold of that bin,
because Guess needs R below `rcrit` and F below `kcrit`, so S = F + R
below their sum. Below that boundary the Guess category has probability
0, and when `rcrit` (or `kcrit`) varies across participants the
posterior explores it: hierarchical Remember/Know/Guess fits can report
divergent transitions near this boundary without biased estimates. If
they do, raise `adapt_delta` (0.95 or 0.99), give `rcrit` a tighter
group-level prior, or check which participants' Guess counts in the top
bins drive it. The Remember/Know model has no such boundary.

The kernel's gradient is wrong, with the value right, on two
measure-zero sets that `init = 0` starts on (a threshold exactly on the
strength mean, and `rcrit == kcrit`); bmm's own initial values avoid
them, and a fit started on them walks off in its first step.

## Sensitivity scales

`dfam` and `drec` are the target means \\\mu_F\\ and \\\mu_R\\ of Wixted
and Mickes (2010), each in units of the lure SD on its own axis, so
estimates can be compared with the values they report. On the
familiarity axis targets and lures both have SD 1, and `dfam` is \\d'\\
there. On the recollection axis the target SD is `exp(sigmar)`, so
`drec` is standardized by the lure SD alone. With `sigmar` fixed at 0 it
is \\d'\\ as well; with `sigmar` estimated it is not the balanced index
\\d_a\\ that
[`sdt_yn()`](https://popov-lab.github.io/bmm/dev/reference/sdt_yn.md)
and
[`sdt_rating()`](https://popov-lab.github.io/bmm/dev/reference/sdt_rating.md)
report as `d`.

This matters when `drec` is compared across conditions whose `sigmar`
differs: a difference in the recollection SD alone can then appear as a
credible difference in `drec`, the problem the sensitivity section of
[`sdt_yn()`](https://popov-lab.github.io/bmm/dev/reference/sdt_yn.md)
describes for \\d'\\. To compare such conditions, compute the balanced
recollection sensitivity from the posterior draws as
`drec / sqrt((1 + exp(2 * sigmar)) / 2)`, or give `sigmar` the same
formula in every condition being compared.

Neither `dfam` nor `drec` is the discriminability of the old/new
decision, which is read off the aggregate strength \\F + R\\ and depends
on both sensitivities, on `sigmar` and on `rho`.

When no Remember/Know split is available (confidence ratings only), use
[`sdt_rating()`](https://popov-lab.github.io/bmm/dev/reference/sdt_rating.md)
instead.

## References

Wixted, J. T., & Mickes, L. (2010). A continuous dual-process model of
remember/know judgments. *Psychological Review*, *117*(4), 1025–1054.
[doi:10.1037/a0020874](https://doi.org/10.1037/a0020874)

Paulewicz, B., & Blaut, A. (2022). The general causal cumulative model
of ordinal response. *PsyArXiv*.
[doi:10.31234/osf.io/e7a3x](https://doi.org/10.31234/osf.io/e7a3x)

Rotello, C. M., Macmillan, N. A., Reeder, J. A., & Wong, M. (2005). The
remember response: Subject to bias, graded, and not a process-pure
indicator of recollection. *Psychonomic Bulletin & Review*, *12*(5),
865–873. [doi:10.3758/BF03196778](https://doi.org/10.3758/BF03196778)

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate a Remember/Know data set (3 new + 3 old confidence levels) for
# 20 subjects: rsdt_cdp() returns the count columns sdt_cdp() expects
dat <- expand.grid(id = 1:20, stimulus = c(0L, 1L))
thresholds <- c(-1.1, -0.5, 0, 0.6, 1.3)
dat <- cbind(dat, rsdt_cdp(nrow(dat), 200, dat$stimulus,
                           dfam = 0.8, drec = 1.0,
                           thresholds = thresholds, rcrit = 0.7,
                           n_new = 3))

model <- sdt_cdp(stimulus = "stimulus", n_new = 3, n_old = 3)

fit <- bmm(
  formula = bmf(dfam ~ 1, drec ~ 1, criterion ~ 1,
                spacing ~ 1, rcrit ~ 1),
  data = dat, model = model, backend = "cmdstanr"
)

# Estimate unequal recollection variance and the F-R correlation
fit_uv <- bmm(
  formula = bmf(dfam ~ 1, drec ~ 1, criterion ~ 1,
                spacing ~ 1, rcrit ~ 1, sigmar ~ 1, rho ~ 1),
  data = dat, model = model, backend = "cmdstanr"
)
} # }
```
