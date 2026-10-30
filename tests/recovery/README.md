# Recovery plans

One file per model: `tests/recovery/<model>.R`. It builds a validation plan with
`bmmtools::validation_plan()` and returns it as its last expression. The plan
says what is simulated, how it is fitted and what counts as recovered; it does
not fit anything when it is sourced.

```r
# tests/recovery/mixture2p.R
bmmtools::validation_plan(
  "mixture2p",
  generator = function(pars, n_trials, model) {
    rvm_reject <- function(n, kappa) {
      out <- numeric(0)
      while (length(out) < n) {
        x <- runif(n, -pi, pi)
        out <- c(out, x[runif(n) < exp(kappa * (cos(x) - 1))])
      }
      out[seq_len(n)]
    }
    memory <- runif(n_trials) < pars$thetat
    y <- ifelse(
      memory,
      rvm_reject(n_trials, pars$kappa) + pars$mu1,
      runif(n_trials, -pi, pi)
    )
    y <- (y + pi) %% (2 * pi) - pi
    stats::setNames(data.frame(y), model$resp_vars$resp_error)
  },
  seed = 2026
)
```

The `generator` takes `(pars, n_trials, model)` and returns one subject's rows
with the model's column names. Without it, a default entry simulates with
bmm's own `r*` function, which does not meet the independent-generator
requirement, so a density bug could cancel out.

The defaults entries `mixture2p` and `ezdm` simulate with bmm's own `r*`
function. Models that have a bmmtools generator adapter but no defaults entry
(for example `mixture3p`) fall back to it silently, so always pass
`generator =`.

A model without a bmmtools defaults entry (0.3.0 has `mixture2p` and `ezdm`)
cannot be built from its name. Pass the model object and `pars`, the
population values on the link scale, as in
`bmmtools::validation_plan(bmm::mixture3p(...), pars = c(kappa = 2, thetat = 1, thetant = 0), generator = ..., sds = ..., grid = ..., seed = 2026)`.
The generator receives `pars` on the natural scale instead (`kappa`, `thetat`
as a probability, `mu1` filled with its fixed value). Add `sds =` and `grid =`
because the hierarchical design needs subject SDs and a subjects-by-trials
grid. The plan must cover a
true-against-recovered comparison, a hierarchical design with several subjects
whose parameters vary, and data from an independent generator, so that a bug in
the density cannot cancel out (see `.github/CONTRIBUTING.md`).

To run a plan and render its card (`local/` is untracked):

```r
plan <- source("tests/recovery/mixture2p.R")$value
result <- bmmtools::run_validation(plan, dir = "local/recovery/mixture2p")
bmmtools::render_card(result, output = "local/recovery/mixture2p/card.html")
```

The run takes hours of sampling on a laptop, and you run it. Until a card site
exists, attach the rendered card.html to the PR (see
`.github/CONTRIBUTING.md`).

This folder is listed in `.Rbuildignore`. Neither `R CMD check` nor CI runs a
plan.
