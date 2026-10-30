# Recovery plans

One file per model: `tests/recovery/<model>.R`. It builds a validation plan with
`bmmtools::validation_plan()` and returns it as its last expression. The plan
says what is simulated, how it is fitted and what counts as recovered; it does
not fit anything when it is sourced.

```r
# tests/recovery/mixture2p.R
bmmtools::validation_plan(
  "mixture2p",
  seed = 2026
)
```

A model without a bmmtools defaults entry (0.3.0 has `mixture2p` and `ezdm`)
needs `pars =`, the population values on the link scale, and usually `sds =`
and `grid =` as well. The plan must cover a true-against-recovered comparison,
a hierarchical design with several subjects, and data from an independent
generator (see `.github/CONTRIBUTING.md`).

To run a plan and render its card (`local/` is untracked):

```r
plan <- source("tests/recovery/mixture2p.R")$value
result <- bmmtools::run_validation(plan, dir = "local/recovery/mixture2p")
bmmtools::render_card(result, output = "local/recovery/mixture2p/card.html")
```

`bmmtools::write_plan()` writes a plan to an `.rds` file. Commit the `.R`
script instead: it can be read in a diff.

This folder is listed in `.Rbuildignore`. Neither `R CMD check` nor CI runs a
plan.
