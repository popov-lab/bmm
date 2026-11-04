# 0006: Trial-level data for the count models

- **Status:** accepted 2026-11-02
- **Date:** 2026-11-02
- **Amends:** [0001](0001-release-scope-1.4-to-1.6.md)

## Context

Item random effects and trial-level covariates are impossible on aggregated counts, so users fit signal detection as a probit GLMM in raw brms. A spike (2026-10-15, `develop`, brms 2.23.0, cmdstanr, 2 chains of 1000 iterations) found that every count model (`sdt_yn`, `sdt_mafc`, `sdt_rating`, `sdt_ranking`, `sdt_cdp`, `mpt`) already accepts one row per trial (`n_trials` of 1s, or one-hot count columns with `nTrials = 1`) and returns the aggregated posterior within Monte Carlo error. Trial-level `sdt_yn` matches the brms probit GLMM's item SD of `d` (0.304 vs 0.295). Missing are a user contract, data checks, documentation, and speed for rating-type models.

Per-row cost is the same for a cell and a trial row: `sdt_yn` 0.90–1.05 µs, `mpt` (two categories) 1.25–1.8 µs, `sdt_rating` 7.8–9.2 µs, brms `bernoulli("probit")` 0.22 µs. `brms::categorical()` was slower than one-hot `multinomial()`, 13.9 vs 11.0 ms per gradient.

GitHub, read 2026-11-02: no open or closed issue asks for trial-level data, item effects or a GLMM workflow. Milestone 1.5.0 has 23 open issues, due 2026-12-15; 1.6.0 has 12, due 2027-02-15. bmm 1.4.1 is on CRAN, `dev-1.5.0` has merged into `develop`, and `dev-1.6.0` is the only integration branch.

## Decision

Decided by GF, 2026-11-02:

- Rating argument: `n_ratings`.
- `mpt()`: new `response` argument.
- `sdt_cdp`: `judgment` and `confidence` arguments.
- 1.4.1 is not a vehicle: it is released.
- Release: T1–T5 aim for 1.5.0, with T1–T3 PRs ready by 2026-12-01. Each merges after an independent agent review and green CI, the summary posted in the PR and 72 hours open (CONTRIBUTING, one active maintainer; new constructor arguments are the exported-API tier). A PR that misses this is retargeted to `dev-1.6.0` (0002). T2a stays in 1.6.0.

**D1. The user declares trial rows at construction, through the arguments that already carry the layout**. No data sniffing, no `*_trial()` constructors. `sdt_yn` and `sdt_mafc` take `n_trials = NULL` with a 0/1 `response`; `sdt_rating` takes one column with values 1..K plus `n_ratings = K`; `sdt_ranking` takes one rank column; `mpt` takes `response =` holding branch names per trial.

**D2. The families stay.** `categorical()` is rejected. A trial-native `sdt_rating_lpmf` is built only if a spike (T2a) measures at least 5x per row. A vectorised `loop = FALSE` likelihood for `sdt_yn` and `sdt_mafc` is T5, gate at least 2x. Threading already works and is documented.

**D3. One shared layer:** `.trial_rows_to_counts()`, a `trial_level` marker in `other_vars` read by `check_data()`, the requirements text and `bmm_data_check()` (which then reports observations per subject and, for items, only min / median / max, and judges sparseness on participant cells), one shared documentation section, and one parity test in `test-conventions.R`. Each help page adds it with one `model_docs()` line.

**Order and targets** (T0–T5; feature PRs target the branch of their release, 0002):

| Step | Content | Target |
|---|---|---|
| T0 | this record, issue draft | `develop` (`.dev/` only) |
| T1 | layer, `sdt_yn`, `sdt_mafc` | `develop` (1.5.0) |
| T4 | probit-GLMM article, after T1 | `develop` (1.5.0) |
| T5 | `loop = FALSE` likelihood, serialised with T1 on `R/model_sdt_yn.R` | `develop` (1.5.0) |
| T2 | `sdt_rating`, `sdt_ranking`, `sdt_cdp`, after T1 | `develop`, else `dev-1.6.0` |
| T3 | `mpt(response =)`, after T1, parallel with T2, never the same file at once | `develop`, else `dev-1.6.0` |
| T2a | rating-family spike, after T2 | `dev-1.6.0` |

T1 and T4 make `sdt_yn()` the drop-in for the probit GLMM and fit 0001's "smaller items" slot. T2 and T3 ride along only if ready in time; 1.5.0 is full.

## Consequences

- 1.5.0 gains up to five items and no new model; 1.6.0 gains T2a and any of T2, T3 that miss 1 December.
- A trial-level `sdt_rating` fit costs about ten times a trial-level `sdt_yn` fit per row; `threads` is the lever until T2a passes.
- `bmm_data_check()` counts rows, not trials, per cell, so aggregated `sdt_*` cells of few rows with many trials are flagged sparse (measured: 2 rows of 50 trials). T1 leaves the aggregated report unchanged; the fix is a separate PR to `develop`.
- `loo()` on trials and cells are different estimands; use `kfold(group = "id")`. Open: `ppp()` needs `participant =` when a fit has two grouping variables; the `loo()` worker error on the trial-level rating fit is unchecked outside `load_all()`.
