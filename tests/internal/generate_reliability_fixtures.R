# Generates the fitted-model fixtures used by tests/testthat/test-reliability.R:
# a trial-level model with persons crossed with sessions (sdm) and an
# aggregated-data model with person effects (m3). Run from the repo root:
#   Rscript tests/internal/generate_reliability_fixtures.R
# Fixtures are git-tracked in tests/testthat/assets/ but excluded from the
# built package via .Rbuildignore; tests skip when they are absent.

devtools::load_all(quiet = TRUE)

fit_args <- list(
  backend = "cmdstanr", chains = 2, iter = 500, warmup = 300,
  refresh = 0, silent = 2, seed = 123
)

save_fixture <- function(fit, name) {
  path <- file.path("tests/testthat/assets", name)
  saveRDS(fit, path)
  cat(name, ":", round(file.size(path) / 1024), "KB\n")
}

# sdm: 30 persons x 2 sessions x 2 conditions x 20 trials. kappa varies over
# persons and person x session (a G-study facet); c differs between conditions
# with person-specific effects (cells, contrasts, retest pairs)
set.seed(10)
n_id <- 30
grid <- expand.grid(session = factor(1:2), cond = factor(c("A", "B")), id = factor(1:n_id))
u_kappa <- stats::rnorm(n_id, 0, 0.35)
u_kappa_session <- stats::rnorm(nrow(grid), 0, 0.2)
u_c <- MASS::mvrnorm(n_id, c(0, 0), matrix(c(0.25, 0.15, 0.15, 0.25), 2))
sdm_data <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  id <- as.integer(grid$id[i])
  cond_b <- grid$cond[i] == "B"
  kappa <- exp(log(8) + u_kappa[id] + u_kappa_session[i])
  c <- exp(log(3) - 0.4 * cond_b + u_c[id, 1 + cond_b])
  data.frame(grid[rep(i, 20), ], y = rsdm(20, c = c, kappa = kappa))
}))
sdm_fit <- do.call(bmm, c(list(
  formula = bmf(
    kappa ~ 1 + (1 | id) + (1 | id:session),
    c ~ 0 + cond + (0 + cond | id)
  ),
  data = sdm_data,
  model = sdm(resp_error = "y")
), fit_args))
save_fixture(sdm_fit, "bmmfit_sdm_reliability.rds")

# m3: aggregated counts per person and condition with person effects on c and a
m3_cats <- c("corr", "other", "dist", "npl")
m3_fit <- do.call(bmm, c(list(
  formula = bmf(
    corr ~ b + a + c, other ~ b + a, dist ~ b + d, npl ~ b,
    c ~ 1 + cond + (1 | ID), a ~ 1 + (1 | ID), d ~ 1
  ),
  data = oberauer_lewandowsky_2019_e1,
  model = m3(
    resp_cats = m3_cats, num_options = paste0("n_", m3_cats),
    choice_rule = "simple", links = list(c = "log", a = "log", d = "log")
  )
), fit_args))
save_fixture(m3_fit, "bmmfit_m3_reliability.rds")
