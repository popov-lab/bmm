# SDT models on real fits (moved from tests/testthat in #422)

skip_if_not_installed("cmdstanr")
skip_if(is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)), "CmdStan is not installed")

ranks4 <- paste0("rank", 1:4)

sim_ranking <- function(n, n_trials, m, d) {
  dat <- data.frame(id = seq_len(n), set_size = rep_len(as.integer(m), n))
  cbind(dat, as.data.frame(rsdt_ranking(n, n_trials, m = m, d = d)))
}

sim_rating <- function(n_subjects, n_trials, d, criterion, n_ratings, spacing = NULL,
                       deltas = NULL, threshold_type = "parsimonious") {
  thr <- .sdt_make_thresholds(criterion, n_ratings, threshold_type, spacing, deltas)
  dat <- expand.grid(id = seq_len(n_subjects), stimulus = c(0L, 1L))
  cbind(dat, as.data.frame(rsdt_rating(nrow(dat), n_trials, dat$stimulus, d, thr)))
}

test_that("sdt_ranking recovers d' and returns a multinomial count array", {
  withr::local_seed(422)
  dat <- sim_ranking(10, 200, m = 4, d = 1.4)
  model <- sdt_ranking(ranks4, m = 4)
  fit <- bmm(bmf(d ~ 1 + (1 | id)), dat, model, backend = "cmdstanr",
             chains = 2, iter = 600, warmup = 300, init = 0.5,
             refresh = 0, silent = 2)

  est <- brms::fixef(fit)["d_Intercept", "Estimate"]
  expect_equal(est, 1.4, tolerance = 0.3)

  pp <- brms::posterior_predict(fit, ndraws = 25)
  expect_length(dim(pp), 3)
  expect_equal(dim(pp)[3], 4L)
  expect_false(any(is.nan(pp)))
})

test_that("sdt_rating posterior_predict draws joint multinomial counts (real fit)", {
  withr::local_seed(801)
  dat <- sim_rating(6, 120, d = 1.4, criterion = 0, n_ratings = 4,
                    spacing = 0.5)
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  fit <- bmm(bmf(d ~ 1, criterion ~ 1, spacing ~ 1), dat, model,
             backend = "cmdstanr", chains = 1, iter = 300, warmup = 150,
             refresh = 0, silent = 2)

  pp <- brms::posterior_predict(fit)
  expect_length(dim(pp), 3)
  row_sums <- apply(pp, c(1, 2), sum)
  expect_true(all(row_sums == 120))
})

test_that("sdt_rating log_distance fits with finite likelihood at K=6 (real fit)", {
  # Regression: at K=6 the log_distance threshold builder has mid=3, so its
  # lower thresholds were filled by a descending Stan for-loop. Stan for-loops
  # only count up, so the loop was skipped and thresholds[1:2] stayed NaN,
  # failing every chain. Earlier tests only used K=4 (mid=2), where the loop ran.
  withr::local_seed(820)
  dat <- sim_rating(6, 150, d = 1.4, criterion = 0, n_ratings = 6,
                    deltas = c(-0.4, -0.4, -0.4, -0.4),
                    threshold_type = "log_distance")
  model <- sdt_rating(paste0("r", 1:6), "stimulus", threshold_type = "log_distance")
  fit <- bmm(bmf(d ~ 1, criterion ~ 1, delta1 ~ 1, delta2 ~ 1,
                 delta3 ~ 1, delta4 ~ 1),
             dat, model, backend = "cmdstanr", chains = 1, iter = 300, warmup = 150,
             refresh = 0, silent = 2)

  expect_true(all(is.finite(brms::log_lik(fit))))
})
