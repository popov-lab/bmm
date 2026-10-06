# Pure computations ------------------------------------------------------------
# Targets come from closed-form conjugate normal posteriors (helper-reliability.R).
# Ratios are compared to 1 because expect_equal()'s tolerance is absolute for
# values below the target.

test_that(".rel_core() recovers the reliability of posterior means with equal trial counts", {
  withr::local_seed(11)
  sim <- conjugate_person_draws(N = 2000, tau2 = 1, sigma2 = 16, K = 30, ndraws = 400)
  core <- .rel_core(sim$draws)
  truth <- conjugate_person_reliability(1, 16, 30)
  expect_lt(abs(core$reliability / truth - 1), 0.03)
  expect_lt(core$lower, core$reliability)
  expect_gt(core$upper, core$reliability)
})

test_that(".rel_core() averages error variances, not posterior SDs, under unequal precision", {
  withr::local_seed(12)
  sim <- conjugate_person_draws(N = 4000, tau2 = 1, sigma2 = 16, K = c(5, 100), ndraws = 400)
  core <- .rel_core(sim$draws)
  truth <- conjugate_person_reliability(1, 16, sim$K)
  expect_lt(abs(core$reliability / truth - 1), 0.03)
})

test_that(".rel_core() returns person-specific reliabilities that follow each person's precision", {
  withr::local_seed(13)
  sim <- conjugate_person_draws(N = 2000, tau2 = 1, sigma2 = 16, K = c(5, 100), ndraws = 400)
  core <- .rel_core(sim$draws)
  expect_length(core$person, 2000)
  # each person's posterior variance carries Monte Carlo error, so compare the
  # average within each trial count
  by_k <- tapply(core$person, sim$K, mean)
  truth <- tapply(sim$lambda, sim$K, mean)
  expect_lt(max(abs(by_k - truth)), 0.05)
  expect_gt(by_k[["100"]], by_k[["5"]] + 0.5)
})

test_that(".rel_project_personwise() matches the closed form after multiplying the trials", {
  withr::local_seed(14)
  sim <- conjugate_person_draws(N = 4000, tau2 = 1, sigma2 = 16, K = c(5, 100), ndraws = 400)
  core <- .rel_core(sim$draws)
  expect_equal(.rel_project_personwise(core, 1)$reliability, core$reliability)
  for (f in c(0.5, 2, 4)) {
    truth <- conjugate_person_reliability(1, 16, f * sim$K)
    expect_lt(abs(.rel_project_personwise(core, f)$reliability / truth - 1), 0.03)
  }
})

test_that(".rel_retest() recovers the true-score correlation of a joint model", {
  withr::local_seed(15)
  N <- 2000
  ndraws <- 300
  r <- 0.6
  prior <- matrix(c(1, r, r, 1), 2)
  s2 <- 0.5
  post_cov <- solve(solve(prior) + diag(1 / s2, 2))
  L <- chol(post_cov)
  true <- matrix(stats::rnorm(2 * N), N) %*% chol(prior)
  theta1 <- theta2 <- matrix(0, ndraws, N)
  for (i in seq_len(N)) {
    y <- true[i, ] + stats::rnorm(2, 0, sqrt(s2))
    mu <- post_cov %*% (y / s2)
    d <- sweep(matrix(stats::rnorm(2 * ndraws), ndraws) %*% L, 2, mu, `+`)
    theta1[, i] <- d[, 1]
    theta2[, i] <- d[, 2]
  }
  rt <- .rel_retest(theta1, theta2)
  expect_lt(abs(stats::median(rt$stability) / r - 1), 0.08)
  expect_lt(abs(stats::median(rt$agreement) / r - 1), 0.08)

  shifted <- .rel_retest(theta1, theta2 + 1)
  expect_equal(shifted$stability, rt$stability)
  expect_lt(stats::median(shifted$agreement), stats::median(rt$agreement) - 0.2)
})

test_that(".rel_gstudy() recovers G, the error variance and the D-study for persons x occasions", {
  withr::local_seed(16)
  sp2 <- 1
  spo2 <- 0.5
  s2 <- 0.8
  sim <- conjugate_facet_draws(N = 1500, sp2 = sp2, spo2 = spo2, s2 = s2, n_o = 3, ndraws = 300)
  interactions <- list(list(facets = "session", var = rep(spo2, 300)))
  g <- .rel_gstudy(sim$universe, interactions, list(), c(session = 3))
  expect_lt(abs(g$g_relative / sim$g_true - 1), 0.03)
  expect_lt(abs(g$s2 / s2 - 1), 0.15)

  same <- .rel_dstudy(g, c(session = 3), 1)
  expect_equal(same$g_relative, g$g_relative)

  truth <- sp2 / (sp2 + spo2 / 6 + s2 / (2 * 6))
  expect_lt(abs(.rel_dstudy(g, c(session = 6), 2)$g_relative / truth - 1), 0.03)
})

test_that(".rel_gstudy() adds facet main effects to the absolute error only", {
  withr::local_seed(17)
  sim <- conjugate_facet_draws(N = 1500, sp2 = 1, spo2 = 0.5, s2 = 0.8, n_o = 3, ndraws = 300)
  interactions <- list(list(facets = "session", var = rep(0.5, 300)))
  mains <- list(list(facets = "session", var = rep(0.3, 300)))
  g <- .rel_gstudy(sim$universe, interactions, mains, c(session = 3))
  expect_lt(abs(g$g_relative / sim$g_true - 1), 0.03)
  truth_abs <- 1 / (1 + (0.5 + 0.8) / 3 + 0.3 / 3)
  expect_lt(abs(g$g_absolute / truth_abs - 1), 0.03)
  expect_lt(g$g_absolute, g$g_relative)
})

test_that(".rel_gstudy() divides the error by the crossing of all person-facet terms", {
  withr::local_seed(18)
  # persons x sessions x stimuli, both crossed with persons: the universe-score
  # error is spo2 / n_o + spk2 / n_k + s2 / (n_o * n_k)
  spo2 <- 0.4
  spk2 <- 0.2
  s2 <- 3
  n <- c(session = 2, stimulus = 10)
  delta <- spo2 / 2 + spk2 / 10 + s2 / 20
  sim <- conjugate_person_draws(N = 3000, tau2 = 1, sigma2 = delta, K = 1, ndraws = 300)
  interactions <- list(
    list(facets = "session", var = rep(spo2, 300)),
    list(facets = "stimulus", var = rep(spk2, 300))
  )
  g <- .rel_gstudy(sim$draws, interactions, list(), n)
  expect_lt(abs(g$s2 / s2 - 1), 0.15)
  truth <- 1 / (1 + spo2 / 4 + spk2 / 30 + s2 / (2 * 4 * 30))
  expect_lt(abs(.rel_dstudy(g, c(session = 4, stimulus = 30), 2)$g_relative / truth - 1), 0.03)
})

test_that(".rel_solve_factor() finds the trial multiplier that reaches a target", {
  project <- function(f) 1 - 0.4 / (1 + f)
  f <- .rel_solve_factor(project, 0.9)
  expect_equal(project(f), 0.9, tolerance = 1e-6)
  expect_equal(.rel_solve_factor(function(f) 0.5 + 0 * f, 0.9), Inf)
})


# Trial counts -----------------------------------------------------------------

test_that("person_trials() returns one positive count per row for every model", {
  skip_on_cran()
  for (case_name in names(stored_frame_cases())) {
    fit <- stored_frame_fit(stored_frame_cases()[[case_name]])
    trials <- person_trials(fit$bmm$model, fit$data, fit$formula)
    expect(
      length(trials) == nrow(fit$data) && all(is.finite(trials) & trials > 0),
      glue::glue("person_trials() fails for the stored-frame case '{case_name}'")
    )
  }
})

test_that("person_trials() reads the counts of aggregated data", {
  skip_on_cran()
  fit <- stored_frame_fit(stored_frame_cases()$m3)
  expect_equal(
    person_trials(fit$bmm$model, fit$data, fit$formula),
    as.numeric(rowSums(fit$data$Y))
  )
  ez <- ezdm("mean_rt", "var_rt", "n_upper", "n_trials", version = "4par")
  expect_equal(person_trials(ez, data.frame(n_trials = c(40, 60)), NULL), c(40, 60))
})


# Fitted models ----------------------------------------------------------------

reliability_fixture <- function(name) {
  path <- test_path("assets", name)
  skip_if_not(file.exists(path), "reliability fixture not available (excluded by .Rbuildignore)")
  restructure(readRDS(path))
}

test_that("bmm_reliability() reports cells, a G-study and person tables for a trial-level model", {
  skip_on_cran()
  fit <- reliability_fixture("bmmfit_sdm_reliability.rds")
  rel <- bmm_reliability(fit)
  expect_s3_class(rel, "bmm_reliability")
  expect_equal(attr(rel, "group"), "id")
  expect_equal(sort(unique(rel$parameter)), c("c", "kappa"))

  kappa <- rel[rel$parameter == "kappa", ]
  expect_equal(nrow(kappa), 1)
  # the G coefficient of the design as fitted is the universe-score reliability
  expect_equal(kappa$g_relative, kappa$reliability)
  expect_true(kappa$g_absolute <= kappa$g_relative)
  expect_equal(kappa$n_session, 2)
  components <- attr(rel, "components")
  expect_equal(components$component, c("id", "id:session", "error"))

  c_rows <- rel[rel$parameter == "c", ]
  expect_equal(as.character(c_rows$cond), c("A", "B"))
  expect_true(all(c_rows$reliability > 0 & c_rows$reliability < 1))
  expect_equal(c_rows$trials_median, c(40, 40))
  expect_true(all(is.na(c_rows$g_relative)))

  persons <- attr(rel, "persons")
  expect_equal(nrow(persons), 30 * 3)
  expect_true(all(persons$reliability <= 1))
})

test_that("bmm_reliability() handles contrasts, retest pairs, scales and fixed parameters", {
  skip_on_cran()
  fit <- reliability_fixture("bmmfit_sdm_reliability.rds")

  contrast <- bmm_reliability(fit, pars = "c", contrast = "cond")
  expect_equal(contrast$contrast, "B - A")
  expect_equal(attr(contrast, "type"), "contrast")

  retest <- bmm_reliability(fit, pars = "c", retest = "cond")
  expect_equal(nrow(retest), 1)
  expect_true(retest$stability > -1 && retest$stability < 1)
  expect_equal(retest$expected_retest, retest$stability * sqrt(retest$reliability_1 * retest$reliability_2))

  native <- bmm_reliability(fit, pars = "kappa", scale = "native")
  expect_false("g_relative" %in% names(native))
  expect_equal(bmm_reliability(fit, pars = "kappa", scale = "sampling")$scale, "link")

  per_session <- bmm_reliability(fit, pars = "kappa", generalizability = FALSE)
  expect_equal(nrow(per_session), 2)
  expect_false("g_relative" %in% names(per_session))

  expect_match(bmm_reliability(fit, pars = "mu")$note, "fixed")
  expect_error(bmm_reliability(fit, pars = "c", contrast = "session"), "not a predictor")
  expect_error(bmm_reliability(fit, group = "nope"), "'group' must be one of")
})

test_that("reliability_for_design() projects trials and facet levels", {
  skip_on_cran()
  fit <- reliability_fixture("bmmfit_sdm_reliability.rds")
  rel <- bmm_reliability(fit)

  design <- reliability_for_design(rel, trials = c(1, 2, 4))
  c_a <- design[design$parameter == "c" & design$cond == "A", ]
  expect_equal(c_a$reliability[1], rel$reliability[rel$parameter == "c" & rel$cond == "A"])
  expect_true(all(diff(c_a$reliability) > 0))
  expect_equal(c_a$trials_per_person, c(40, 80, 160))

  kappa <- reliability_for_design(rel, session = c(2, 6))
  kappa <- kappa[kappa$parameter == "kappa", ]
  expect_equal(kappa$reliability[1], rel$g_relative[rel$parameter == "kappa"])
  expect_gt(kappa$reliability[2], kappa$reliability[1])

  target <- reliability_for_design(rel, target = 0.9)
  c_target <- target[target$parameter == "c", ]
  for (i in seq_len(nrow(c_target))) {
    at <- reliability_for_design(rel, trials = c_target$trials[i])
    at <- at[at$parameter == "c" & at$cond == c_target$cond[i], ]
    expect_equal(at$reliability, 0.9, tolerance = 1e-6)
  }
  kept <- reliability_for_design(rel[rel$parameter == "c", ])
  expect_equal(unique(kept$parameter), "c")
  expect_error(reliability_for_design(rel[0, ]), "no rows")
  expect_error(reliability_for_design(rel, stimulus = 4), "not a facet")
  expect_error(reliability_for_design(rel, trials = 0), "positive")
  expect_error(
    reliability_for_design(bmm_reliability(fit, pars = "c", retest = "cond")),
    "not test-retest"
  )
})

test_that("bmm_reliability() works on aggregated data and counts its trials", {
  skip_on_cran()
  fit <- reliability_fixture("bmmfit_m3_reliability.rds")
  rel <- bmm_reliability(fit, pars = c("c", "a"))
  expect_equal(attr(rel, "group"), "ID")
  c_rows <- rel[rel$parameter == "c", ]
  expect_equal(nrow(c_rows), length(unique(fit$data$cond)))
  first <- fit$data[fit$data$cond == c_rows$cond[1], ]
  expect_equal(c_rows$trials_median[1], stats::median(tapply(rowSums(first$Y), first$ID, sum)))
  expect_true(all(is.finite(rel$snr)))
  # a ~ 1 pools all conditions, so every trial of a person counts
  expect_equal(
    rel$trials_median[rel$parameter == "a"],
    stats::median(tapply(rowSums(fit$data$Y), fit$data$ID, sum))
  )
})

test_that("bmm_reliability() rejects bad arguments at the boundary", {
  skip_on_cran()
  fit <- reliability_fixture("bmmfit_sdm_reliability.rds")
  expect_error(bmm_reliability(fit, contrast = c("cond", "session")), "one predictor")
  expect_error(bmm_reliability(fit, contrast = "cond", retest = "cond"), "not both")
  expect_error(bmm_reliability(fit, generalizability = "yes"), "TRUE or FALSE")
  expect_error(bmm_reliability(fit, re_formula = NA), "remove 're_formula'")
  expect_error(bmm_reliability(fit, prob = 1.2), "'prob'")
  expect_error(bmm_reliability(fit, ndraws = 1), "'ndraws'")
  expect_error(bmm_reliability(fit, pars = "nope"), "Unknown parameter")
})

test_that("bmm_reliability() needs at least three persons per cell", {
  skip_on_cran()
  fit <- reliability_fixture("bmmfit_sdm_reliability.rds")
  nd <- data.frame(cond = factor("A", levels = c("A", "B")), id = factor(1:2, levels = 1:30))
  expect_error(bmm_reliability(fit, pars = "c", newdata = nd), "at least 3")
})

test_that("generalizability needs the link scale, also for design projections", {
  skip_on_cran()
  fit <- reliability_fixture("bmmfit_sdm_reliability.rds")
  native <- bmm_reliability(fit, pars = "kappa", scale = "native")
  expect_match(native$note, "scale = 'link'")
  expect_null(attr(native, "components"))
  expect_error(reliability_for_design(native, session = 4), "link scale only")
})

test_that("print.bmm_reliability() labels the scale and explains the columns", {
  skip_on_cran()
  fit <- reliability_fixture("bmmfit_sdm_reliability.rds")
  expect_output(print(bmm_reliability(fit)), "scale: link.*g_relative / g_absolute")
  expect_output(print(bmm_reliability(fit, pars = "c", contrast = "cond")), "B - A")
  expect_output(print(bmm_reliability(fit, pars = "c", retest = "cond")), "stability: correlation")
})
