# Test Ranking SDT model specification and integration

############################################################################# !
# MODEL CONSTRUCTOR TESTS                                                ####
############################################################################# !

test_that("sdt_ranking model can be created", {
  model <- sdt_ranking("observed", rank = "rank", m = 4)
  expect_s3_class(model, "bmmodel")
  expect_s3_class(model, "sdt")
  expect_s3_class(model, "sdt_ranking")
})

test_that("sdt_ranking has correct parameters (only dprime, no criterion)", {
  model <- sdt_ranking("observed", rank = "rank", m = 4)
  expect_true("dprime" %in% names(model$parameters))
  expect_true("mu" %in% names(model$parameters))
  expect_false("criterion" %in% names(model$parameters))
  expect_equal(model$fixed_parameters$mu, 0)
})

test_that("sdt_ranking stores m and rank correctly", {
  model <- sdt_ranking("observed", rank = "rank", m = 4)
  expect_equal(model$other_vars$m, 4L)
  expect_equal(model$other_vars$rank, "rank")

  model2 <- sdt_ranking("observed", rank = "pos", m = 3)
  expect_equal(model2$other_vars$m, 3L)
  expect_equal(model2$other_vars$rank, "pos")
})

test_that("sdt_ranking has correct links", {
  model <- sdt_ranking("observed", rank = "rank", m = 4)
  expect_equal(model$links$mu, "identity")
  expect_equal(model$links$dprime, "identity")
})

test_that("sdt_ranking accepts custom links", {
  model <- sdt_ranking("observed", rank = "rank", m = 4,
                       links = list(dprime = "log"))
  expect_equal(model$links$dprime, "log")
})

test_that("sdt_ranking requires m >= 2", {
  expect_error(sdt_ranking("observed", rank = "rank", m = 1))
})

test_that("sdt_ranking requires rank argument", {
  expect_error(sdt_ranking("observed", m = 4))
})

test_that("sdt_ranking only supports gumbel_min and normal", {
  expect_silent(sdt_ranking("observed", rank = "rank",
                            m = 4, dist = "gumbel_min"))
  expect_silent(sdt_ranking("observed", rank = "rank",
                            m = 4, dist = "normal"))
  expect_error(sdt_ranking("observed", rank = "rank",
                           m = 4, dist = "logistic"))
})

test_that("sdt_ranking with dist='normal' has sdratio as overridable fixed parameter", {
  model <- sdt_ranking("observed", rank = "rank",
                       m = 4, dist = "normal")
  expect_true("sdratio" %in% names(model$parameters))
  expect_equal(model$fixed_parameters$sdratio, 0)
})

test_that("sdt_ranking with dist='gumbel_min' does not have sdratio", {
  model <- sdt_ranking("observed", rank = "rank",
                       m = 4, dist = "gumbel_min")
  expect_false("sdratio" %in% names(model$parameters))
})

test_that("sdt_ranking has init_ranges with mu and dprime", {
  model <- sdt_ranking("observed", rank = "rank", m = 4)
  expect_false(is.null(model$init_ranges))
  expect_true(all(c("mu", "dprime") %in% names(model$init_ranges)))
  expect_equal(model$init_ranges$dprime, c(0.5, 1.5))
})

test_that("sdt_ranking normal has sdratio in init_ranges", {
  model <- sdt_ranking("observed", rank = "rank", m = 4,
                       dist = "normal")
  expect_true("sdratio" %in% names(model$init_ranges))
  expect_equal(model$init_ranges$sdratio, c(-0.3, 0.3))
})


############################################################################# !
# R-SIDE PROBABILITY TESTS                                               ####
############################################################################# !

test_that(".ranking_prob_r gives valid probabilities for Gumbel", {
  p <- bmm:::.ranking_prob_r(1.5, 1L, 4L, dist = "gumbel_min")
  expect_true(p > 0 && p < 1)
})

test_that(".ranking_all_probs_r sums to 1 for Gumbel", {
  probs <- bmm:::.ranking_all_probs_r(1.5, 4L, dist = "gumbel_min")
  expect_equal(sum(probs), 1.0, tolerance = 1e-10)
})

test_that(".ranking_all_probs_r at d'=0 gives uniform for Gumbel", {
  probs <- bmm:::.ranking_all_probs_r(0, 4L, dist = "gumbel_min")
  expect_equal(probs, rep(0.25, 4), tolerance = 1e-6)
})

test_that(".ranking_all_probs_r is monotone: rank 1 increases with dprime", {
  p_low <- bmm:::.ranking_all_probs_r(0.5, 4L, dist = "gumbel_min")
  p_high <- bmm:::.ranking_all_probs_r(2.0, 4L, dist = "gumbel_min")
  expect_true(p_high[1] > p_low[1])
})

test_that(".ranking_all_probs_r works for m=3", {
  probs <- bmm:::.ranking_all_probs_r(1.0, 3L, dist = "gumbel_min")
  expect_length(probs, 3)
  expect_equal(sum(probs), 1.0, tolerance = 1e-10)
})

test_that(".ranking_all_probs_r works for Gaussian UV-SDT", {
  probs <- bmm:::.ranking_all_probs_r(1.5, 4L, dist = "normal", sdratio = 0)
  expect_equal(sum(probs), 1.0, tolerance = 1e-6)
  expect_true(probs[1] > probs[4])
})


############################################################################# !
# DISTRIBUTION FUNCTION TESTS                                            ####
############################################################################# !

test_that("rsdt_ranking generates valid ranking data", {
  set.seed(42)
  dat <- rsdt_ranking(n_per_cell = 100, n_subjects = 5, dprime = 1.5, m = 4)

  expect_true("id" %in% colnames(dat))
  expect_true("rank" %in% colnames(dat))
  expect_true("maxrank" %in% colnames(dat))
  expect_true("observed" %in% colnames(dat))
  expect_equal(nrow(dat), 20)
  expect_true(all(dat$observed >= 0))
  expect_true(all(dat$maxrank == 4))
  expect_true(all(dat$rank %in% 1:4))
})

test_that("rsdt_ranking validates inputs", {
  expect_error(rsdt_ranking(n_per_cell = 50, n_subjects = 1, dprime = 1, m = NULL))
  expect_error(rsdt_ranking(n_per_cell = 50, n_subjects = 1, dprime = 1, m = 1))
})

test_that("dsdt_ranking computes valid ranking densities", {
  dens <- dsdt_ranking(counts = c(40, 30, 20, 10), dprime = 1.5, m = 4)
  expect_true(dens > 0)

  log_dens <- dsdt_ranking(counts = c(40, 30, 20, 10), dprime = 1.5, m = 4,
                           log = TRUE)
  expect_equal(log(dens), log_dens, tolerance = 1e-10)
})


############################################################################# !
# CHECK_DATA TESTS                                                       ####
############################################################################# !

test_that("check_data validates ranking data correctly", {
  model <- sdt_ranking("observed", rank = "rank", m = 4)
  dat <- rsdt_ranking(n_per_cell = 50, n_subjects = 3, dprime = 1.5, m = 4)
  formula <- bmf(dprime ~ 1)

  result <- check_data(model, dat, formula)
  expect_true("rank_pos" %in% colnames(result))
  expect_true("max_rank" %in% colnames(result))
  expect_equal(unique(result$max_rank), 4L)
})

test_that("check_data errors on missing rank column for ranking", {
  model <- sdt_ranking("observed", rank = "rank", m = 4)
  dat <- data.frame(id = 1:5, observed = 50, wrong_col = 1:5)
  formula <- bmf(dprime ~ 1)
  expect_error(check_data(model, dat, formula))
})

test_that("check_data errors on invalid rank values for ranking", {
  model <- sdt_ranking("observed", rank = "rank", m = 4)
  dat <- data.frame(observed = c(10, 10, 10), rank = c(0, 2, 5))
  formula <- bmf(dprime ~ 1)
  expect_error(check_data(model, dat, formula))
})


############################################################################# !
# FORMULA CONSTRUCTION TESTS                                              ####
############################################################################# !

test_that("bmf2bf.sdt_ranking produces correct formula for ranking", {
  model <- sdt_ranking("observed", rank = "rank", m = 4)
  bf <- bmf2bf(model, bmf(dprime ~ 1))

  expect_s3_class(bf, "brmsformula")
  bf_str <- deparse(bf$formula, width.cutoff = 500)
  expect_true(grepl("observed", bf_str))
  expect_true(grepl("vint.*rank_pos.*max_rank", bf_str))
  expect_false(grepl("trials", bf_str))
})


############################################################################# !
# STAN CODE GENERATION TESTS                                              ####
############################################################################# !

test_that("sdt_ranking produces valid stancode (Gumbel, m=4)", {
  dat <- rsdt_ranking(n_per_cell = 50, n_subjects = 3, dprime = 1.5, m = 4)
  model <- sdt_ranking("observed", rank = "rank", m = 4)
  formula <- bmf(dprime ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
  expect_true(grepl("sdt_ranking_lpmf", code))
  expect_true(grepl("lgamma", code))
})

test_that("sdt_ranking produces valid stancode (Gumbel, m=3)", {
  dat <- rsdt_ranking(n_per_cell = 50, n_subjects = 3, dprime = 1.5, m = 3)
  model <- sdt_ranking("observed", rank = "rank", m = 3)
  formula <- bmf(dprime ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
})

test_that("sdt_ranking stancode works with predictors on dprime", {
  set.seed(42)
  dat <- rsdt_ranking(n_per_cell = 50, n_subjects = 5, dprime = 1.5, m = 4)
  dat$condition <- rep(c("A", "B"), length.out = nrow(dat))
  model <- sdt_ranking("observed", rank = "rank", m = 4)
  formula <- bmf(dprime ~ condition)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
})

test_that("sdt_ranking stancode works with random effects on dprime", {
  dat <- rsdt_ranking(n_per_cell = 50, n_subjects = 5, dprime = 1.5, m = 4)
  model <- sdt_ranking("observed", rank = "rank", m = 4)
  formula <- bmf(dprime ~ 1 + (1 | id))
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
})
