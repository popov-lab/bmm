# Test SDAI model specification and integration

############################################################################# !
# MODEL CONSTRUCTOR TESTS                                                ####
############################################################################# !

test_that("sdai model can be created with default arguments", {
  resp_cols <- c("i1_c1", "i1_c2", "i1_c3", "i1_c4",
                 "i2_c1", "i2_c2", "i2_c3", "i2_c4")
  expect_silent(sdai(resp = resp_cols))
})

test_that("sdai model can be created with both distributions", {
  resp_cols <- c("i1_c1", "i1_c2", "i1_c3", "i1_c4",
                 "i2_c1", "i2_c2", "i2_c3", "i2_c4")
  expect_silent(sdai(resp = resp_cols, dist = "gumbel_min"))
  expect_silent(sdai(resp = resp_cols, dist = "normal"))
})

test_that("sdai model has correct class structure", {
  resp_cols <- c("i1_c1", "i1_c2", "i1_c3", "i1_c4",
                 "i2_c1", "i2_c2", "i2_c3", "i2_c4")
  model <- sdai(resp = resp_cols, dist = "gumbel_min")
  expect_s3_class(model, "bmmodel")
  expect_s3_class(model, "sdai")
  expect_s3_class(model, "sdai_gumbel")

  model_g <- sdai(resp = resp_cols, dist = "normal")
  expect_s3_class(model_g, "sdai_gaussian")
})

test_that("sdai model parameters are correctly defined", {
  resp_cols <- c("i1_c1", "i1_c2", "i1_c3", "i1_c4",
                 "i2_c1", "i2_c2", "i2_c3", "i2_c4")
  model <- sdai(resp = resp_cols, n_confidence = 4)
  expect_true("mu" %in% names(model$parameters))
  expect_true("dprime" %in% names(model$parameters))
  expect_true("criterion" %in% names(model$parameters))
  expect_equal(model$fixed_parameters$mu, 0)
  # n_confidence=4 -> n_thresholds=3 -> n_margins=2
  expect_true("crla" %in% names(model$parameters))
  expect_true("crlb" %in% names(model$parameters))
  expect_true("crha" %in% names(model$parameters))
  expect_true("crhb" %in% names(model$parameters))
})

test_that("sdai gaussian model includes sdratio", {
  resp_cols <- c("i1_c1", "i1_c2", "i1_c3", "i1_c4",
                 "i2_c1", "i2_c2", "i2_c3", "i2_c4")
  model <- sdai(resp = resp_cols, dist = "normal")
  expect_true("sdratio" %in% names(model$parameters))
})

test_that("sdai model has correct default links", {
  resp_cols <- c("i1_c1", "i1_c2", "i1_c3", "i1_c4",
                 "i2_c1", "i2_c2", "i2_c3", "i2_c4")
  model <- sdai(resp = resp_cols, n_confidence = 4)
  expect_equal(model$links$mu, "identity")
  expect_equal(model$links$dprime, "identity")
  expect_equal(model$links$criterion, "identity")
  expect_equal(model$links$crla, "log")
  expect_equal(model$links$crha, "log")
})

test_that("sdai model rejects invalid arguments", {
  resp_cols <- c("i1_c1", "i1_c2", "i1_c3", "i1_c4",
                 "i2_c1", "i2_c2", "i2_c3", "i2_c4")
  expect_error(sdai(resp = resp_cols, n_sources = 1))
  expect_error(sdai(resp = resp_cols, n_confidence = 1))
  expect_error(sdai(resp = resp_cols, dist = "cauchy"))
})

test_that("sdai model with n_confidence=2 has no margin parameters", {
  resp_cols <- c("i1_c1", "i1_c2", "i2_c1", "i2_c2")
  model <- sdai(resp = resp_cols, n_confidence = 2)
  expect_false("crla" %in% names(model$parameters))
  expect_false("crha" %in% names(model$parameters))
  expect_equal(model$other_vars$n_margins, 0)
})

test_that("sdai with 3 sources creates correct structure", {
  resp_cols <- paste0("i", rep(1:3, each = 4), "_c", rep(1:4, 3))
  model <- sdai(resp = resp_cols, n_sources = 3, n_confidence = 4)
  expect_equal(model$other_vars$n_sources, 3)
  expect_equal(model$other_vars$n_confidence, 4)
})


############################################################################# !
# CHECK_DATA TESTS                                                       ####
############################################################################# !

test_that("sdai check_data validates response columns", {
  resp_cols <- c("i1_c1", "i1_c2", "i1_c3", "i1_c4",
                 "i2_c1", "i2_c2", "i2_c3", "i2_c4")
  model <- sdai(resp = resp_cols)
  dat <- data.frame(
    i1_c1 = c(10, 5), i1_c2 = c(20, 15), i1_c3 = c(15, 20), i1_c4 = c(5, 10),
    i2_c1 = c(8, 12), i2_c2 = c(18, 22), i2_c3 = c(22, 18), i2_c4 = c(12, 8)
  )
  formula <- bmf(dprime ~ 1, criterion ~ 1, crla ~ 1, crha ~ 1,
                 crlb ~ 1, crhb ~ 1)
  result <- check_data(model, dat, formula)
  expect_true(!is.null(attr(result, "resp_cols")))
  expect_true(!is.null(attr(result, "resp_mat")))
})

test_that("sdai check_data rejects missing columns", {
  resp_cols <- c("i1_c1", "i1_c2", "i1_c3", "i1_c4",
                 "i2_c1", "i2_c2", "i2_c3", "i2_c4")
  model <- sdai(resp = resp_cols)
  dat <- data.frame(i1_c1 = 10, i1_c2 = 20)
  formula <- bmf(dprime ~ 1, criterion ~ 1, crla ~ 1, crha ~ 1,
                 crlb ~ 1, crhb ~ 1)
  expect_error(check_data(model, dat, formula), "not found")
})

test_that("sdai check_data rejects negative counts", {
  resp_cols <- c("i1_c1", "i1_c2", "i1_c3", "i1_c4",
                 "i2_c1", "i2_c2", "i2_c3", "i2_c4")
  model <- sdai(resp = resp_cols)
  dat <- data.frame(
    i1_c1 = c(-1), i1_c2 = c(20), i1_c3 = c(15), i1_c4 = c(5),
    i2_c1 = c(8), i2_c2 = c(18), i2_c3 = c(22), i2_c4 = c(12)
  )
  formula <- bmf(dprime ~ 1, criterion ~ 1, crla ~ 1, crha ~ 1,
                 crlb ~ 1, crhb ~ 1)
  expect_error(check_data(model, dat, formula), "negative")
})

test_that("sdai check_data rejects wrong number of columns", {
  resp_cols <- c("i1_c1", "i1_c2", "i1_c3")
  model <- sdai(resp = resp_cols, n_sources = 2, n_confidence = 4)
  dat <- data.frame(i1_c1 = 10, i1_c2 = 20, i1_c3 = 15)
  formula <- bmf(dprime ~ 1, criterion ~ 1, crla ~ 1, crha ~ 1,
                 crlb ~ 1, crhb ~ 1)
  expect_error(check_data(model, dat, formula), "Expected 8")
})


############################################################################# !
# R-SIDE PROBABILITY COMPUTATION TESTS                                   ####
############################################################################# !

test_that("Gumbel tabs probabilities sum to 1", {
  thresholds <- matrix(c(-1, 0, 1), nrow = 1)
  probs <- .sdai_tabs_probs_gumbel(thresholds)
  expect_equal(sum(probs), 1, tolerance = 1e-6)
})

test_that("Gumbel tabs probabilities sum to 1 with multiple draws", {
  thresholds <- matrix(c(-1, 0, 1, -0.5, 0.5, 1.5), nrow = 2, byrow = TRUE)
  probs <- .sdai_tabs_probs_gumbel(thresholds)
  expect_equal(rowSums(probs), c(1, 1), tolerance = 1e-6)
})

test_that("Gumbel tpres probabilities sum to 1 across all cells", {
  thresholds <- matrix(c(-1, 0, 1), nrow = 1)
  dprime <- 1.5
  probs <- .sdai_tpres_probs_gumbel(thresholds, dprime)
  # All 2*n_cats probabilities form one simplex (one multinomial)
  expect_equal(sum(probs[1, ]), 1, tolerance = 1e-4)
})

test_that("Gaussian tabs probabilities sum to 1", {
  thresholds <- matrix(c(-1, 0, 1), nrow = 1)
  probs <- .sdai_tabs_probs_gaussian(thresholds)
  expect_equal(sum(probs), 1, tolerance = 1e-6)
})

test_that("Gaussian tpres probabilities sum to 1 across all cells", {
  thresholds <- matrix(c(-1, 0, 1), nrow = 1)
  dprime <- 1.5
  sdratio <- 1.0
  probs <- .sdai_tpres_probs_gaussian(thresholds, dprime, sdratio)
  expect_equal(sum(probs[1, ]), 1, tolerance = 1e-4)
})

test_that("Gumbel pfun returns valid probabilities", {
  p <- .sdai_pfun_gumbel(-100, 100, 0, 1)
  expect_true(p > 0 && p <= 1)
})

test_that("SDAI thresholds are ordered correctly", {
  criterion <- c(0)
  margins_l <- matrix(c(0.5, 0.3), nrow = 1)
  margins_h <- matrix(c(0.5, 0.3), nrow = 1)
  thres <- .sdai_thresholds(criterion, margins_l, margins_h)
  expect_true(all(diff(thres[1, ]) > 0))
})

test_that("SDAI thresholds with single threshold equals criterion", {
  criterion <- c(0.5)
  thres <- .sdai_thresholds(criterion, NULL, NULL)
  expect_equal(thres[1, 1], 0.5)
})

test_that("higher dprime increases source discrimination in Gumbel", {
  thresholds <- matrix(c(-1, 0, 1), nrow = 1)
  probs_low <- .sdai_tpres_probs_gumbel(thresholds, dprime = 0.5)
  probs_high <- .sdai_tpres_probs_gumbel(thresholds, dprime = 2.0)
  n_cats <- 4
  # With higher dprime, source 2 (g1=dprime, g2=0) concentrates in the

  # highest category (col 2*n_cats) where it dominates
  expect_true(probs_high[1, 2 * n_cats] > probs_low[1, 2 * n_cats])
})


############################################################################# !
# STAN CODE GENERATION TESTS                                             ####
############################################################################# !

test_that("Gumbel lpmf generation produces valid Stan code", {
  code <- .sdai_generate_gumbel_lpmf(n_confidence = 4, n_sources = 2,
                                      n_margins = 2)
  expect_true(grepl("sdai_gumbel_lpmf", code))
  expect_true(grepl("multinomial_lpmf", code))
  expect_true(grepl("sdai_pfun", code))
  expect_true(grepl("sdai_p_hit_gumbel", code))
  expect_true(grepl("thres\\[", code))
})

test_that("Gaussian lpmf generation produces valid Stan code", {
  code <- .sdai_generate_gaussian_lpmf(n_confidence = 4, n_sources = 2,
                                        n_margins = 2)
  expect_true(grepl("sdai_gaussian_lpmf", code))
  expect_true(grepl("multinomial_lpmf", code))
  expect_true(grepl("integrate_1d", code))
  expect_true(grepl("sdai_int_inst_uvg", code))
  expect_true(grepl("sdratio", code))
})

test_that("Gumbel lpmf with n_confidence=2 has no margin parameters", {
  code <- .sdai_generate_gumbel_lpmf(n_confidence = 2, n_sources = 2,
                                      n_margins = 0)
  expect_false(grepl("crl", code))
  expect_false(grepl("crh", code))
  expect_true(grepl("thres\\[1\\] = criterion", code))
})

test_that("Generated lpmf has correct number of vint args", {
  # 2 sources x 4 confidence = 8 total, minus 1 for y = 7 vint args
  code <- .sdai_generate_gumbel_lpmf(n_confidence = 4, n_sources = 2,
                                      n_margins = 2)
  expect_true(grepl("int y7", code))
  expect_false(grepl("int y8", code))
})


############################################################################# !
# CONFIGURE_MODEL / PIPELINE TESTS                                       ####
############################################################################# !

test_that("sdai Gumbel produces valid stancode", {
  resp_cols <- c("i1_c1", "i1_c2", "i1_c3", "i1_c4",
                 "i2_c1", "i2_c2", "i2_c3", "i2_c4")
  model <- sdai(resp = resp_cols, n_sources = 2, n_confidence = 4)
  dat <- data.frame(
    i1_c1 = c(10, 5), i1_c2 = c(20, 15), i1_c3 = c(15, 20), i1_c4 = c(5, 10),
    i2_c1 = c(8, 12), i2_c2 = c(18, 22), i2_c3 = c(22, 18), i2_c4 = c(12, 8)
  )
  formula <- bmf(dprime ~ 1, criterion ~ 1, crla ~ 1, crha ~ 1,
                 crlb ~ 1, crhb ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
  expect_true(grepl("sdai_gumbel", code))
  expect_true(grepl("sdai_pfun", code))
})

test_that("sdai Gaussian produces valid stancode", {
  resp_cols <- c("i1_c1", "i1_c2", "i1_c3", "i1_c4",
                 "i2_c1", "i2_c2", "i2_c3", "i2_c4")
  model <- sdai(resp = resp_cols, n_sources = 2, n_confidence = 4,
                dist = "normal")
  dat <- data.frame(
    i1_c1 = c(10, 5), i1_c2 = c(20, 15), i1_c3 = c(15, 20), i1_c4 = c(5, 10),
    i2_c1 = c(8, 12), i2_c2 = c(18, 22), i2_c3 = c(22, 18), i2_c4 = c(12, 8)
  )
  formula <- bmf(dprime ~ 1, criterion ~ 1, sdratio ~ 1,
                 crla ~ 1, crha ~ 1, crlb ~ 1, crhb ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
  expect_true(grepl("sdai_gaussian", code))
  expect_true(grepl("integrate_1d", code))
})

test_that("sdai with n_confidence=2 produces valid stancode", {
  resp_cols <- c("i1_c1", "i1_c2", "i2_c1", "i2_c2")
  model <- sdai(resp = resp_cols, n_sources = 2, n_confidence = 2)
  dat <- data.frame(
    i1_c1 = c(30, 20), i1_c2 = c(20, 30),
    i2_c1 = c(25, 35), i2_c2 = c(35, 25)
  )
  formula <- bmf(dprime ~ 1, criterion ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
  expect_true(grepl("sdai_gumbel", code))
})


############################################################################# !
# CLOSURE FACTORY TESTS                                                  ####
############################################################################# !

test_that("sdai_make_log_lik creates a function", {
  ll_fn <- .sdai_make_log_lik(4, 2, 2, c("a", "b"), "gumbel")
  expect_true(is.function(ll_fn))
  expect_equal(names(formals(ll_fn)), c("i", "prep"))
})

test_that("sdai_make_posterior_predict creates a function", {
  pp_fn <- .sdai_make_posterior_predict(4, 2, 2, c("a", "b"), "gumbel")
  expect_true(is.function(pp_fn))
  expect_equal(names(formals(pp_fn)), c("i", "prep", "..."))
})
