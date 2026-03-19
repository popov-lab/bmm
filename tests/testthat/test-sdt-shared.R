# Test SDT shared infrastructure

############################################################################# !
# CDF HELPER TESTS                                                       ####
############################################################################# !

test_that(".sdt_cdf returns correct values for known inputs", {
  # Normal: Phi(0) = 0.5
  expect_equal(bmm:::.sdt_cdf(0, "normal"), 0.5)
  # Logistic: inv_logit(0) = 0.5
  expect_equal(bmm:::.sdt_cdf(0, "logistic"), 0.5)
  # Gumbel_min at 0: exp(-exp(0)) = exp(-1)
  expect_equal(bmm:::.sdt_cdf(0, "gumbel_min"), exp(-1), tolerance = 1e-10)
  # Gumbel_max at 0: 1 - exp(-exp(0)) = 1 - exp(-1)
  expect_equal(bmm:::.sdt_cdf(0, "gumbel_max"), 1 - exp(-1), tolerance = 1e-10)
})

test_that(".sdt_cdf is vectorized over eta", {
  eta <- c(-1, 0, 1)
  result <- bmm:::.sdt_cdf(eta, "normal")
  expect_length(result, 3)
  expect_equal(result, pnorm(eta))
})


############################################################################# !
# DPRIME AND CRITERION TESTS                                              ####
############################################################################# !

test_that("sdt_dprime computes correct values for normal distribution", {
  hr <- pnorm(1)
  far <- pnorm(-1)
  expect_equal(sdt_dprime(hr, far, dist = "normal"), 2, tolerance = 1e-10)
})

test_that("sdt_criterion computes correct values for normal distribution", {
  hr <- pnorm(1)
  far <- pnorm(-1)
  expect_equal(sdt_criterion(hr, far, dist = "normal"), 0, tolerance = 1e-10)
})

test_that("sdt_dprime validates input", {
  expect_error(sdt_dprime(0, 0.5), "between 0 and 1")
  expect_error(sdt_dprime(1, 0.5), "between 0 and 1")
  expect_error(sdt_dprime(0.5, 0), "between 0 and 1")
})


############################################################################# !
# COMBINE_SDT_RESPONSE TESTS                                              ####
############################################################################# !

test_that("combine_sdt_response maps categories correctly", {
  stim <- c(0, 0, 0, 1, 1, 1)
  conf <- c(3, 2, 1, 1, 2, 3)
  result <- combine_sdt_response(stim, conf, n_levels = 3)
  expect_equal(result, c(1, 2, 3, 4, 5, 6))
})

test_that("combine_sdt_response validates inputs", {
  expect_error(combine_sdt_response(c(0, 2), c(1, 1), n_levels = 3),
               "stimulus must be 0")
  expect_error(combine_sdt_response(c(0, 1), c(0, 1), n_levels = 3),
               "confidence must be between")
})

test_that("combine_sdt_response handles K=4 (2 confidence levels)", {
  stim <- c(0, 0, 1, 1)
  conf <- c(2, 1, 1, 2)
  result <- combine_sdt_response(stim, conf, n_levels = 2)
  expect_equal(result, c(1, 2, 3, 4))
})


############################################################################# !
# THRESHOLD CONSTRUCTION TESTS                                            ####
############################################################################# !

test_that(".sdt_make_thresholds parsimonious computes logit-spaced thresholds", {
  thr <- .sdt_make_thresholds(criterion = 0, n_ratings = 4,
                              threshold_type = "parsimonious", spacing = 0)
  expected <- log(c(1, 2, 3) / c(3, 2, 1))
  expect_equal(thr, expected, tolerance = 1e-8)

  thr_shifted <- .sdt_make_thresholds(criterion = 1, n_ratings = 4,
                                      threshold_type = "parsimonious", spacing = 0)
  expect_equal(thr_shifted, expected + 1, tolerance = 1e-8)

  thr_pos <- .sdt_make_thresholds(criterion = 0, n_ratings = 6,
                                  threshold_type = "parsimonious", spacing = 0.5)
  expect_true(all(diff(thr_pos) > 0))
})

test_that(".sdt_make_thresholds log_ratio computes ratio-scaled thresholds", {
  thr <- .sdt_make_thresholds(criterion = 0, n_ratings = 4,
                              threshold_type = "log_ratio", deltas = c(0, 0))
  expect_equal(thr, c(-1, 0, 1), tolerance = 1e-8)

  thr6 <- .sdt_make_thresholds(criterion = 0, n_ratings = 6,
                               threshold_type = "log_ratio",
                               deltas = c(0, 0, 0, 0))
  expect_equal(thr6, c(-2, -1, 0, 1, 2), tolerance = 1e-8)

  thr_shifted <- .sdt_make_thresholds(criterion = 1, n_ratings = 4,
                                      threshold_type = "log_ratio",
                                      deltas = c(0, 0))
  expect_equal(thr_shifted, c(0, 1, 2), tolerance = 1e-8)

  thr_asym <- .sdt_make_thresholds(criterion = 0, n_ratings = 6,
                                   threshold_type = "log_ratio",
                                   deltas = c(0.3, -0.2, 0.5, 0.1))
  expect_true(all(diff(thr_asym) > 0))

  thr_k4 <- .sdt_make_thresholds(criterion = 0, n_ratings = 4,
                                 threshold_type = "log_ratio",
                                 deltas = c(log(2), 0))
  expect_equal(thr_k4, c(-2, 0, 1), tolerance = 1e-8)
})


############################################################################# !
# VECTORIZATION TESTS                                                      ####
############################################################################# !

test_that(".sdt_eta is vectorized over all arguments", {
  eta <- bmm:::.sdt_eta(1.5, 0.2, c(0, 1))
  expect_length(eta, 2)
  expect_equal(eta[1], -1.5/2 - 0.2)
  expect_equal(eta[2], 1.5/2 - 0.2)

  eta <- bmm:::.sdt_eta(c(1.0, 2.0), 0.2, c(0, 1))
  expect_length(eta, 2)
  expect_equal(eta[1], -1.0/2 - 0.2)
  expect_equal(eta[2], 2.0/2 - 0.2)

  eta <- bmm:::.sdt_eta(1.5, 0.2, c(0, 1), sdratio = c(1.0, 1.3))
  expect_length(eta, 2)
  expect_equal(eta[1], -1.5/2 - 0.2)
  expect_equal(eta[2], (1.5/2 - 0.2) / 1.3)
})

test_that(".sdt_eta EV and UV give same results for noise trials", {
  eta_ev <- bmm:::.sdt_eta(1.5, 0.2, 0, sdratio = 1)
  eta_uv <- bmm:::.sdt_eta(1.5, 0.2, 0, sdratio = 1.5)
  expect_equal(eta_ev, eta_uv)
})
