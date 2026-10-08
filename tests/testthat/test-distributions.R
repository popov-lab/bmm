test_that("sdm distribution functions run without errors", {
  n <- 10
  res <- dsdm(runif(n, -pi, pi), mu = 1, c = 3, kappa = 1:n)
  expect_true(length(res) == n)
  res <- psdm(runif(n, -pi, pi), mu = rnorm(n), c = 0:(n - 1), kappa = 0:(n - 1))
  expect_true(length(res) == n)
  res <- rsdm(n, mu = rnorm(n), c = 0:(n - 1), kappa = 0:(n - 1))
  expect_true(length(res) == n)

  x <- runif(n, -pi, pi)
  res <- dsdm(x, mu = 1, c = 3, kappa = 1:n)
  res_log <- dsdm(x, mu = 1, c = 3, kappa = 1:n, log = TRUE)
  expect_true(all.equal(res, exp(res_log)))
})

test_that("dsdm integrates to 1", {
  expect_equal(integrate(dsdm, -pi, pi, mu = 0, c = 3, kappa = 3)$value, 1, tolerance = 1e-6)
})


test_that("psdm is between 0 and 1", {
  res <- psdm(runif(1000, -pi, pi), mu = runif(1000, -pi, pi), c = 3, kappa = 3)
  expect_true(all(res >= 0) && all(res <= 1))
})

test_that("psdm returns 0 for q == -pi, 0.5 for q = mu, and 1 for q almost pi, when mu == 0", {
  expect_equal(psdm(-pi, mu = 0, c = 3, kappa = 3), 0)
  expect_equal(psdm(0, mu = 0, c = 3, kappa = 3), 0.5)
  expect_equal(psdm(pi - 0.00000000001, mu = 0, c = 3, kappa = 3), 1)
})


test_that("rsdm returns values between -pi and pi", {
  res <- rsdm(1000, mu = 0, c = 3, kappa = 3)
  expect_true(all(res >= -pi) && all(res <= pi))

  res <- rsdm(1000, mu = 0, c = 3, kappa = 3)
  expect_true(all(res >= -pi) && all(res <= pi))
})

test_that("rsdm draws each value from its own parameter values", {
  withr::local_seed(445)
  n <- 20000
  group <- sample(rep(1:2, n / 2))
  pars <- data.frame(mu = c(0, 2), c = c(1, 10), kappa = c(2, 30))
  y <- rsdm(n, pars$mu[group], pars$c[group], pars$kappa[group])

  for (g in 1:2) {
    err <- y[group == g] - pars$mu[g]
    ref_cos <- integrate(function(x) cos(x) * dsdm(x, 0, pars$c[g], pars$kappa[g]), -pi, pi)$value
    expect_lt(abs(mean(cos(err)) - ref_cos), 5 * sd(cos(err)) / sqrt(length(err)))
    expect_lt(abs(mean(sin(err))), 5 * sd(sin(err)) / sqrt(length(err)))
  }
})

test_that("rsdm recycles parameters to n", {
  draw <- function(...) withr::with_seed(445, rsdm(10, ...))
  expect_identical(
    draw(mu = c(0, 2), c = c(2, 10), kappa = c(3, 30)),
    draw(mu = rep_len(c(0, 2), 10), c = rep_len(c(2, 10), 10), kappa = rep_len(c(3, 30), 10))
  )
})

test_that("rsdm returns finite draws when the density peak overflows", {
  withr::local_seed(445)
  y <- rsdm(10, c = 100, kappa = 400)
  expect_true(all(is.finite(y)))
  expect_lt(max(abs(y)), 0.05)
})

test_that("rejection_sampling takes arguments of length n per draw", {
  withr::local_seed(445)
  n <- 20000
  group <- sample(rep(1:2, n / 2))
  shape <- c(1, 9)
  # Beta(shape, 1): density shape * x^(shape - 1) peaks at shape, mean shape / (shape + 1)
  y <- rejection_sampling(
    n,
    f = function(x, group, shape) shape[group] * x^(shape[group] - 1),
    max_f = shape[group],
    proposal_fun = stats::runif,
    group = group, shape = shape
  )

  for (g in 1:2) {
    yg <- y[group == g]
    expect_lt(abs(mean(yg) - shape[g] / (shape[g] + 1)), 5 * sd(yg) / sqrt(length(yg)))
  }

  # with n = 1 every argument belongs to the single draw and is passed whole
  expect_length(rejection_sampling(1, function(x, g) g(x), 1, stats::runif, g = stats::dunif), 1)
})

test_that("rejection_sampling validates n and max_f", {
  expect_error(rejection_sampling(2.5, stats::dunif, 1, stats::runif), "whole number")
  expect_error(rejection_sampling(5, stats::dunif, 0, stats::runif), "max_f")
  expect_error(rejection_sampling(5, stats::dunif, c(1, 2), stats::runif), "max_f")
})

test_that("rejection_sampling errors instead of looping forever", {
  # a regression here hangs; fail the test instead of timing out the CI job
  setTimeLimit(elapsed = 20, transient = TRUE)
  withr::defer(setTimeLimit(elapsed = Inf))
  expect_error(rejection_sampling(5, stats::dunif, Inf, stats::runif), "max_f")
  expect_error(rsdm(5, mu = NA), "NA")
  expect_error(rejection_sampling(5, function(x) 0 * x, 1, stats::runif), "accepted")
  expect_error(rejection_sampling(5, function(x) rep(1, length(x)), 1, function(n) rep(NA_real_, n)), "NA")
})

test_that("conversion between sdm parametrizations works", {
  kappa <- rnorm(100, 5, 1)
  c_b <- rnorm(100, 5, 1)
  c_se <- c_bessel2sqrtexp(c_b, kappa)
  c_b2 <- c_sqrtexp2bessel(c_se, kappa)
  expect_equal(c_b, c_b2)
})

test_that("dsdm parametrization conversion returns accurate results", {
  y <- seq(-pi, pi, length.out = 100)
  kappa <- rnorm(100, 5, 1)
  c_b <- rnorm(100, 5, 1)
  c_se <- c_bessel2sqrtexp(c_b, kappa)
  d1 <- dsdm(y, 0, c_b, kappa, parametrization = "bessel")
  d2 <- dsdm(y, 0, c_se, kappa, parametrization = "sqrtexp")
  expect_equal(d1, d2)
})

test_that("dmixture2p integrates to 1", {
  expect_equal(integrate(dmixture2p, -pi, pi,
    mu = runif(1, min = -pi, pi),
    kappa = runif(1, min = 1, max = 20),
    p_mem = runif(1, min = 0, max = 1)
  )$value, 1, tolerance = 1e-6)
})

test_that("dmixture3p integrates to 1", {
  expect_equal(integrate(dmixture3p, -pi, pi,
    mu = runif(3, min = -pi, pi),
    kappa = runif(1, min = 1, max = 20),
    p_mem = runif(1, min = 0, max = 0.6),
    p_nt = runif(1, min = 0, max = 0.3)
  )$value, 1, tolerance = 1e-6)
})

test_that("dimm integrates to 1", {
  expect_equal(integrate(dimm, -pi, pi,
    mu = runif(3, min = -pi, pi),
    dist = c(0, runif(2, min = 0.1, max = pi)),
    kappa = runif(1, min = 1, max = 20),
    c = runif(1, min = 0, max = 3),
    a = runif(1, min = 0, max = 1),
    s = runif(1, min = 1, max = 20),
    b = 0
  )$value, 1, tolerance = 1e-6)
})

test_that("rmixture2p returns values between -pi and pi", {
  res <- rmixture2p(500,
    mu = runif(1, min = -pi, pi),
    kappa = runif(1, min = 1, max = 20),
    p_mem = runif(1, min = 0, max = 1)
  )
  expect_true(all(res >= -pi) && all(res <= pi))
})

test_that("rmixture2p draws each value from its own parameter values", {
  withr::local_seed(445)
  n <- 20000
  group <- sample(rep(1:2, n / 2))
  pars <- data.frame(mu = c(0, 2), kappa = c(2, 30), p_mem = c(0.5, 0.9))
  y <- rmixture2p(n, pars$mu[group], pars$kappa[group], pars$p_mem[group])

  for (g in 1:2) {
    err <- y[group == g] - pars$mu[g]
    ref_cos <- pars$p_mem[g] * besselI(pars$kappa[g], 1, TRUE) / besselI(pars$kappa[g], 0, TRUE)
    expect_lt(abs(mean(cos(err)) - ref_cos), 5 * sd(cos(err)) / sqrt(length(err)))
    expect_lt(abs(mean(sin(err))), 5 * sd(sin(err)) / sqrt(length(err)))
  }
})

test_that("rmixture2p recycles parameters to n", {
  draw <- function(...) withr::with_seed(445, rmixture2p(10, ...))
  expect_identical(
    draw(mu = c(0, 2), kappa = c(2, 30), p_mem = c(0.5, 0.9)),
    draw(mu = rep_len(c(0, 2), 10), kappa = rep_len(c(2, 30), 10), p_mem = rep_len(c(0.5, 0.9), 10))
  )
})

test_that("rmixture3p returns values between -pi and pi", {
  res <- rmixture3p(500,
    mu = runif(3, min = -pi, pi),
    kappa = runif(1, min = 1, max = 20),
    p_mem = runif(1, min = 0, max = 0.6),
    p_nt = runif(1, min = 0, max = 0.3)
  )
  expect_true(all(res >= -pi) && all(res <= pi))
})

test_that("rimm returns values between -pi and pi", {
  res <- rimm(500,
    mu = runif(3, min = -pi, pi),
    dist = c(0, runif(2, min = 0.1, max = pi)),
    kappa = runif(1, min = 1, max = 20),
    c = runif(1, min = 0, max = 3),
    a = runif(1, min = 0, max = 1),
    s = runif(1, min = 1, max = 20),
    b = 0
  )
  expect_true(all(res >= -pi) && all(res <= pi))
})

test_that("dm3 requires custom act_funs to be specified", {
  model <- m3(
    resp_cats = c("corr", "other", "dist", "npl"),
    num_options = c("n_corr", "n_other", "n_dist", "n_npl"),
    choice_rule = "simple",
    version = "custom"
  )
  expect_error(
    dm3(x = c(10, 10, 10, 10), pars = c(a = 1, b = 1, c = 1, f = 1), m3_model = model),
    "Activation functions can only be generated"
  )
})

test_that("dm3 works for a simple m3 model", {
  model <- m3(
    resp_cats = c("corr", "other", "npl"),
    num_options = c(1, 4, 5),
    choice_rule = "simple",
    version = "ss"
  )
  dens <- dm3(x = c(20, 10, 10), pars = c(a = 1, b = 1, c = 2), m3_model = model)
  expect_type(dens, "double")
  expect_length(dens, 1)
  # compare with # compare with lgamma(size + 1) + sum(x * log(prob) - lgamma(x + 1))
  expect_equal(
    dens,
    lgamma(41) - lgamma(21) - lgamma(11) - lgamma(11) +
      sum(log(c(1 + 1 + 2, (1 + 1) * 4, 1 * 5) / sum(c(1 + 1 + 2, (1 + 1) * 4, 1 * 5))) * c(20, 10, 10))
  )
})

test_that("dm3 works for a complex span m3 model", {
  model <- m3(
    resp_cats = c("corr", "dist_context", "other", "dist_other", "npl"),
    num_options = c(1, 10, 4, 10, 5),
    choice_rule = "simple",
    version = "cs"
  )
  dens <- dm3(x = c(20, 5, 10, 5, 10), pars = c(a = 1, b = 1, c = 2, f = 0), m3_model = model)
  expect_type(dens, "double")
  expect_length(dens, 1)
  # compare with lgamma(size + 1) + sum(x * log(prob) - lgamma(x + 1))
  expect_equal(
    dens,
    lgamma(51) - lgamma(21) - 2 * lgamma(11) - 2 * lgamma(6) +
      sum(
        log(
          c(1 + 1 + 2, 1 * 10, (1 + 1) * 4, 1 * 10, 1 * 5)
          / sum(c(1 + 1 + 2, 1 * 10, (1 + 1) * 4, 1 * 10, 1 * 5))
        )
        * c(20, 5, 10, 5, 10)
      )
  )
})

test_that("dm3 works for a custom m3 model", {
  model <- m3(
    resp_cats = c("correct", "lures", "nonpresented"),
    num_options = c(1, 4, 5),
    choice_rule = "simple",
    version = "custom"
  )
  act_funs <- bmf(
    correct ~ background + item + binding,
    lures ~ background + item,
    nonpresented ~ background
  )
  dens <- dm3(
    x = c(20, 10, 10), pars = c(background = 1, item = 1, binding = 2),
    m3_model = model, act_funs = act_funs
  )
  expect_type(dens, "double")
  expect_length(dens, 1)
  # compare with lgamma(size + 1) + sum(x * log(prob) - lgamma(x + 1))
  expect_equal(
    dens,
    lgamma(41) - lgamma(21) - lgamma(11) - lgamma(11) +
      sum(log(c(1 + 1 + 2, (1 + 1) * 4, 1 * 5) / sum(c(1 + 1 + 2, (1 + 1) * 4, 1 * 5))) * c(20, 10, 10))
  )
})

test_that("rm3 works for a simple m3 model", {
  model <- m3(
    resp_cats = c("corr", "other", "npl"),
    num_options = c(1, 4, 50),
    choice_rule = "simple",
    version = "ss"
  )
  res <- rm3(n = 10, size = 100, pars = c(a = 1, b = 1, c = 2), m3_model = model)
  expect_type(res, "integer")
  expect_true("matrix" %in% class(res))
  expect_true(nrow(res) == 10 && ncol(res) == 3)
  expect_true(all(rowSums(res) == 100))
  expect_equal(colnames(res), model$resp_vars$resp_cats)
  expect_true(median(res[, "npl"]) > median(res))
})

test_that("rm3 works for a complexspan m3 model", {
  model <- m3(
    resp_cats = c("corr", "dist_context", "other", "dist_other", "npl"),
    num_options = c(1, 10, 4, 100, 5),
    choice_rule = "simple",
    version = "cs"
  )
  res <- rm3(n = 10, size = 100, pars = c(a = 1, b = 1, c = 2, f = 0), m3_model = model)
  expect_type(res, "integer")
  expect_true("matrix" %in% class(res))
  expect_true(nrow(res) == 10 && ncol(res) == 5)
  expect_true(all(rowSums(res) == 100))
  expect_equal(colnames(res), model$resp_vars$resp_cats)
  expect_true(median(res[, "dist_other"]) > median(res))
})

# DDM distribution tests ----
test_that("dddm runs without errors with scalar inputs", {
  res <- dddm(0.5, 1, drift = 2, bound = 1.5, ndt = 0.3)
  expect_type(res, "double")
  expect_length(res, 1)
  expect_false(is.na(res))
  expect_false(is.infinite(res))
})

test_that("dddm handles numeric response codes correctly", {
  res_upper <- dddm(0.5, 1, drift = 2, bound = 1.5, ndt = 0.3)
  res_lower <- dddm(0.5, 0, drift = 2, bound = 1.5, ndt = 0.3)
  expect_type(res_upper, "double")
  expect_type(res_lower, "double")
  expect_false(identical(res_upper, res_lower))
})

test_that("dddm handles character response codes correctly", {
  res_numeric <- dddm(0.5, 1, drift = 2, bound = 1.5, ndt = 0.3)
  res_char <- dddm(0.5, "upper", drift = 2, bound = 1.5, ndt = 0.3)
  expect_equal(res_numeric, res_char)
})

test_that("dddm log parameter works correctly", {
  res <- dddm(0.5, 1, drift = 2, bound = 1.5, ndt = 0.3, log = FALSE)
  res_log <- dddm(0.5, 1, drift = 2, bound = 1.5, ndt = 0.3, log = TRUE)
  expect_equal(log(res), res_log)
})

test_that("dddm vectorizes correctly with scalar rt/response and vector parameters", {
  # This is the log_lik case: single observation, multiple posterior samples
  rt_single <- 0.5
  resp_single <- 1
  drift_vec <- rnorm(10, 2, 0.3)
  bound_vec <- rep(1.5, 10)
  ndt_vec <- rep(0.3, 10)
  
  res <- dddm(rt_single, resp_single, drift_vec, bound_vec, ndt_vec)
  expect_length(res, 10)
  expect_false(any(is.na(res)))
  expect_false(any(is.infinite(res)))
})

test_that("dddm vectorizes correctly with vector rt/response and scalar parameters", {
  rt_vec <- c(0.4, 0.5, 0.6)
  resp_vec <- c(1, 1, 0)
  drift_scalar <- 2.0
  
  res <- dddm(rt_vec, resp_vec, drift_scalar, bound = 1.5, ndt = 0.3)
  expect_length(res, 3)
  expect_false(any(is.na(res)))
})

test_that("dddm vectorizes correctly with vector rt/response and vector parameters", {
  n <- 5
  rt_vec <- runif(n, 0.3, 0.8)
  resp_vec <- sample(0:1, n, replace = TRUE)
  drift_vec <- rnorm(n, 2, 0.3)
  
  res <- dddm(rt_vec, resp_vec, drift_vec, bound = 1.5, ndt = 0.3)
  expect_length(res, n)
  expect_false(any(is.na(res)))
})

test_that("dddm rejects negative RTs", {
  expect_error(
    dddm(-0.5, 1, drift = 2, bound = 1.5, ndt = 0.3),
    "Negative RTs are not allowed"
  )
})

test_that("dddm rejects invalid response codes", {
  expect_error(
    dddm(0.5, 2, drift = 2, bound = 1.5, ndt = 0.3),
    "Invalid numeric responses"
  )
  expect_error(
    dddm(0.5, "invalid", drift = 2, bound = 1.5, ndt = 0.3),
    "Invalid responses"
  )
})

test_that("dddm rejects mismatched rt and response lengths", {
  expect_error(
    dddm(c(0.5, 0.6), 1, drift = 2, bound = 1.5, ndt = 0.3),
    "Different number of RTs and responses"
  )
})

test_that("rddm returns data frame with correct structure", {
  res <- rddm(10, drift = 2, bound = 1.5, ndt = 0.3)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 10)
  expect_true(all(c("rt", "response") %in% names(res)))
  expect_true(all(res$rt > 0))
  expect_true(all(res$response %in% c(0, 1)))
})

test_that("rddm handles vector parameters correctly", {
  n <- 5
  drift_vec <- rnorm(n, 2, 0.3)
  res <- rddm(n, drift = drift_vec, bound = 1.5, ndt = 0.3)
  expect_equal(nrow(res), n)
})

test_that("rddm samples from correct boundaries based on drift", {
  # Positive drift should favor upper boundary (response = 1)
  res_pos <- rddm(1000, drift = 3, bound = 1.5, ndt = 0.3)
  expect_gt(mean(res_pos$response == 1), 0.8)
  
  # Negative drift should favor lower boundary (response = 0)
  res_neg <- rddm(1000, drift = -3, bound = 1.5, ndt = 0.3)
  expect_gt(mean(res_neg$response == 0), 0.8)
})

test_that("rddm respects non-decision time", {
  res <- rddm(100, drift = 2, bound = 1.5, ndt = 0.3)
  expect_true(all(res$rt >= 0.3))
})

test_that("rddm starting point bias affects response proportions", {
  # Bias toward upper boundary
  res_upper <- rddm(1000, drift = 0, bound = 1.5, ndt = 0.3, zr = 0.7)
  expect_gt(mean(res_upper$response == 1), 0.6)
  
  # Bias toward lower boundary
  res_lower <- rddm(1000, drift = 0, bound = 1.5, ndt = 0.3, zr = 0.3)
  expect_gt(mean(res_lower$response == 0), 0.6)
})

# -----------------------------------------------------------------------------
# cswald distribution function tests
# -----------------------------------------------------------------------------

test_that("dcswald runs without errors and returns correct output", {
  n <- 10
  rt <- runif(n, 0.4, 1.5)
  response <- sample(c(0, 1), n, replace = TRUE)

  # simple version
  res <- dcswald(rt, response,
    drift = 2, bound = 1.5, ndt = 0.3,
    version = "simple"
  )
  expect_length(res, n)
  expect_type(res, "double")

  # crisk version
  res_crisk <- dcswald(rt, response,
    drift = 2, bound = 1.5, ndt = 0.3,
    zr = 0.5, version = "crisk"
  )
  expect_length(res_crisk, n)
  expect_type(res_crisk, "double")

  # test log = FALSE
  res_exp <- dcswald(rt, response,
    drift = 2, bound = 1.5, ndt = 0.3,
    version = "simple", log = FALSE
  )
  expect_true(all(res_exp >= 0))
  expect_equal(res_exp, exp(res))
})

test_that("dcswald handles vectorized parameters", {
  n <- 5
  rt <- runif(n, 0.4, 1.5)
  response <- sample(c(0, 1), n, replace = TRUE)

  # vectorized drift
  res <- dcswald(rt, response,
    drift = 1:n, bound = 1.5, ndt = 0.3,
    version = "simple"
  )
  expect_length(res, n)

  # all parameters vectorized
  res <- dcswald(rt, response,
    drift = 1:n, bound = seq(1, 2, length.out = n),
    ndt = seq(0.1, 0.3, length.out = n), version = "simple"
  )
  expect_length(res, n)
})

test_that("rcswald generates valid output", {
  n <- 100
  out <- rcswald(n, drift = 2, bound = 1.5, ndt = 0.3)

  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), n)
  expect_true(all(c("rt", "response") %in% names(out)))
  expect_true(all(out$rt > 0))
  expect_true(all(out$response %in% c(0, 1)))
})

test_that("rcswald handles vectorized parameters", {
  n <- 10
  out <- rcswald(n,
    drift = runif(n, 1, 3), bound = runif(n, 1, 2),
    ndt = runif(n, 0.1, 0.3), zr = runif(n, 0.3, 0.7)
  )
  expect_equal(nrow(out), n)
})

test_that("pcswald returns probabilities between 0 and 1 for simple version", {
  q <- seq(0.4, 2, length.out = 20)
  p <- pcswald(q,
    response = 1, drift = 2, bound = 1.5, ndt = 0.3,
    version = "simple"
  )

  expect_length(p, 20)
  expect_true(all(p >= 0, na.rm = TRUE))
  expect_true(all(p <= 1, na.rm = TRUE))
  # CDF should be monotonically increasing
  expect_true(all(diff(p) >= 0))
})

test_that("pcswald returns probabilities between 0 and 1 for crisk version", {
  q <- seq(0.4, 2, length.out = 20)

  # upper boundary
  p_upper <- pcswald(q,
    response = 1, drift = 2, bound = 1.5, ndt = 0.3,
    version = "crisk"
  )
  expect_true(all(p_upper >= 0))
  expect_true(all(p_upper <= 1))

  # lower boundary
  p_lower <- pcswald(q,
    response = 0, drift = 2, bound = 1.5, ndt = 0.3,
    version = "crisk"
  )
  expect_true(all(p_lower >= 0))
  expect_true(all(p_lower <= 1))
})

test_that("pcswald warns for response=0 in simple version", {
  expect_warning(
    pcswald(1.0,
      response = 0, drift = 2, bound = 1.5, ndt = 0.3,
      version = "simple"
    ),
    "CDF for response=0 is not well-defined"
  )
})

test_that("pcswald handles log.p and lower.tail arguments", {
  q <- 1.0
  p <- pcswald(q, response = 1, drift = 2, bound = 1.5, ndt = 0.3)

  # log.p = TRUE
  p_log <- pcswald(q,
    response = 1, drift = 2, bound = 1.5, ndt = 0.3,
    log.p = TRUE
  )
  expect_equal(p_log, log(p))

  # lower.tail = FALSE
  p_upper <- pcswald(q,
    response = 1, drift = 2, bound = 1.5, ndt = 0.3,
    lower.tail = FALSE
  )
  expect_equal(p_upper, 1 - p)
})

test_that("qcswald returns valid quantiles for simple version", {
  p <- c(0.1, 0.5, 0.9)
  q <- qcswald(p,
    response = 1, drift = 2, bound = 1.5, ndt = 0.3,
    version = "simple"
  )

  expect_length(q, 3)
  expect_true(all(q > 0.3)) # all quantiles > ndt
  expect_true(all(diff(q) > 0)) # monotonically increasing
})

test_that("qcswald returns valid quantiles for crisk version", {
  p <- c(0.1, 0.5, 0.9)
  q <- qcswald(p,
    response = 1, drift = 2, bound = 1.5, ndt = 0.3,
    version = "crisk"
  )

  expect_length(q, 3)
  expect_true(all(q > 0.3)) # all quantiles > ndt
})

test_that("qcswald warns for response=0 in simple version", {
  expect_warning(
    qcswald(0.5,
      response = 0, drift = 2, bound = 1.5, ndt = 0.3,
      version = "simple"
    ),
    "Quantile for response=0 is not well-defined"
  )
})

test_that("pcswald and qcswald are inverses (round-trip)", {
  p_orig <- c(0.1, 0.3, 0.5, 0.7, 0.9)

  # simple version
  q <- qcswald(p_orig,
    response = 1, drift = 2, bound = 1.5, ndt = 0.3,
    version = "simple"
  )
  p_back <- pcswald(q,
    response = 1, drift = 2, bound = 1.5, ndt = 0.3,
    version = "simple"
  )
  expect_equal(p_orig, p_back, tolerance = 1e-4)

  # crisk version
  q_crisk <- qcswald(p_orig,
    response = 1, drift = 2, bound = 1.5, ndt = 0.3,
    version = "crisk"
  )
  p_back_crisk <- pcswald(q_crisk,
    response = 1, drift = 2, bound = 1.5,
    ndt = 0.3, version = "crisk"
  )
  expect_equal(p_orig, p_back_crisk, tolerance = 1e-4)
})

test_that("cswald parameter validation works", {
  # negative bound
  expect_error(
    dcswald(1, 1, drift = 2, bound = -1, ndt = 0.3),
    "boundary"
  )

  # negative ndt
  expect_error(
    dcswald(1, 1, drift = 2, bound = 1.5, ndt = -0.1),
    "non-decision time"
  )

  expect_error(
    dcswald(1, 1, drift = 2, bound = 1.5, ndt = 0.3, zr = 1.5),
    "starting point"
  )

  # negative s
  expect_error(
    dcswald(1, 1, drift = 2, bound = 1.5, ndt = 0.3, s = -1),
    "diffusion constant"
  )
})

test_that("dcswald errors when rt < ndt", {
  expect_error(
    dcswald(rt = 0.2, response = 1, drift = 2, bound = 1.5, ndt = 0.3),
    "smaller than the non-decision time"
  )
})

test_that("cswald simple and crisk give similar results for correct responses", {
  # For high drift and unbiased starting point (zr=0.5), the simple version

  # should give similar densities to crisk for correct responses (response=1)
  rt <- seq(0.4, 1.5, by = 0.1)
  drift <- 3
  bound_simple <- 1.5 # simple version: distance to correct boundary
  bound_crisk <- 3.0 # crisk version: total boundary separation (2x simple)

  # Get densities for correct responses
  dens_simple <- dcswald(rt,
    response = 1, drift = drift, bound = bound_simple,
    ndt = 0.2, version = "simple", log = FALSE
  )
  dens_crisk <- dcswald(rt,
    response = 1, drift = drift, bound = bound_crisk,
    ndt = 0.2, zr = 0.5, version = "crisk", log = FALSE
  )

  # Densities should be correlated (similar shape) though not identical
  # because crisk accounts for competing accumulator
  cor_value <- cor(dens_simple, dens_crisk)
  expect_true(cor_value > 0.95)
})

test_that("qcswald adaptive bounds work for slow drift", {
  # With very slow drift, RTs can be very long
  # The adaptive bounds should handle this
  q <- qcswald(p = 0.5, response = 1, drift = 0.1, bound = 2, ndt = 0.3)
  expect_true(is.finite(q))
  expect_true(q > 0.3) # greater than ndt

  # Verify round-trip
  p_back <- pcswald(q, response = 1, drift = 0.1, bound = 2, ndt = 0.3)
  expect_equal(p_back, 0.5, tolerance = 0.01)
})
test_that("dcswald handles extreme drift values", {
  # Very high drift - fast accurate responses
  expect_no_error(
    dcswald(rt = 0.5, response = 1, drift = 10, bound = 1.5, ndt = 0.2)
  )

  # Very low drift - slow responses
  expect_no_error(
    dcswald(rt = 2.0, response = 1, drift = 0.1, bound = 1.5, ndt = 0.2)
  )

  # Negative drift for crisk version
  expect_no_error(
    dcswald(
      rt = 0.8, response = 0, drift = -2, bound = 1.5, ndt = 0.2,
      zr = 0.5, version = "crisk"
    )
  )
})

test_that("dcswald handles extreme bound values", {
  # Very small bound - fast responses
  expect_no_error(
    dcswald(rt = 0.3, response = 1, drift = 2, bound = 0.5, ndt = 0.1)
  )

  # Very large bound - slow responses
  expect_no_error(
    dcswald(rt = 5.0, response = 1, drift = 2, bound = 5, ndt = 0.2)
  )
})

test_that("dcswald handles boundary zr values for crisk version", {
  # zr near lower boundary
  expect_no_error(
    dcswald(
      rt = 0.8, response = 1, drift = 2, bound = 1.5, ndt = 0.2,
      zr = 0.01, version = "crisk"
    )
  )

  # zr near upper boundary
  expect_no_error(
    dcswald(
      rt = 0.8, response = 1, drift = 2, bound = 1.5, ndt = 0.2,
      zr = 0.99, version = "crisk"
    )
  )

  # zr exactly at boundaries should error
  expect_error(
    dcswald(
      rt = 0.8, response = 1, drift = 2, bound = 1.5, ndt = 0.2,
      zr = 0, version = "crisk"
    ),
    "between 0 and 1"
  )

  expect_error(
    dcswald(
      rt = 0.8, response = 1, drift = 2, bound = 1.5, ndt = 0.2,
      zr = 1, version = "crisk"
    ),
    "between 0 and 1"
  )
})

test_that("rcswald generates valid data with extreme parameters", {
  # High drift - should produce mostly correct responses
  dat_high <- rcswald(n = 100, drift = 5, bound = 1.5, ndt = 0.2)
  expect_true(mean(dat_high$response == 1) > 0.9)
  expect_true(all(dat_high$rt > 0.2)) # all RTs > ndt

  # Low drift - should produce more errors
  dat_low <- rcswald(n = 100, drift = 0.5, bound = 1.5, ndt = 0.2)
  expect_true(mean(dat_low$response == 0) > 0.1) # some errors

  # Extreme zr values
  dat_biased_upper <- rcswald(n = 100, drift = 1, bound = 2, ndt = 0.2, zr = 0.9)
  expect_true(mean(dat_biased_upper$response == 1) > 0.7) # biased toward upper

  dat_biased_lower <- rcswald(n = 100, drift = 1, bound = 2, ndt = 0.2, zr = 0.1)
  expect_true(mean(dat_biased_lower$response == 0) > 0.3) # biased toward lower
})

test_that("pcswald and qcswald handle extreme parameters", {
  # pcswald with high drift should give CDF close to 1 quickly
  p_high <- pcswald(q = 1.0, response = 1, drift = 5, bound = 1, ndt = 0.2)
  expect_true(p_high > 0.99)

  # qcswald should return valid quantiles for correct responses
  q_50 <- qcswald(p = 0.5, response = 1, drift = 2, bound = 1.5, ndt = 0.3)
  expect_true(q_50 > 0.3) # greater than ndt

  # Very small p should give RT close to ndt
  q_small <- qcswald(p = 0.01, response = 1, drift = 2, bound = 1.5, ndt = 0.3)
  expect_true(q_small > 0.3 && q_small < 0.5)
})

test_that("dcswald error-response likelihood stays finite in the upper tail (#376)", {
  # The simple version routes response = 0 through the shifted-Wald survival.
  # The old naive log(1 - exp(cdf)) cancelled to NaN/-Inf for late RTs.
  ll <- dcswald(c(5, 10, 20, 40),
    response = 0, drift = 2, bound = 1, ndt = 0.3, version = "simple"
  )
  expect_true(all(is.finite(ll)))
})

test_that("dcswald error-response survival equals 1 - CDF in mid-range (#376)", {
  rt <- c(0.4, 0.6, 0.9, 1.3)
  ll_error <- dcswald(rt,
    response = 0, drift = 2, bound = 1, ndt = 0.3, version = "simple"
  )
  cdf <- pcswald(rt,
    response = 1, drift = 2, bound = 1, ndt = 0.3, version = "simple"
  )
  expect_equal(ll_error, log(1 - cdf))
})

test_that("rm3 works without providing b parameter", {
  model <- m3(
    resp_cats = c("corr", "other", "npl"),
    num_options = c(1, 4, 5),
    choice_rule = "simple",
    version = "ss"
  )
  # b parameter should be added automatically
  res <- rm3(n = 10, size = 100, pars = c(a = 1, c = 2), m3_model = model)
  expect_type(res, "integer")
  expect_true("matrix" %in% class(res))
  expect_true(nrow(res) == 10 && ncol(res) == 3)
  expect_true(all(rowSums(res) == 100))
  expect_equal(colnames(res), model$resp_vars$resp_cats)
})

test_that("rm3 works with full bmmformula", {
  model <- m3(
    resp_cats = c("corr", "other", "npl"),
    num_options = c(1, 4, 5),
    choice_rule = "simple",
    version = "ss"
  )
  # Full formula with both activation formulas and parameter formulas
  full_formula <- bmf(
    corr ~ b + a + c,
    other ~ b + a,
    npl ~ b,
    a ~ 1,
    c ~ 1
  )
  res <- rm3(
    n = 10, size = 100, pars = c(a = 1, c = 2),
    m3_model = model, act_funs = full_formula
  )
  expect_type(res, "integer")
  expect_true("matrix" %in% class(res))
  expect_true(nrow(res) == 10 && ncol(res) == 3)
  expect_true(all(rowSums(res) == 100))
  expect_equal(colnames(res), model$resp_vars$resp_cats)
})

test_that("dm3 works without providing b parameter", {
  model <- m3(
    resp_cats = c("corr", "other", "npl"),
    num_options = c(1, 4, 5),
    choice_rule = "simple",
    version = "ss"
  )
  # b parameter should be added automatically
  dens <- dm3(x = c(20, 10, 10), pars = c(a = 1, c = 2), m3_model = model)
  expect_type(dens, "double")
  expect_length(dens, 1)
})

test_that("dm3 works with full bmmformula", {
  model <- m3(
    resp_cats = c("corr", "other", "npl"),
    num_options = c(1, 4, 5),
    choice_rule = "simple",
    version = "ss"
  )
  full_formula <- bmf(
    corr ~ b + a + c,
    other ~ b + a,
    npl ~ b,
    a ~ 1,
    c ~ 1
  )
  dens <- dm3(
    x = c(20, 10, 10), pars = c(a = 1, c = 2),
    m3_model = model, act_funs = full_formula
  )
  expect_type(dens, "double")
  expect_length(dens, 1)
})

test_that("dm3 matches activations and num_options by category, not by position", {
  model <- m3(resp_cats = c("corr", "other", "npl"), num_options = c(1, 4, 5), choice_rule = "simple")
  pars <- c(a = 1, b = 0.1, c = 2)
  expected <- dm3(c(20, 10, 10), pars, model, bmf(corr ~ b + a + c, other ~ b + a, npl ~ b))
  expect_equal(dm3(c(20, 10, 10), pars, model, bmf(npl ~ b, other ~ b + a, corr ~ b + a + c)), expected)

  model$other_vars$num_options <- c(npl = 5, corr = 1, other = 4)
  expect_equal(dm3(c(20, 10, 10), pars, model, bmf(corr ~ b + a + c, other ~ b + a, npl ~ b)), expected)
  expect_error(dm3(c(20, 10, 10), pars, model, bmf(corr ~ b + a + c, npl ~ b)), "'other'")
})

test_that("dm3 and rm3 use Gaussian-rule probabilities for choice_rule = 'gaussian'", {
  model <- m3(resp_cats = c("corr", "other", "npl"), num_options = c(1, 4, 10),
              choice_rule = "gaussian", version = "ss")
  probs <- exp(vapply(1:3, function(k) m3_gauss_logp(k, 2.5, 1, 0, 1, 4, 10), numeric(1)))
  expect_equal(sum(probs), 1, tolerance = 1e-6)
  expect_equal(dm3(c(60, 25, 15), c(a = 1, c = 1.5), model),
               stats::dmultinom(c(60, 25, 15), prob = probs, log = TRUE), tolerance = 1e-8)
  # the softmax rule gives different probabilities for the same activations
  softmax <- m3(resp_cats = c("corr", "other", "npl"), num_options = c(1, 4, 10), version = "ss")
  expect_gt(abs(dm3(c(60, 25, 15), c(a = 1, c = 1.5), softmax) - dm3(c(60, 25, 15), c(a = 1, c = 1.5), model)), 0.1)
  withr::local_seed(1)
  draws <- rm3(n = 2000, size = 1, pars = c(a = 1, c = 1.5), m3_model = model)
  expect_lt(max(abs(colMeans(draws) - probs)), 0.04)
})

test_that("rm3 errors when full formula has no activation functions", {
  model <- m3(
    resp_cats = c("corr", "other", "npl"),
    num_options = c(1, 4, 5),
    choice_rule = "simple",
    version = "ss"
  )
  # Formula with only parameter formulas, no activation formulas
  wrong_formula <- bmf(
    a ~ 1,
    c ~ 1
  )
  expect_error(
    rm3(
      n = 10, size = 100, pars = c(a = 1, c = 2),
      m3_model = model, act_funs = wrong_formula
    ),
    "No activation formulas found"
  )
})

test_that("rm3 with unpack=TRUE returns named vector for n=1", {
  model <- m3(
    resp_cats = c("corr", "other", "npl"),
    num_options = c(1, 4, 5),
    choice_rule = "simple",
    version = "ss"
  )

  # Test unpack=TRUE
  result <- rm3(
    n = 1, size = 100, pars = c(a = 1, c = 2),
    m3_model = model, unpack = TRUE
  )

  expect_type(result, "integer")
  expect_false("matrix" %in% class(result))
  expect_true("integer" %in% class(result))
  expect_equal(names(result), c("corr", "other", "npl"))
  expect_equal(sum(result), 100)
})

test_that("rm3 with unpack=TRUE still returns matrix for n>1", {
  model <- m3(
    resp_cats = c("corr", "other", "npl"),
    num_options = c(1, 4, 5),
    choice_rule = "simple",
    version = "ss"
  )

  # unpack should be ignored when n > 1
  result <- rm3(
    n = 5, size = 100, pars = c(a = 1, c = 2),
    m3_model = model, unpack = TRUE
  )

  expect_true("matrix" %in% class(result))
  expect_equal(dim(result), c(5, 3))
  expect_equal(colnames(result), c("corr", "other", "npl"))
})


# Tests for EZDM distribution functions ----------------------------------------

test_that("dezdm 3par runs without errors and returns correct output", {
  ll <- dezdm(
    mean_rt = 0.5, var_rt = 0.02, n_upper = 80, n_trials = 100,
    drift = 2, bound = 1.5, ndt = 0.3
  )
  expect_type(ll, "double")
  expect_length(ll, 1)
  expect_true(is.finite(ll))
  expect_true(ll < 0)
})

test_that("dezdm 4par runs without errors and returns correct output", {
  ll <- dezdm(
    mean_rt = c(0.45, 0.55), var_rt = c(0.018, 0.025),
    n_upper = 80, n_trials = 100,
    drift = 2, bound = 1.5, ndt = 0.3, zr = 0.55, version = "4par"
  )
  expect_type(ll, "double")
  expect_length(ll, 1)
  expect_true(is.finite(ll))
  expect_true(ll < 0)
})

test_that("dezdm log and non-log outputs are consistent", {
  args <- list(
    mean_rt = 0.5, var_rt = 0.02, n_upper = 80, n_trials = 100,
    drift = 2, bound = 1.5, ndt = 0.3
  )
  ll_log <- do.call(dezdm, c(args, log = TRUE))
  ll_nolog <- do.call(dezdm, c(args, log = FALSE))
  expect_equal(exp(ll_log), ll_nolog)
})

test_that("rezdm 3par returns correct output structure", {
  res <- rezdm(n = 100, n_trials = 50, drift = 2, bound = 1.5, ndt = 0.3)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 100)
  expect_equal(ncol(res), 4)
  expect_true(all(c("mean_rt", "var_rt", "n_upper", "n_trials") %in% names(res)))
})

test_that("rezdm 4par returns correct output structure", {
  res <- rezdm(
    n = 100, n_trials = 50, drift = 2, bound = 1.5, ndt = 0.3,
    zr = 0.55, version = "4par"
  )
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 100)
  expect_equal(ncol(res), 6)
  expect_true(all(
    c(
      "mean_rt_upper", "mean_rt_lower", "var_rt_upper",
      "var_rt_lower", "n_upper", "n_trials"
    ) %in% names(res)
  ))
})

test_that("rezdm 3par returns plausible values", {
  res <- rezdm(n = 500, n_trials = 100, drift = 2, bound = 1.5, ndt = 0.3)

  # mean RT should be greater than ndt
  expect_true(all(res$mean_rt > 0.3))

  # variance should be positive
  expect_true(all(res$var_rt > 0))

  # n_upper should be between 0 and n_trials

  expect_true(all(res$n_upper >= 0))
  expect_true(all(res$n_upper <= 100))

  # n_trials should all be 100
  expect_true(all(res$n_trials == 100))
})

test_that("rezdm 4par returns plausible values", {
  set.seed(123)
  res <- rezdm(
    n = 500, n_trials = 100, drift = 2, bound = 1.5, ndt = 0.3,
    zr = 0.5, version = "4par"
  )

  # mean RT should be positive (where not NA)
  # Note: individual samples can occasionally be < ndt due to sampling variation
  expect_true(all(res$mean_rt_upper[!is.na(res$mean_rt_upper)] > 0))
  expect_true(all(res$mean_rt_lower[!is.na(res$mean_rt_lower)] > 0))

  # on average, mean RT should be greater than ndt
  expect_true(mean(res$mean_rt_upper, na.rm = TRUE) > 0.3)
  expect_true(mean(res$mean_rt_lower, na.rm = TRUE) > 0.3)

  # variance should be positive (where not NA)
  expect_true(all(res$var_rt_upper[!is.na(res$var_rt_upper)] > 0))
  expect_true(all(res$var_rt_lower[!is.na(res$var_rt_lower)] > 0))

  # n_upper should be between 0 and n_trials
  expect_true(all(res$n_upper >= 0))
  expect_true(all(res$n_upper <= 100))
})

test_that("dezdm validates parameters correctly", {
  # bound must be positive
  expect_error(
    dezdm(
      mean_rt = 0.5, var_rt = 0.02, n_upper = 80, n_trials = 100,
      drift = 2, bound = -1, ndt = 0.3
    ),
    "bound must be positive"
  )

  # ndt must be positive
  expect_error(
    dezdm(
      mean_rt = 0.5, var_rt = 0.02, n_upper = 80, n_trials = 100,
      drift = 2, bound = 1.5, ndt = -0.1
    ),
    "ndt must be positive"
  )

  # n_upper cannot exceed n_trials
  expect_error(
    dezdm(
      mean_rt = 0.5, var_rt = 0.02, n_upper = 150, n_trials = 100,
      drift = 2, bound = 1.5, ndt = 0.3
    ),
    "n_upper cannot exceed n_trials"
  )

  # n_trials must be larger than 2
  expect_error(
    dezdm(
      mean_rt = 0.5, var_rt = 0.02, n_upper = 1, n_trials = 1,
      drift = 2, bound = 1.5, ndt = 0.3
    ),
    "n_trials must be larger than 2"
  )

  expect_error(
    dezdm(
      mean_rt = 0.5, var_rt = 0.02, n_upper = 2, n_trials = 2,
      drift = 2, bound = 1.5, ndt = 0.3
    ),
    "n_trials must be larger than 2"
  )

  # the logit-scale binomial would return a finite density for fractional counts
  expect_error(
    dezdm(
      mean_rt = 0.5, var_rt = 0.02, n_upper = 80.5, n_trials = 100,
      drift = 2, bound = 1.5, ndt = 0.3
    ),
    "must be whole numbers"
  )
  expect_error(
    dezdm(
      mean_rt = 0.5, var_rt = 0.02, n_upper = 80, n_trials = 100.5,
      drift = 2, bound = 1.5, ndt = 0.3
    ),
    "must be whole numbers"
  )

  # version must be valid
  expect_error(
    dezdm(
      mean_rt = 0.5, var_rt = 0.02, n_upper = 80, n_trials = 100,
      drift = 2, bound = 1.5, ndt = 0.3, version = "5par"
    ),
    "should be one of"
  )

  # 4par requires length-2 vectors
  expect_error(
    dezdm(
      mean_rt = 0.5, var_rt = 0.02, n_upper = 80, n_trials = 100,
      drift = 2, bound = 1.5, ndt = 0.3, version = "4par"
    ),
    "must be length 2"
  )

  # zr must be between 0 and 1 for 4par
  expect_error(
    dezdm(
      mean_rt = c(0.5, 0.6), var_rt = c(0.02, 0.03),
      n_upper = 80, n_trials = 100,
      drift = 2, bound = 1.5, ndt = 0.3, zr = 1.5, version = "4par"
    ),
    "zr must be between 0 and 1"
  )
})

test_that("rezdm 4par handles cells with fewer than 2 responses at a boundary", {
  withr::local_seed(1)
  res <- rezdm(n = 200, n_trials = 3, drift = 0, bound = 1, ndt = 0.3,
               version = "4par")
  expect_identical(is.na(res$mean_rt_upper), res$n_upper < 2)
  expect_identical(is.na(res$var_rt_upper), res$n_upper < 2)
  expect_identical(is.na(res$mean_rt_lower), res$n_trials - res$n_upper < 2)
})

test_that("rezdm validates parameters correctly", {
  # n must be single integer
  expect_error(
    rezdm(n = c(10, 20), n_trials = 100, drift = 2, bound = 1.5, ndt = 0.3),
    "n must be a single integer"
  )

  # n_trials must be larger than 2
  expect_error(
    rezdm(n = 10, n_trials = 1, drift = 2, bound = 1.5, ndt = 0.3),
    "n_trials must be larger than 2"
  )

  expect_error(
    rezdm(n = 10, n_trials = 2, drift = 2, bound = 1.5, ndt = 0.3),
    "n_trials must be larger than 2"
  )
})

test_that("dezdm handles edge case with near-zero drift", {
  # near-zero drift should not cause errors
  ll <- dezdm(
    mean_rt = 0.5, var_rt = 0.02, n_upper = 50, n_trials = 100,
    drift = 1e-8, bound = 1.5, ndt = 0.3
  )
  expect_true(is.finite(ll))
})

test_that("dezdm 4par handles edge case with near-zero drift", {
  # near-zero drift should not cause NaN in pC computation
  # Test with different zr values to ensure pC -> zr as drift -> 0

  # Test with zr = 0.5 (symmetric starting point)
  ll1 <- dezdm(
    mean_rt = c(0.5, 0.5), var_rt = c(0.02, 0.02),
    n_upper = 50, n_trials = 100,
    drift = 1e-8, bound = 1.5, ndt = 0.3, zr = 0.5, version = "4par"
  )
  expect_true(is.finite(ll1))

  # Test with zr = 0.3 (bias toward lower boundary)
  ll2 <- dezdm(
    mean_rt = c(0.5, 0.5), var_rt = c(0.02, 0.02),
    n_upper = 30, n_trials = 100,
    drift = 1e-8, bound = 1.5, ndt = 0.3, zr = 0.3, version = "4par"
  )
  expect_true(is.finite(ll2))

  # Test with zr = 0.7 (bias toward upper boundary)
  ll3 <- dezdm(
    mean_rt = c(0.5, 0.5), var_rt = c(0.02, 0.02),
    n_upper = 70, n_trials = 100,
    drift = 1e-8, bound = 1.5, ndt = 0.3, zr = 0.7, version = "4par"
  )
  expect_true(is.finite(ll3))

  # Test with exactly zero drift
  ll4 <- dezdm(
    mean_rt = c(0.5, 0.5), var_rt = c(0.02, 0.02),
    n_upper = 50, n_trials = 100,
    drift = 0, bound = 1.5, ndt = 0.3, zr = 0.5, version = "4par"
  )
  expect_true(is.finite(ll4))
})

test_that("dezdm 4par handles edge cases with few responses at boundary", {
  # when n_upper = 1, only binomial contributes (no mean/var for upper)
  ll <- dezdm(
    mean_rt = c(NA, 0.55), var_rt = c(NA, 0.025),
    n_upper = 1, n_trials = 100,
    drift = 2, bound = 1.5, ndt = 0.3, zr = 0.3, version = "4par"
  )
  expect_true(is.finite(ll))

  # when n_lower = 1 (n_upper = 99)
  ll <- dezdm(
    mean_rt = c(0.45, NA), var_rt = c(0.018, NA),
    n_upper = 99, n_trials = 100,
    drift = 2, bound = 1.5, ndt = 0.3, zr = 0.7, version = "4par"
  )
  expect_true(is.finite(ll))
})

# Changing a count to make a boundary sparse would change the binomial term
# too, so the reference keeps the counts and assembles the density from its
# parts: the binomial from the Wiener absorption probability, written out
# here, and each boundary's RT terms from the sampling distribution given in
# ?ezdm_dist, fed with the cumulants that the tests further down check against
# high-precision references. What is under test is which terms dezdm() adds;
# the first row, with every summary present, checks the reference itself.
test_that("dezdm 4par leaves out the RT terms of a boundary without summaries (#430)", {
  drift <- 1.2
  bound <- 1.4
  ndt <- 0.3
  zr <- 0.6
  n_upper <- 18
  n_trials <- 30
  moments <- .ezdm_moments_4par(drift, bound, zr, 1)
  mean_rt <- ndt + c(moments$mdt_upper * 1.03, moments$mdt_lower * 0.97)
  var_rt <- c(moments$vrt_upper * 1.1, moments$vrt_lower * 0.9)

  p_upper <- (1 - exp(-2 * drift * zr * bound)) / (1 - exp(-2 * drift * bound))
  binomial <- stats::dbinom(n_upper, n_trials, p_upper, log = TRUE)
  rt_terms <- function(mean_rt, var_rt, n, mdt, vrt, k3, k4) {
    W <- k4 / n + 2 * vrt^2 / (n - 1)
    stats::dgamma(var_rt, shape = vrt^2 / W, rate = vrt / W, log = TRUE) +
      stats::dnorm(mean_rt, mean = ndt + mdt + k3 / n / W * (var_rt - vrt),
                   sd = sqrt(vrt / n - (k3 / n)^2 / W), log = TRUE)
  }
  upper <- rt_terms(mean_rt[1], var_rt[1], n_upper, moments$mdt_upper,
                    moments$vrt_upper, moments$k3_upper, moments$k4_upper)
  lower <- rt_terms(mean_rt[2], var_rt[2], n_trials - n_upper,
                    moments$mdt_lower, moments$vrt_lower, moments$k3_lower,
                    moments$k4_lower)

  # rows: all present; no upper summaries; no lower summaries; none; a lower
  # variance without its mean; a lower mean without its variance
  mean_rt_obs <- matrix(c(
    mean_rt, NA, mean_rt[2], mean_rt[1], NA, NA, NA,
    mean_rt[1], NA, mean_rt
  ), ncol = 2, byrow = TRUE)
  var_rt_obs <- matrix(c(
    var_rt, NA, var_rt[2], var_rt[1], NA, NA, NA,
    var_rt, var_rt[1], NA
  ), ncol = 2, byrow = TRUE)
  ll <- dezdm(mean_rt_obs, var_rt_obs, n_upper, n_trials, drift = drift,
              bound = bound, ndt = ndt, zr = zr, version = "4par")

  expect_equal(
    ll, binomial + c(upper + lower, lower, upper, 0, upper, upper),
    tolerance = 1e-10
  )
})

test_that("generated data from rezdm has reasonable density under dezdm", {
  # generate data from known parameters
  set.seed(123)
  params <- list(drift = 2, bound = 1.5, ndt = 0.3, s = 1)
  sim_data <- rezdm(n = 1, n_trials = 100, drift = 2, bound = 1.5, ndt = 0.3)

  # evaluate density at generated point
  ll <- dezdm(
    mean_rt = sim_data$mean_rt,
    var_rt = sim_data$var_rt,
    n_upper = sim_data$n_upper,
    n_trials = sim_data$n_trials,
    drift = 2, bound = 1.5, ndt = 0.3
  )

  # density should be finite and reasonable (not extremely low)
  expect_true(is.finite(ll))
  expect_true(ll > -100)
})

# Vectorization tests for EZDM functions ----------------------------------------

test_that("dezdm 3par is vectorized over observations", {
  # single observation values
  ll1 <- dezdm(
    mean_rt = 0.5, var_rt = 0.02, n_upper = 80, n_trials = 100,
    drift = 2, bound = 1.5, ndt = 0.3
  )
  ll2 <- dezdm(
    mean_rt = 0.55, var_rt = 0.025, n_upper = 75, n_trials = 100,
    drift = 2, bound = 1.5, ndt = 0.3
  )

  # vectorized call
  ll_vec <- dezdm(
    mean_rt = c(0.5, 0.55),
    var_rt = c(0.02, 0.025),
    n_upper = c(80, 75),
    n_trials = c(100, 100),
    drift = 2, bound = 1.5, ndt = 0.3
  )

  expect_length(ll_vec, 2)
  expect_equal(ll_vec[1], ll1)
  expect_equal(ll_vec[2], ll2)
})

test_that("dezdm 3par recycles parameters correctly", {
  # generate test data
  set.seed(42)
  sim_data <- rezdm(n = 5, n_trials = 100, drift = 2, bound = 1.5, ndt = 0.3)

  # vectorized call with scalar parameters
  ll_vec <- dezdm(
    mean_rt = sim_data$mean_rt,
    var_rt = sim_data$var_rt,
    n_upper = sim_data$n_upper,
    n_trials = sim_data$n_trials,
    drift = 2, bound = 1.5, ndt = 0.3
  )

  expect_length(ll_vec, 5)
  expect_true(all(is.finite(ll_vec)))

  # loop-based verification
  ll_loop <- sapply(seq_len(nrow(sim_data)), function(i) {
    dezdm(
      mean_rt = sim_data$mean_rt[i],
      var_rt = sim_data$var_rt[i],
      n_upper = sim_data$n_upper[i],
      n_trials = sim_data$n_trials[i],
      drift = 2, bound = 1.5, ndt = 0.3
    )
  })

  expect_equal(ll_vec, ll_loop)
})

test_that("dezdm 4par works with matrix inputs", {
  # generate test data
  set.seed(42)
  sim_data <- rezdm(
    n = 5, n_trials = 100, drift = 2, bound = 1.5, ndt = 0.3,
    zr = 0.55, version = "4par"
  )

  # create matrices for mean_rt and var_rt
  mean_rt_mat <- cbind(sim_data$mean_rt_upper, sim_data$mean_rt_lower)
  var_rt_mat <- cbind(sim_data$var_rt_upper, sim_data$var_rt_lower)

  # vectorized call
  ll_vec <- dezdm(
    mean_rt = mean_rt_mat,
    var_rt = var_rt_mat,
    n_upper = sim_data$n_upper,
    n_trials = sim_data$n_trials,
    drift = 2, bound = 1.5, ndt = 0.3, zr = 0.55, version = "4par"
  )

  expect_length(ll_vec, 5)
  expect_true(all(is.finite(ll_vec)))

  # loop-based verification
  ll_loop <- sapply(seq_len(nrow(sim_data)), function(i) {
    dezdm(
      mean_rt = c(sim_data$mean_rt_upper[i], sim_data$mean_rt_lower[i]),
      var_rt = c(sim_data$var_rt_upper[i], sim_data$var_rt_lower[i]),
      n_upper = sim_data$n_upper[i],
      n_trials = sim_data$n_trials[i],
      drift = 2, bound = 1.5, ndt = 0.3, zr = 0.55, version = "4par"
    )
  })

  expect_equal(ll_vec, ll_loop)
})

test_that("dezdm 4par accepts length-2 vectors for single observation", {
  # length-2 vectors (backward compatibility)
  ll <- dezdm(
    mean_rt = c(0.45, 0.55),
    var_rt = c(0.018, 0.025),
    n_upper = 80, n_trials = 100,
    drift = 2, bound = 1.5, ndt = 0.3, zr = 0.55, version = "4par"
  )

  expect_length(ll, 1)
  expect_true(is.finite(ll))
})

test_that("dezdm 3par handles varying n_trials correctly", {
  # different n_trials for each observation
  ll_vec <- dezdm(
    mean_rt = c(0.5, 0.55, 0.52),
    var_rt = c(0.02, 0.025, 0.022),
    n_upper = c(80, 40, 60),
    n_trials = c(100, 50, 75),
    drift = 2, bound = 1.5, ndt = 0.3
  )

  expect_length(ll_vec, 3)
  expect_true(all(is.finite(ll_vec)))
})


# Tests for EZCDM distribution functions ---------------------------------------

test_that("dezcdm runs without errors and returns a finite log-likelihood", {
  ll <- dezcdm(
    mean_angle = 0.1, var_angle = 0.2, mean_rt = 0.9, var_rt = 0.1,
    n_trials = 100, driftrate = 2, bound = 1.5, ndt = 0.3
  )
  expect_type(ll, "double")
  expect_length(ll, 1)
  expect_true(is.finite(ll))
})

test_that("dezcdm log = FALSE is exp of log = TRUE", {
  args <- list(
    mean_angle = 0.1, var_angle = 0.4, mean_rt = 0.6, var_rt = 0.05,
    n_trials = 10, driftrate = 1, bound = 1, ndt = 0.3
  )
  expect_equal(
    do.call(dezcdm, c(args, log = FALSE)),
    exp(do.call(dezcdm, c(args, log = TRUE)))
  )
})

test_that("dezcdm is vectorized over observations and over parameter draws", {
  obs <- data.frame(
    mean_angle = c(0.1, -0.3, 2), var_angle = c(0.2, 0.5, 0.9),
    mean_rt = c(0.8, 1, 1.2), var_rt = c(0.1, 0.2, 0.3), n_trials = c(50, 100, 20)
  )
  ll_vec <- dezcdm(
    obs$mean_angle, obs$var_angle, obs$mean_rt, obs$var_rt, obs$n_trials,
    driftrate = 2, driftangle = 0.2, bound = 1.5, ndt = 0.3
  )
  ll_loop <- vapply(seq_len(nrow(obs)), function(i) {
    dezcdm(
      obs$mean_angle[i], obs$var_angle[i], obs$mean_rt[i], obs$var_rt[i], obs$n_trials[i],
      driftrate = 2, driftangle = 0.2, bound = 1.5, ndt = 0.3
    )
  }, numeric(1))
  expect_equal(ll_vec, ll_loop)

  draws <- c(1, 2, 3)
  ll_draws <- dezcdm(0.1, 0.2, 0.8, 0.1, 50, driftrate = draws, bound = 1.5, ndt = 0.3)
  ll_draws_loop <- vapply(draws, function(v) {
    dezcdm(0.1, 0.2, 0.8, 0.1, 50, driftrate = v, bound = 1.5, ndt = 0.3)
  }, numeric(1))
  expect_equal(ll_draws, ll_draws_loop)
})

test_that("dezcdm angle term is the exact von Mises log-likelihood of the trials", {
  withr::local_seed(2024)
  pars <- list(driftrate = 1.7, driftangle = 0.4, bound = 1.2, ndt = 0.25)
  kappa <- pars$bound * pars$driftrate
  theta <- brms::rvon_mises(80, pars$driftangle, kappa)
  summary <- .circular_summary(theta)
  mean_rt <- 0.9
  var_rt <- 0.2
  n <- length(theta)

  moments <- .ezcdm_moments(pars$driftrate, pars$bound, pars$ndt)
  rt <- .ez_rt_terms(moments$VRT, moments$k3, moments$k4, n)
  rt_terms <- stats::dgamma(var_rt, rt$shape, rt$rate, log = TRUE) +
    stats::dnorm(mean_rt, moments$MRT + rt$slope * (var_rt - moments$VRT), rt$sd, log = TRUE)
  angle_term <- dezcdm(
    summary$mean_angle, summary$var_angle, mean_rt, var_rt, n,
    driftrate = pars$driftrate, driftangle = pars$driftangle,
    bound = pars$bound, ndt = pars$ndt
  ) - rt_terms

  expect_equal(
    angle_term,
    sum(brms::dvon_mises(theta, pars$driftangle, kappa, log = TRUE)) + n * log(2 * pi)
  )
})

test_that(".ezcdm_moments matches the first-passage Laplace transform", {
  # E[exp(-sT)] = I0(a v) / I0(a sqrt(2 s + v^2)); moments by central differences
  laplace <- function(s, a, v) {
    exp(a * v - a * sqrt(2 * s + v^2)) *
      besselI(a * v, 0, TRUE) / besselI(a * sqrt(2 * s + v^2), 0, TRUE)
  }
  h <- 1e-4
  for (p in list(c(a = 1, v = 2), c(a = 2, v = 1), c(a = 1.5, v = 0.5), c(a = 3, v = 3))) {
    m1 <- -(laplace(h, p["a"], p["v"]) - laplace(-h, p["a"], p["v"])) / (2 * h)
    m2 <- (laplace(h, p["a"], p["v"]) - 2 + laplace(-h, p["a"], p["v"])) / h^2
    moments <- .ezcdm_moments(p["v"], p["a"], ndt = 0)
    expect_equal(unname(moments$MRT), unname(m1), tolerance = 1e-6)
    expect_equal(unname(moments$VRT), unname(m2 - m1^2), tolerance = 1e-5)
  }
})

test_that(".ezcdm_moments reaches the zero-drift limits", {
  moments <- .ezcdm_moments(driftrate = 1e-6, bound = 1.3, ndt = 0.2)
  expect_equal(moments$MRT, 0.2 + 1.3^2 / 2)
  expect_equal(moments$VRT, 1.3^4 / 8)
})

test_that(".ezcdm_moments small-kappa series matches the Bessel expressions", {
  kappa <- 0.0099
  R_direct <- besselI(kappa, 1, TRUE) / besselI(kappa, 0, TRUE)
  moments <- .ezcdm_moments(driftrate = kappa, bound = 1, ndt = 0)
  expect_equal(moments$R, R_direct, tolerance = 1e-9)
  expect_equal(moments$VRT / (1 / kappa)^2, R_direct^2 - 1 + 2 * R_direct / kappa, tolerance = 1e-6)
})

test_that(".ezcdm_moments k4 matches high-precision references in every branch", {
  # 60-digit mpmath values of the fourth derivative at 0 of
  # K(l) = log I0(a v) - log I0(a sqrt(v^2 + 2 l)); kappa = a v covers the
  # small-kappa series, the closed form and the asymptotic series
  ref <- data.frame(
    a = c(0.5, 1, 1, 1.5, 1.8, 3, 2, 4, 2.5),
    v = c(0.02, 0.05, 0.5, 2, 1.5, 3, 10, 50, 320),
    k4 = c(
      0.0003356701670180080061, 0.085789222748270863285, 0.072568765521756091338,
      0.052584022947964774569, 0.36764406706077275354, 0.016629217455392289592,
      2.7527454971494603793e-6, 7.618390443994564079e-11, 1.0892093592819363706e-16
    )
  )
  moments <- .ezcdm_moments(ref$v, ref$a, ndt = 0)
  expect_equal(moments$k4, ref$k4, tolerance = 1e-10)
})

test_that(".ezcdm_moments k4 reaches the zero-drift limit", {
  # log I0(z) = z^2/4 - z^4/64 + z^6/576 - 11 z^8/49152 + ...  gives k4 = 11 a^8 / 128
  moments <- .ezcdm_moments(driftrate = 1e-6, bound = 1.3, ndt = 0.2)
  expect_equal(moments$k4, 11 * 1.3^8 / 128)
})

test_that(".ezcdm_moments k4 is continuous across its series branches", {
  for (kappa in c(0.25, 100)) {
    below <- .ezcdm_moments(driftrate = kappa * (1 - 1e-12), bound = 1, ndt = 0)$k4
    above <- .ezcdm_moments(driftrate = kappa * (1 + 1e-12), bound = 1, ndt = 0)$k4
    expect_lt(abs(above / below - 1), 1e-9)
  }
})

test_that(".ezcdm_moments k3 matches high-precision references in every branch", {
  # 60-digit mpmath values of minus the third derivative at 0 of
  # K(l) = log I0(a v) - log I0(a sqrt(v^2 + 2 l)); kappa = a v from 0.01 to 5000
  # covers the small-kappa series, the closed form and the asymptotic series
  ref <- data.frame(
    a = c(0.5, 1, 0.4, 1, 1.8, 1.5, 3, 2, 1.2, 4, 2.5, 2),
    v = c(0.02, 0.05, 0.5, 0.5, 1.5, 2, 3, 10, 50, 50, 320, 2500),
    k3 = c(
      0.001302016196980727444719, 0.08322600416498629405255, 0.0003343895036682916188472,
      0.073456137602591188782, 0.2548940767246537981923, 0.05976258840288373705414,
      0.03119349519173859704521, 0.00005589774778980381167081, 1.126194453909258777992e-8,
      3.814339514451903375109e-8, 2.231446701613348675612e-12, 6.142361446350825587212e-17
    )
  )
  moments <- .ezcdm_moments(ref$v, ref$a, ndt = 0)
  expect_equal(moments$k3, ref$k3, tolerance = 1e-11)
})

test_that(".ezcdm_moments k3 reaches the zero-drift limit", {
  # log I0(z) = z^2/4 - z^4/64 + z^6/576 - ...  gives k3 = a^6 / 12
  moments <- .ezcdm_moments(driftrate = 1e-6, bound = 1.3, ndt = 0.2)
  expect_equal(moments$k3, 1.3^6 / 12)
})

test_that(".ezcdm_moments k3 is continuous across its series branches", {
  for (kappa in c(0.25, 100)) {
    below <- .ezcdm_moments(driftrate = kappa * (1 - 1e-12), bound = 1, ndt = 0)$k3
    above <- .ezcdm_moments(driftrate = kappa * (1 + 1e-12), bound = 1, ndt = 0)$k3
    expect_lt(abs(above / below - 1), 1e-10)
  }
})

test_that(".ezcdm_moments large-kappa expansion matches 50-digit Bessel references", {
  # log I0(kappa) and I1/I0 from mpmath with 50 significant digits
  ref <- data.frame(
    kappa = c(1e3, 1e4, 2e5, 1e6),
    log_I0 = c(995.627308889869464671467764481, 9994.47590378143230100450870026,
               199992.97802576903180290160416, 999992.1733063128132527062308),
    R = c(0.999499874874804280198918174348, 0.999949998749874980464686451826,
          0.999997499996874984374877928418, 0.999999499999874999874999804687)
  )
  moments <- .ezcdm_moments(driftrate = ref$kappa, bound = 1, ndt = 0)
  expect_equal(moments$log_I0, ref$log_I0, tolerance = 1e-14)
  expect_equal(moments$R, ref$R, tolerance = 1e-14)

  below <- .ezcdm_moments(driftrate = 1000 * (1 - 1e-12), bound = 1, ndt = 0)
  above <- .ezcdm_moments(driftrate = 1000 * (1 + 1e-12), bound = 1, ndt = 0)
  expect_lt(abs(above$R / below$R - 1), 1e-12)
  expect_lt(abs(above$VRT / below$VRT - 1), 1e-9)
})

test_that(".ezcdm_moments stays finite at extreme drift rates", {
  # driftrate = 0 stands in for exp() of a very negative linear predictor
  tiny <- .ezcdm_moments(driftrate = c(0, 1e-300, 1e-200), bound = 1.5, ndt = 0.3)
  expect_true(all(vapply(tiny, function(x) all(is.finite(x)), logical(1))))
  expect_equal(tiny$VRT, rep(1.5^4 / 8, 3))
  expect_equal(tiny$MRT, rep(0.3 + 1.5^2 / 2, 3))

  huge <- .ezcdm_moments(driftrate = c(1e6, 1e8), bound = 1, ndt = 0.3)
  expect_true(all(vapply(huge, function(x) all(is.finite(x)), logical(1))))
  expect_true(is.finite(dezcdm(0.1, 0.2, 0.9, 0.1, 100, driftrate = 1e-300, bound = 1, ndt = 0.3)))
  expect_true(is.finite(dezcdm(0.1, 0.2, 0.3, 1e-12, 100, driftrate = 1e3, bound = 1e3, ndt = 0.3)))
})

test_that("the skew bounds quoted in ?ezcdm_dist hold over the kappa range", {
  kappa <- 10^seq(-4, 5, length.out = 400)
  # driftrate = bound = sqrt(kappa) fixes tau = 1; the ratios do not depend on tau
  moments <- .ezcdm_moments(driftrate = sqrt(kappa), bound = sqrt(kappa), ndt = 0)
  excess <- moments$k4 / moments$VRT^2
  # large-n limits: Var(var_rt) relative to the normal-RT value, and Cor(mean_rt, var_rt)
  variance_ratio <- 1 + excess / 2
  correlation <- moments$k3 / sqrt(moments$VRT * (moments$k4 + 2 * moments$VRT^2))

  expect_true(all(excess >= 0))
  expect_equal(max(3 + excess), 8.5, tolerance = 1e-3)
  expect_lte(max(3 + excess), 8.5)
  expect_equal(max(variance_ratio), 3.75, tolerance = 1e-3)
  expect_lte(max(variance_ratio), 3.75)
  expect_gt(max(correlation), 0.68)
  expect_lt(max(correlation), 0.69)
})

test_that("the conditional sd of mean_rt is positive over the kappa range", {
  kappa <- 10^seq(-8, 5, length.out = 300)
  moments <- .ezcdm_moments(driftrate = kappa, bound = 1, ndt = 0)
  for (n in c(3, 10, 1000)) {
    rt <- .ez_rt_terms(moments$VRT, moments$k3, moments$k4, n)
    expect_true(all(is.finite(rt$sd) & rt$sd^2 > 0.5 * moments$VRT / n))
  }
})

test_that("dezcdm with k3 = 0 is the product of independent mean_rt and var_rt terms", {
  moments_full <- .ezcdm_moments
  local_mocked_bindings(.ezcdm_moments = function(...) {
    moments <- moments_full(...)
    moments$k3 <- 0 * moments$k3
    moments
  })
  obs <- data.frame(mean_angle = c(0.1, -0.4), var_angle = c(0.2, 0.5), mean_rt = c(0.9, 1.4),
                    var_rt = c(0.1, 0.35), n_trials = c(100, 20))
  pars <- list(driftrate = c(2, 0.7), driftangle = 0.1, bound = c(1.5, 2.2), ndt = 0.3)
  ll <- dezcdm(obs$mean_angle, obs$var_angle, obs$mean_rt, obs$var_rt, obs$n_trials,
               driftrate = pars$driftrate, driftangle = pars$driftangle,
               bound = pars$bound, ndt = pars$ndt)

  moments <- moments_full(pars$driftrate, pars$bound, pars$ndt)
  W <- moments$k4 / obs$n_trials + 2 * moments$VRT^2 / (obs$n_trials - 1)
  independent <- obs$n_trials * (moments$kappa * (1 - obs$var_angle) *
    cos(obs$mean_angle - pars$driftangle) - moments$log_I0) +
    stats::dnorm(obs$mean_rt, moments$MRT, sqrt(moments$VRT / obs$n_trials), log = TRUE) +
    stats::dgamma(obs$var_rt, moments$VRT^2 / W, moments$VRT / W, log = TRUE)
  expect_equal(ll, independent)
})

test_that("dezcdm is continuous across the small-kappa branch", {
  ll <- function(driftrate) {
    dezcdm(0.1, 0.9, 0.8, 0.125, 100, driftrate = driftrate, bound = 1, ndt = 0.3)
  }
  below <- ll(0.01 * (1 - 1e-9))
  above <- ll(0.01 * (1 + 1e-9))
  expect_true(is.finite(below) && is.finite(above))
  expect_lt(abs(below - above), 1e-6)
})

test_that("dezcdm validates parameters", {
  args <- list(
    mean_angle = 0.1, var_angle = 0.2, mean_rt = 0.9, var_rt = 0.1,
    n_trials = 100, driftrate = 2, bound = 1.5, ndt = 0.3
  )
  bad <- function(...) do.call(dezcdm, utils::modifyList(args, list(...)))
  expect_error(bad(driftrate = 0), "driftrate must be positive")
  expect_error(bad(bound = -1), "bound must be positive")
  expect_error(bad(ndt = -0.1), "ndt must be non-negative")
  expect_error(bad(n_trials = 2), "n_trials must be larger than 2")
  expect_error(bad(var_angle = 1.1), "var_angle must be between 0 and 1")
  expect_error(bad(var_angle = -0.1), "var_angle must be between 0 and 1")
  expect_error(bad(var_rt = 0), "var_rt must be positive")
})

test_that("rezcdm returns plausible summary statistics", {
  withr::local_seed(123)
  sim <- rezcdm(n = 200, n_trials = 50, driftrate = 2, driftangle = 0.5, bound = 1.5, ndt = 0.3)

  expect_s3_class(sim, "data.frame")
  expect_equal(nrow(sim), 200)
  expect_named(sim, c("mean_angle", "var_angle", "mean_rt", "var_rt", "n_trials"))
  expect_true(all(sim$mean_rt > 0.3))
  expect_true(all(sim$var_angle >= 0 & sim$var_angle <= 1))
  expect_true(all(sim$var_rt > 0))
  expect_true(all(abs(sim$mean_angle) <= pi))
  expect_true(all(sim$n_trials == 50))
})

test_that("rezcdm reproduces the moments used by dezcdm", {
  withr::local_seed(99)
  pars <- list(driftrate = 1.5, driftangle = -0.6, bound = 2, ndt = 0.2)
  n_trials <- 40
  sim <- do.call(rezcdm, c(list(n = 4000, n_trials = n_trials), pars))
  moments <- .ezcdm_moments(pars$driftrate, pars$bound, pars$ndt)

  expect_equal(mean(sim$mean_rt), moments$MRT, tolerance = 0.01)
  expect_equal(var(sim$mean_rt), moments$VRT / n_trials, tolerance = 0.05)
  expect_equal(mean(sim$var_rt), moments$VRT, tolerance = 0.02)
  expect_equal(
    var(sim$var_rt),
    moments$k4 / n_trials + 2 * moments$VRT^2 / (n_trials - 1),
    tolerance = 0.1
  )
  expect_equal(cov(sim$mean_rt, sim$var_rt), moments$k3 / n_trials, tolerance = 0.1)
  resultant <- mean((1 - sim$var_angle) * cos(sim$mean_angle - pars$driftangle))
  expect_equal(atan2(mean(sin(sim$mean_angle)), mean(cos(sim$mean_angle))), pars$driftangle, tolerance = 0.02)
  expect_equal(resultant, moments$R, tolerance = 0.02)
})

test_that("rezcdm does not truncate mean_rt, and values below ndt stay rare", {
  withr::local_seed(11)
  sim <- rezcdm(n = 20000, n_trials = 3, driftrate = 0.03, bound = 0.03, ndt = 0.3)
  below <- mean(sim$mean_rt < 0.3)
  expect_gt(below, 0)
  expect_lt(below, 0.005)
})

test_that("rezcdm recycles parameters and n_trials across replicates", {
  withr::local_seed(1)
  sim <- rezcdm(n = 4, n_trials = c(10, 20), driftrate = c(1, 3), bound = 1, ndt = c(0.2, 0.4))
  expect_equal(sim$n_trials, c(10, 20, 10, 20))
  expect_true(all(sim$mean_rt > c(0.2, 0.4, 0.2, 0.4)))
})

test_that("rezcdm validates parameters", {
  expect_error(rezcdm(n = c(1, 2), n_trials = 10, driftrate = 1, bound = 1, ndt = 0.3), "n must be")
  expect_error(rezcdm(n = 5, n_trials = 2, driftrate = 1, bound = 1, ndt = 0.3), "n_trials must be larger than 2")
  expect_error(rezcdm(n = 5, n_trials = 10.5, driftrate = 1, bound = 1, ndt = 0.3), "whole numbers")
  expect_error(rezcdm(n = 2, n_trials = c(10, NA), driftrate = 1, bound = 1, ndt = 0.3), "without missing values")
  expect_error(rezcdm(n = 5, n_trials = 10, driftrate = -1, bound = 1, ndt = 0.3), "driftrate must be positive")
  expect_error(rezcdm(n = 5, n_trials = 10, driftrate = 1, bound = 0, ndt = 0.3), "bound must be positive")
  expect_error(rezcdm(n = 5, n_trials = 10, driftrate = 1, bound = 1, ndt = -0.3), "ndt must be non-negative")
})

# Tests for the ezdm decision-time cumulants (issue #407) ----------------------

# Reference values from local/ezdm/k34_derivation.py: the cumulants of the
# boundary-conditional first-passage time, evaluated at 60+ digits from the
# closed forms with enough guard digits for the cancellation, and cross-checked
# against numerical differentiation of the CGF itself (agreement 7e-59).
ezdm_cumulant_references <- function() {
  # read.csv rather than a tribble: tibble is not a declared dependency, and the
  # 17-digit decimals round-trip the doubles exactly
  ref <- utils::read.csv(text = "drift,bound,zr,s,boundary,MDT,VRT,k3,k4
0,1.5,0.5,1.0,upper,0.5625,0.2109375,0.18984375,0.25934012276785714
1e-08,1.5,0.5,1.0,upper,0.56249999999999999,0.21093749999999999,0.18984374999999999,0.25934012276785712
0.01,1.0,0.5,1.0,upper,0.24999791668749979,0.041665833345981972,0.016666160724536864,0.01011863757652396
0.2,1.5,0.5,1.0,upper,0.5583188760874424,0.20719186958160496,0.1847501123616437,0.25009358517694253
0.5,1.5,0.5,1.0,upper,0.53753609752617892,0.18908944624838941,0.16079771187077391,0.2077761272089794
2.0,1.0,0.5,1.0,upper,0.19039853898894122,0.021351238396358676,0.0060181161652498743,0.0025955878267946837
5.0,3.0,0.5,1.4,upper,0.29971538162917291,0.023326995904847908,0.0053841925909205833,0.0020413082600528851
100.0,1.5,0.5,1.0,upper,0.0075,7.5e-7,2.25e-10,1.125e-13
0.05,1.5,0.8,1.0,upper,0.2698340611020102,0.13266229314104025,0.14203028421087174,0.21614906104371065
0.05,1.5,0.8,1.0,lower,0.71971935057529668,0.22439914724951283,0.19251972387929096,0.2597655971356454
1.0,1.5,0.8,1.0,upper,0.21774203644294689,0.08200902526658194,0.071383919789698674,0.089504207727222578
1.0,1.5,0.8,1.0,lower,0.62736556037724539,0.15310222810508356,0.10517772494355196,0.11488864600247718
2.0,0.5,0.95,1.4,upper,0.004013234648407409,0.0001250812152839004,8.439155091234556e-6,7.9963594687025653e-7
2.0,0.5,0.95,1.4,lower,0.041690711753955124,0.0006885741617546384,3.252941531017657e-5,2.4237919043289123e-6
20.0,3.0,0.6,1.0,upper,0.060000000000000003,0.00015000000000000001,1.1250000000000001e-6,1.4062500000000001e-8
1500.0,1.0,0.6,1.0,upper,0.00026666666666666668,1.1851851851851853e-10,1.580246913580247e-16,3.5116598079561044e-22")
  ref$b <- ifelse(ref$boundary == "upper", ref$zr, 1 - ref$zr) * ref$bound
  ref
}

test_that(".ezdm_cumulants matches high-precision references in both branches", {
  ref <- ezdm_cumulant_references()
  got <- .ezdm_cumulants(ref$b, ref$bound, ref$drift^2 / ref$s^4, ref$s)

  # row-wise ratios: expect_equal() pools its tolerance over the vector, and the
  # references span 20 orders of magnitude
  for (moment in c("MDT", "VRT", "k3", "k4")) {
    expect_equal(got[[moment]] / ref[[moment]], rep(1, nrow(ref)), tolerance = 1e-10, info = moment)
  }
})

test_that(".ezdm_cumulants is continuous across the series/closed-form seam", {
  for (bound in c(0.8, 1.5, 3)) {
    for (zr in c(0.5, 0.7, 0.9)) {
      w_seam <- (0.7 / bound)^2
      below <- .ezdm_cumulants(zr * bound, bound, w_seam * (1 - 1e-12), 1)
      above <- .ezdm_cumulants(zr * bound, bound, w_seam * (1 + 1e-12), 1)
      expect_lt(max(abs(unlist(above) / unlist(below) - 1)), 1e-9)
    }
  }
})

test_that(".ezdm_cumulants reproduces the zero-drift limits", {
  # k_n(w = 0) = (-1)^n (2/s^2)^n a_n n! (b^2n - b0^2n); for the symmetric start
  # point these collapse to the 3par values below
  a <- 1.3
  zero_drift <- .ezdm_cumulants(a / 2, a, 0, 1)
  expect_equal(zero_drift$MDT, a^2 / 4)
  expect_equal(zero_drift$VRT, a^4 / 24)
  expect_equal(zero_drift$k3, a^6 / 60)
  expect_equal(zero_drift$k4, 17 * a^8 / 1680)

  # excess kurtosis 17 * 576 / 1680 at zero drift, i.e. kurtosis 8.83
  expect_equal(zero_drift$k4 / zero_drift$VRT^2, 17 * 576 / 1680)
})

test_that(".ezdm_cumulants reproduces the 4par zero-drift limits", {
  bound <- 1.6
  zr <- 0.7
  z <- bound / 2
  x0 <- zr * bound - z
  zero_drift <- .ezdm_cumulants(c(z + x0, z - x0), bound, 0, 1)

  expect_equal(
    zero_drift$MDT,
    c(4 * z^2 - (z + x0)^2, 4 * z^2 - (z - x0)^2) / 3
  )
  expect_equal(
    zero_drift$VRT,
    c(32 * z^4 - 2 * (z + x0)^4, 32 * z^4 - 2 * (z - x0)^4) / 45
  )
})

test_that(".ezdm_cumulants agrees with the eigenfunction series of the first-passage density", {
  # A referee that shares nothing with the log-sinh CGF: the large-time series
  # f(t) ~ sum_k k sin(k pi x / a) exp(-c_k t), c_k = v^2 / 2 + k^2 pi^2 / (2 a^2),
  # for s = 1 and x the distance to the boundary that is hit, integrated term by
  # term (int t^m exp(-c t) dt = m! / c^(m + 1)). The m = 0 sum converges like
  # 1 / k, so sum sin(k phi) / k = (pi - phi) / 2 is split off.
  eigen_cumulants <- function(x, a, v, terms = 4e5) {
    k <- seq_len(terms)
    phi <- pi * x / a
    beta_sq <- (v * a / pi)^2
    sin_k <- sin(k * phi)
    c_k <- v^2 / 2 + k^2 * pi^2 / (2 * a^2)
    mass <- 2 * a^2 / pi^2 *
      ((pi - phi) / 2 - beta_sq * sum(sin_k / (k * (k^2 + beta_sq))))
    mu <- vapply(1:4, \(m) factorial(m) * sum(k * sin_k / c_k^(m + 1)), numeric(1)) / mass
    c(
      MDT = mu[1],
      VRT = mu[2] - mu[1]^2,
      k3 = mu[3] - 3 * mu[2] * mu[1] + 2 * mu[1]^3,
      k4 = mu[4] - 4 * mu[3] * mu[1] - 3 * mu[2]^2 + 12 * mu[2] * mu[1]^2 - 6 * mu[1]^4
    )
  }

  # series branch (bound * |drift| / s^2 = 0.42), zero drift, and the closed
  # forms at s != 1 (3.28); (drift, bound, s) is (drift / s, bound / s, 1)
  cells <- data.frame(drift = c(0.3, 0, -1.5), bound = 1.4, s = c(1, 1, 0.8))
  zr <- 0.65
  for (i in seq_len(nrow(cells))) {
    drift <- cells$drift[i]
    bound <- cells$bound[i]
    s <- cells$s[i]
    got <- .ezdm_cumulants(c(zr, 1 - zr) * bound, bound, drift^2 / s^4, s)
    for (side in 1:2) {
      hit_distance <- c(1 - zr, zr)[side] * bound
      expected <- eigen_cumulants(hit_distance / s, bound / s, drift / s)
      # as ratios, so an error in k4 alone is not pooled with the larger MDT
      expect_equal(
        vapply(got, `[`, numeric(1), side) / expected, c(MDT = 1, VRT = 1, k3 = 1, k4 = 1),
        tolerance = 1e-7
      )
    }
  }
})

test_that(".ezdm_cumulants stays finite at extreme drift", {
  # cosh(t)/sinh(t) is NaN above t = 710 and t^4 * csch^2(t) is Inf * 0 there;
  # the old .ezdm_moments_4par() returned NaN for every moment at drift = 1500
  extreme <- .ezdm_cumulants(0.6, 1, c(0, 1e-300, 1e-8, 1, 1e6, 1e12)^2, 1)
  expect_true(all(is.finite(unlist(extreme))))

  moments <- .ezdm_moments_4par(drift = 1500, bound = 1, zr = 0.6, s = 1)
  expect_true(all(is.finite(unlist(moments))))
  expect_true(all(is.finite(unlist(.ezdm_moments_3par(1500, 1, 1)))))
})

test_that(".ez_rt_terms keeps a positive conditional variance across the ezdm grid", {
  # n = 2 is the smallest count the 4par likelihood admits at one boundary
  grid <- expand.grid(
    drift = 10^seq(-3, 2, length.out = 25), bound = c(0.5, 1, 2, 4),
    zr = c(0.5, 0.7, 0.9), s = c(0.3, 1, 3), n = c(2, 3, 5, 10, 50, 200, 1000)
  )
  moments <- .ezdm_cumulants(
    grid$zr * grid$bound, grid$bound, grid$drift^2 / grid$s^4, grid$s
  )
  rt <- .ez_rt_terms(moments$VRT, moments$k3, moments$k4, grid$n)

  ratio <- rt$sd^2 / (moments$VRT / grid$n)
  expect_true(all(is.finite(ratio)))
  # the conditional sd never drops far below the marginal one: corr^2 < 0.7
  expect_gt(min(ratio), 0.3)
  expect_lte(max(ratio), 1)
})

test_that(".ez_rt_terms reduces to the independent normal and chi-square without skew", {
  rt <- .ez_rt_terms(VRT = 0.3, k3 = 0, k4 = 0, n_trials = 25)

  expect_equal(rt$shape, (25 - 1) / 2)
  expect_equal(rt$rate, (25 - 1) / (2 * 0.3))
  expect_equal(rt$slope, 0)
  expect_equal(rt$sd, sqrt(0.3 / 25))
})

test_that(".ez_rt_terms reproduces the exact moments of mean_rt and var_rt", {
  # Var(sample variance) = k4 / n + 2 VRT^2 / (n - 1) and
  # Cov(sample mean, sample variance) = k3 / n, for any distribution
  n <- c(2, 5, 25, 400)
  moments <- .ezdm_moments_4par(drift = 1.3, bound = 1.4, zr = 0.7, s = 1)
  for (side in c("upper", "lower")) {
    VRT <- moments[[paste0("vrt_", side)]]
    k3 <- moments[[paste0("k3_", side)]]
    k4 <- moments[[paste0("k4_", side)]]
    rt <- .ez_rt_terms(VRT, k3, k4, n)
    var_of_var <- rt$shape / rt$rate^2

    expect_equal(rt$shape / rt$rate, rep(VRT, 4))
    expect_equal(var_of_var, k4 / n + 2 * VRT^2 / (n - 1))
    expect_equal(rt$slope * var_of_var, k3 / n)
    expect_equal(rt$sd^2 + rt$slope^2 * var_of_var, VRT / n)
  }
})

test_that(".ez_rt_terms matches the sampling moments of simulated diffusion summaries", {
  # a referee that shares no formula with bmm: the empirical moments of the
  # mean and variance of n decision times drawn by rtdists
  skip_on_cran()
  skip_if_not_installed("rtdists")
  withr::local_seed(407)

  n <- 10
  drift <- 1.3
  bound <- 1.4
  decision_times <- matrix(
    rtdists::rdiffusion(n * 1e5, a = bound, v = drift, t0 = 0, z = bound / 2)$rt,
    nrow = n
  )
  mean_rt <- colMeans(decision_times)
  var_rt <- apply(decision_times, 2, stats::var)

  moments <- .ezdm_moments_3par(drift, bound, 1)
  rt <- .ez_rt_terms(moments$VRT, moments$k3, moments$k4, n)
  var_of_var <- rt$shape / rt$rate^2

  # Monte Carlo error is about 1%; the pre-#407 terms miss var(var_rt) by 3.6x
  # and cov(mean_rt, var_rt) entirely. Compared as ratios, because
  # expect_equal() treats `tolerance` as absolute when the target is smaller
  # than it, and these moments are ~0.005
  expect_equal(var_of_var / stats::var(var_rt), 1, tolerance = 0.05)
  expect_equal(rt$slope * var_of_var / stats::cov(mean_rt, var_rt), 1, tolerance = 0.05)
  expect_equal((rt$sd^2 + rt$slope^2 * var_of_var) / stats::var(mean_rt), 1, tolerance = 0.05)
})

test_that("dezdm 4par uses a boundary's RT summaries only from two responses on", {
  lpdf <- function(n_upper, mean_rt, var_rt) {
    dezdm(
      mean_rt = mean_rt, var_rt = var_rt, n_upper = n_upper, n_trials = 30,
      drift = 1, bound = 1.5, ndt = 0.25, zr = 0.6, version = "4par"
    )
  }

  expect_equal(lpdf(1, c(NA, 0.8), c(NA, 0.1)), lpdf(1, c(0.7, 0.8), c(0.05, 0.1)))
  expect_equal(lpdf(29, c(0.7, NA), c(0.05, NA)), lpdf(29, c(0.7, 0.8), c(0.05, 0.1)))
  expect_true(is.finite(lpdf(1, c(NA, 0.8), c(NA, 0.1))))
  expect_true(is.finite(lpdf(2, c(0.7, 0.8), c(0.05, 0.1))))
  expect_false(lpdf(2, c(0.7, 0.8), c(0.05, 0.1)) == lpdf(2, c(0.9, 0.8), c(0.05, 0.1)))
})

test_that(".ezdm_pc is stable in both drift directions", {
  # the old exp(2 k z) form returned NaN from |drift| * bound / s^2 ~ 710 at
  # negative drift
  drift <- c(-1500, -400, -1, 0, 1, 400, 1500)
  pC <- .ezdm_moments_4par(drift = drift, bound = 1, zr = 0.6, s = 1)$pC

  expect_true(all(is.finite(pC)))
  expect_true(all(pC >= 0 & pC <= 1))
  expect_equal(pC[drift == 0], 0.6)
  expect_true(all(diff(pC) >= 0))

  # 3par is the same function at zr = 0.5
  expect_equal(
    .ezdm_moments_3par(drift, 1, 1)$pC,
    .ezdm_moments_4par(drift, bound = 1, zr = 0.5, s = 1)$pC
  )
})

test_that(".ezdm_pc keeps its slope in drift through zero", {
  # a constant b / b0 at k = 0 has the right value but no gradient, which is
  # what a sampler started at drift = 0 would see; d pC / d k = r (b0 - b) there
  b <- 0.45
  b0 <- 1.5
  h <- 1e-7
  slope <- (.ezdm_pc(b, b0, h) - .ezdm_pc(b, b0, -h)) / (2 * h)
  expect_equal(slope / (b / b0 * (b0 - b)), 1, tolerance = 1e-6)
  expect_equal(.ezdm_pc(b, b0, 0), b / b0)

  # continuous across the switch between the series and the expm1 ratio
  k_switch <- 1e-4 / (2 * b0)
  around <- .ezdm_pc(b, b0, c(-1, -1, 1, 1) * k_switch * (1 + c(1, -1, -1, 1) * 1e-9))
  expect_equal(around[1] / around[2], 1, tolerance = 1e-12)
  expect_equal(around[3] / around[4], 1, tolerance = 1e-12)

  # negative drift mirrors positive drift at the mirrored start point
  expect_equal(.ezdm_pc(b, b0, -1.3), 1 - .ezdm_pc(b0 - b, b0, 1.3))
})

test_that("the 3par density takes its binomial on the logit scale drift * bound / s^2", {
  # the identity the Stan likelihood and .dezdm_3par() both rely on
  drift <- c(-3, -0.4, 0, 0.7, 2.5)
  bound <- c(0.6, 1.2, 2, 1.5, 0.9)
  s <- c(1, 0.7, 1.3, 2, 1)
  expect_equal(stats::qlogis(.ezdm_moments_3par(drift, bound, s)$pC), drift * bound / s^2)

  # and the density agrees with the probability-scale binomial where that is exact
  moments <- .ezdm_moments_3par(1.2, 1.5, 1)
  with_logit <- dezdm(0.6, 0.05, n_upper = 40, n_trials = 50, drift = 1.2, bound = 1.5, ndt = 0.25)
  rt <- .ez_rt_terms(moments$VRT, moments$k3, moments$k4, 50)
  with_prob <- stats::dbinom(40, 50, moments$pC, log = TRUE) +
    stats::dgamma(0.05, rt$shape, rt$rate, log = TRUE) +
    stats::dnorm(0.6, 0.25 + moments$MDT + rt$slope * (0.05 - moments$VRT), rt$sd, log = TRUE)
  expect_equal(with_logit, with_prob)

  # pC rounds to 1 here, and the old form returned -Inf
  expect_true(is.finite(dezdm(0.3, 1e-4, n_upper = 49, n_trials = 50, drift = 60, bound = 1.5, ndt = 0.25)))
})

test_that("the 4par density is the same when drift, zr and the boundaries are mirrored", {
  # log(1 - pC) taken of a pC that had rounded to 1 broke this from
  # |drift| * bound / s^2 of about 18 (drift = 10 here is 31), and returned -Inf
  # on one side only from about 37
  grid <- expand.grid(
    drift = c(-25, -10, -5, -1.2, -0.3, 0.3, 1.2, 5, 10, 25), zr = c(0.2, 0.7),
    n_upper = c(0, 1, 59, 60)
  )
  moments <- .ezdm_moments_4par(grid$drift, 1.5, grid$zr, 0.7)
  mean_rt <- 0.25 + cbind(moments$mdt_upper * 1.03, moments$mdt_lower * 0.97)
  var_rt <- cbind(moments$vrt_upper * 0.9, moments$vrt_lower * 1.2)

  as_coded <- dezdm(mean_rt, var_rt, grid$n_upper, 60,
    drift = grid$drift, bound = 1.5, ndt = 0.25, zr = grid$zr, s = 0.7, version = "4par"
  )
  mirrored <- dezdm(mean_rt[, 2:1], var_rt[, 2:1], 60 - grid$n_upper, 60,
    drift = -grid$drift, bound = 1.5, ndt = 0.25, zr = 1 - grid$zr, s = 0.7, version = "4par"
  )

  expect_true(all(is.finite(as_coded)))
  expect_lt(max(abs(mirrored / as_coded - 1)), 1e-11)
})

test_that(".ezdm_logit_pc is the logit of .ezdm_pc and the 3par logit at zr = 0.5", {
  bound <- 1.5
  k <- c(-40, -20, -5, 5, 20, 40) / bound
  expect_lt(max(abs(.ezdm_logit_pc(bound / 2, bound / 2, k) / (k * bound) - 1)), 1e-13)

  # where pC is exact, on both sides of the switch between the series and the
  # closed form
  k <- c(-3, -1e-3, -1e-6, 0, 1e-6, 1e-3, 3)
  logit <- .ezdm_logit_pc(0.45, bound - 0.45, k)
  expect_lt(max(abs(logit / stats::qlogis(.ezdm_pc(0.45, bound, k)) - 1)), 1e-13)

  # pC rounds to 1 here and qlogis() of it is Inf
  expect_equal(.ezdm_logit_pc(0.45, bound - 0.45, 30), 2 * 30 * 0.45)
})
