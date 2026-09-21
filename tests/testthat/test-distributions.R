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
