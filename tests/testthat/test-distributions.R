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

  # zr out of range (note: actual error message has typo "startin point")
  expect_error(
    dcswald(1, 1, drift = 2, bound = 1.5, ndt = 0.3, zr = 1.5),
    "startin"
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
