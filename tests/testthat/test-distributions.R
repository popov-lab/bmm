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
  expect_equal(integrate(dsdm, -pi, pi, mu = 0, c = 3, kappa = 3)$value, 1)
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
  )$value, 1)
})

test_that("dmixture3p integrates to 1", {
  expect_equal(integrate(dmixture3p, -pi, pi,
    mu = runif(3, min = -pi, pi),
    kappa = runif(1, min = 1, max = 20),
    p_mem = runif(1, min = 0, max = 0.6),
    p_nt = runif(1, min = 0, max = 0.3)
  )$value, 1)
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
  )$value, 1)
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
    "No activation functions"
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
  res <- dcswald(rt, response, drift = 2, bound = 1.5, ndt = 0.3,
                 version = "simple")
  expect_length(res, n)
  expect_type(res, "double")

  # crisk version
  res_crisk <- dcswald(rt, response, drift = 2, bound = 1.5, ndt = 0.3,
                       zr = 0.5, version = "crisk")
  expect_length(res_crisk, n)
  expect_type(res_crisk, "double")

  # test log = FALSE
  res_exp <- dcswald(rt, response, drift = 2, bound = 1.5, ndt = 0.3,
                     version = "simple", log = FALSE)
  expect_true(all(res_exp >= 0))
  expect_equal(res_exp, exp(res))
})

test_that("dcswald handles vectorized parameters", {
  n <- 5
  rt <- runif(n, 0.4, 1.5)
  response <- sample(c(0, 1), n, replace = TRUE)

  # vectorized drift
  res <- dcswald(rt, response, drift = 1:n, bound = 1.5, ndt = 0.3,
                 version = "simple")
  expect_length(res, n)

  # all parameters vectorized
  res <- dcswald(rt, response, drift = 1:n, bound = seq(1, 2, length.out = n),
                 ndt = seq(0.1, 0.3, length.out = n), version = "simple")
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
  out <- rcswald(n, drift = runif(n, 1, 3), bound = runif(n, 1, 2),
                 ndt = runif(n, 0.1, 0.3), zr = runif(n, 0.3, 0.7))
  expect_equal(nrow(out), n)
})

test_that("pcswald returns probabilities between 0 and 1 for simple version", {
  q <- seq(0.4, 2, length.out = 20)
  p <- pcswald(q, response = 1, drift = 2, bound = 1.5, ndt = 0.3,
               version = "simple")

  expect_length(p, 20)
  expect_true(all(p >= 0, na.rm = TRUE))
  expect_true(all(p <= 1, na.rm = TRUE))
  # CDF should be monotonically increasing
  expect_true(all(diff(p) >= 0))
})

test_that("pcswald returns probabilities between 0 and 1 for crisk version", {
  q <- seq(0.4, 2, length.out = 20)

  # upper boundary
  p_upper <- pcswald(q, response = 1, drift = 2, bound = 1.5, ndt = 0.3,
                     version = "crisk")
  expect_true(all(p_upper >= 0))
  expect_true(all(p_upper <= 1))

  # lower boundary
  p_lower <- pcswald(q, response = 0, drift = 2, bound = 1.5, ndt = 0.3,
                     version = "crisk")
  expect_true(all(p_lower >= 0))
  expect_true(all(p_lower <= 1))
})

test_that("pcswald warns for response=0 in simple version", {
  expect_warning(
    pcswald(1.0, response = 0, drift = 2, bound = 1.5, ndt = 0.3,
            version = "simple"),
    "CDF for response=0 is not well-defined"
  )
})

test_that("pcswald handles log.p and lower.tail arguments", {
  q <- 1.0
  p <- pcswald(q, response = 1, drift = 2, bound = 1.5, ndt = 0.3)

  # log.p = TRUE
  p_log <- pcswald(q, response = 1, drift = 2, bound = 1.5, ndt = 0.3,
                   log.p = TRUE)
  expect_equal(p_log, log(p))

  # lower.tail = FALSE
  p_upper <- pcswald(q, response = 1, drift = 2, bound = 1.5, ndt = 0.3,
                     lower.tail = FALSE)
  expect_equal(p_upper, 1 - p)
})

test_that("qcswald returns valid quantiles for simple version", {
  p <- c(0.1, 0.5, 0.9)
  q <- qcswald(p, response = 1, drift = 2, bound = 1.5, ndt = 0.3,
               version = "simple")

  expect_length(q, 3)
  expect_true(all(q > 0.3))  # all quantiles > ndt
  expect_true(all(diff(q) > 0))  # monotonically increasing
})

test_that("qcswald returns valid quantiles for crisk version", {
  p <- c(0.1, 0.5, 0.9)
  q <- qcswald(p, response = 1, drift = 2, bound = 1.5, ndt = 0.3,
               version = "crisk")

  expect_length(q, 3)
  expect_true(all(q > 0.3))  # all quantiles > ndt
})

test_that("qcswald warns for response=0 in simple version", {
  expect_warning(
    qcswald(0.5, response = 0, drift = 2, bound = 1.5, ndt = 0.3,
            version = "simple"),
    "Quantile for response=0 is not well-defined"
  )
})

test_that("pcswald and qcswald are inverses (round-trip)", {
  p_orig <- c(0.1, 0.3, 0.5, 0.7, 0.9)

  # simple version
  q <- qcswald(p_orig, response = 1, drift = 2, bound = 1.5, ndt = 0.3,
               version = "simple")
  p_back <- pcswald(q, response = 1, drift = 2, bound = 1.5, ndt = 0.3,
                    version = "simple")
  expect_equal(p_orig, p_back, tolerance = 1e-4)

  # crisk version
  q_crisk <- qcswald(p_orig, response = 1, drift = 2, bound = 1.5, ndt = 0.3,
                     version = "crisk")
  p_back_crisk <- pcswald(q_crisk, response = 1, drift = 2, bound = 1.5,
                          ndt = 0.3, version = "crisk")
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
