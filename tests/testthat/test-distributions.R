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

test_that("dddm handles variability parameters", {
  res_no_var <- dddm(0.5, 1, drift = 2, bound = 1.5, ndt = 0.3, 
                     sdrift = 0, sndt = 0, szr = 0)
  res_with_var <- dddm(0.5, 1, drift = 2, bound = 1.5, ndt = 0.3, 
                       sdrift = 0.5, sndt = 0.05, szr = 0.1)
  expect_type(res_with_var, "double")
  expect_false(identical(res_no_var, res_with_var))
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
    "Invalid responses passed"
  )
})

test_that("dddm rejects mismatched rt and response lengths", {
  expect_error(
    dddm(c(0.5, 0.6), 1, drift = 2, bound = 1.5, ndt = 0.3),
    "Different number of rts and responses"
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

test_that("rddm handles variability parameters", {
  res <- rddm(10, drift = 2, bound = 1.5, ndt = 0.3,
              sdrift = 0.5, sndt = 0.05, szr = 0.1)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 10)
})

test_that("rddm starting point bias affects response proportions", {
  # Bias toward upper boundary
  res_upper <- rddm(1000, drift = 0, bound = 1.5, ndt = 0.3, zr = 0.7)
  expect_gt(mean(res_upper$response == 1), 0.6)
  
  # Bias toward lower boundary
  res_lower <- rddm(1000, drift = 0, bound = 1.5, ndt = 0.3, zr = 0.3)
  expect_gt(mean(res_lower$response == 0), 0.6)
})
