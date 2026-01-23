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
    c("mean_rt_upper", "mean_rt_lower", "var_rt_upper",
      "var_rt_lower", "n_upper", "n_trials") %in% names(res)
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
  # drift must be positive
  expect_error(
    dezdm(mean_rt = 0.5, var_rt = 0.02, n_upper = 80, n_trials = 100,
          drift = -1, bound = 1.5, ndt = 0.3),
    "drift must be positive"
  )

  # bound must be positive
  expect_error(
    dezdm(mean_rt = 0.5, var_rt = 0.02, n_upper = 80, n_trials = 100,
          drift = 2, bound = -1, ndt = 0.3),
    "bound must be positive"
  )

  # ndt must be positive
  expect_error(
    dezdm(mean_rt = 0.5, var_rt = 0.02, n_upper = 80, n_trials = 100,
          drift = 2, bound = 1.5, ndt = -0.1),
    "ndt must be positive"
  )

  # n_upper cannot exceed n_trials
  expect_error(
    dezdm(mean_rt = 0.5, var_rt = 0.02, n_upper = 150, n_trials = 100,
          drift = 2, bound = 1.5, ndt = 0.3),
    "n_upper cannot exceed n_trials"
  )
  
  # n_trials must be larger than 2
  expect_error(
    dezdm(mean_rt = 0.5, var_rt = 0.02, n_upper = 1, n_trials = 1,
          drift = 2, bound = 1.5, ndt = 0.3),
    "n_trials must be larger than 2"
  )
  
  expect_error(
    dezdm(mean_rt = 0.5, var_rt = 0.02, n_upper = 2, n_trials = 2,
          drift = 2, bound = 1.5, ndt = 0.3),
    "n_trials must be larger than 2"
  )

  # version must be valid
  expect_error(
    dezdm(mean_rt = 0.5, var_rt = 0.02, n_upper = 80, n_trials = 100,
          drift = 2, bound = 1.5, ndt = 0.3, version = "5par"),
    "version must be either"
  )

  # 4par requires length-2 vectors
  expect_error(
    dezdm(mean_rt = 0.5, var_rt = 0.02, n_upper = 80, n_trials = 100,
          drift = 2, bound = 1.5, ndt = 0.3, version = "4par"),
    "mean_rt must be length 2"
  )

  # zr must be between 0 and 1 for 4par
  expect_error(
    dezdm(mean_rt = c(0.5, 0.6), var_rt = c(0.02, 0.03),
          n_upper = 80, n_trials = 100,
          drift = 2, bound = 1.5, ndt = 0.3, zr = 1.5, version = "4par"),
    "zr must be between 0 and 1"
  )
})

test_that("rezdm validates parameters correctly", {
  # drift must be positive
  expect_error(
    rezdm(n = 10, n_trials = 100, drift = -1, bound = 1.5, ndt = 0.3),
    "drift must be positive"
  )

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
  ll1 <- dezdm(mean_rt = 0.5, var_rt = 0.02, n_upper = 80, n_trials = 100,
               drift = 2, bound = 1.5, ndt = 0.3)
  ll2 <- dezdm(mean_rt = 0.55, var_rt = 0.025, n_upper = 75, n_trials = 100,
               drift = 2, bound = 1.5, ndt = 0.3)

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
