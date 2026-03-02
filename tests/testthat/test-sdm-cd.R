test_that("sdm simple CD constructor works", {
  m <- sdm(response = "resp", probe = "probe", target = "target",
           task = "cd")
  expect_s3_class(m, "bmmodel")
  expect_s3_class(m, "change_detection")
  expect_s3_class(m, "sdm")
  expect_s3_class(m, "sdm_simple")
  expect_s3_class(m, "sdm_simple_cd")
  expect_false(inherits(m, "circular"))
  expect_equal(m$resp_vars$response, "resp")
  expect_true(m$void_mu)
  expect_equal(m$fixed_parameters$beta, 0)
})

test_that("sdm DE constructor still works", {
  m <- sdm(resp_error = "y")
  expect_s3_class(m, "circular")
  expect_s3_class(m, "sdm_simple_de")
  expect_false(inherits(m, "change_detection"))
})

test_that("sdm CD validates required arguments", {
  expect_error(sdm(task = "cd"), "response")
  expect_error(sdm(response = "r", task = "cd"), "probe")
  expect_error(sdm(response = "r", probe = "p", task = "cd"), "target")
})

test_that("SDM normalization is approximately translation-invariant in mu", {
  # exact invariance holds mathematically; small differences are trapezoidal
  # quadrature discretization artifacts from grid/peak alignment
  x_grid <- seq(-pi, pi, length.out = 101)
  dx <- x_grid[2] - x_grid[1]
  for (c_val in c(1, 4, 8)) {
    for (kappa_val in c(1, 5, 20)) {
      Z0 <- sum(exp(bmm:::.dsdm_numer_sqrtexp(x_grid, mu = 0,
                                                c_val, kappa_val, log = TRUE))) * dx
      Z1 <- sum(exp(bmm:::.dsdm_numer_sqrtexp(x_grid, mu = pi / 3,
                                                c_val, kappa_val, log = TRUE))) * dx
      expect_equal(Z0, Z1, tolerance = 1e-3,
                   label = sprintf("Z(mu=0) vs Z(mu=pi/3) for c=%g, kappa=%g",
                                   c_val, kappa_val))
    }
  }
})

test_that("dsdm_cd returns valid probabilities", {
  p <- dsdm_cd(1, probe = 0, c = 4, kappa = 3)
  expect_true(p >= 0 && p <= 1)
})

test_that("dsdm_cd: same probe has lower P(change)", {
  p_same <- dsdm_cd(1, probe = 0, c = 4, kappa = 3)
  p_diff <- dsdm_cd(1, probe = pi, c = 4, kappa = 3)
  expect_lt(p_same, p_diff)
})

test_that("rsdm_cd returns binary values", {
  set.seed(42)
  r <- rsdm_cd(100, probe = 1.0, c = 4, kappa = 3)
  expect_true(all(r %in% c(0, 1)))
})

test_that("dsdm_cd with 51 points matches high-resolution reference", {
  # use 201-point quadrature as reference
  ref_dsdm_cd <- function(response, probe, c, kappa, mu = 0, beta = 0) {
    x_grid <- seq(-pi, pi, length.out = 201)
    ln_ret <- bmm:::.dsdm_numer_sqrtexp(x_grid, mu = mu, c, kappa, log = TRUE)
    ln_same <- bmm:::.dsdm_numer_sqrtexp(x_grid, mu = probe + mu, c, kappa,
                                          log = TRUE)
    exp_ln_ret <- exp(ln_ret)
    llr <- ln_ret - ln_same
    w <- stats::plogis(5 * (llr - beta))
    p <- sum(w * exp_ln_ret) / sum(exp_ln_ret)
    max(min(p, 1 - 1e-10), 1e-10)
  }

  params <- expand.grid(
    c = c(1, 4, 8),
    kappa = c(1, 5, 20),
    probe = c(0, pi / 4, pi)
  )
  for (i in seq_len(nrow(params))) {
    p <- params[i, ]
    val <- dsdm_cd(1, probe = p$probe, c = p$c, kappa = p$kappa)
    ref <- ref_dsdm_cd(1, probe = p$probe, c = p$c, kappa = p$kappa)
    expect_equal(val, ref, tolerance = 0.02,
                 label = sprintf("c=%g, kappa=%g, probe=%g", p$c, p$kappa,
                                 p$probe))
  }
})

test_that("sdm simple CD pipeline runs with mock backend", {
  dat <- data.frame(
    resp = c(rep(0, 50), rep(1, 50)),
    probe = runif(100, -pi, pi),
    target = runif(100, -pi, pi)
  )

  m <- sdm(response = "resp", probe = "probe", target = "target",
           task = "cd")
  f <- bmf(c ~ 1, kappa ~ 1)

  mock_fit <- bmm(f, dat, m, backend = "mock", mock_fit = 1, rename = FALSE)
  expect_equal(mock_fit$fit, 1)
  expect_type(mock_fit$bmm, "list")
})
