library(testthat)

test_that("k2sd works", {
  # Test format
  kappa <- runif(10)
  expect_length(k2sd(kappa), 10)
  expect_type(k2sd(kappa), "double")
  expect_length(k2sd(2), 1)
  expect_type(k2sd(2), "double")

  # Test edge cases
  expect_equal(k2sd(0), Inf)
  expect_equal(k2sd(Inf), 0)

  # Test vector of edge cases
  expect_equal(k2sd(c(0, Inf)), c(Inf, 0))

  # Test known values (compared to pre-computed results)
  expect_equal(k2sd(1), 1.270088, tolerance = 1e-6)
  expect_equal(k2sd(10), 0.3248638, tolerance = 1e-6)

  # Test NA handling
  expect_true(is.na(k2sd(NA)))
  expect_true(all(is.na(k2sd(c(1, NA, 3))[2])))

  # Test invalid inputs
  expect_error(k2sd("a"))
  expect_error(k2sd(NULL))
})

test_that("c parameter conversions work", {
  # Test basic conversion
  c_sqrtexp <- 4
  kappa <- 3
  c_bessel <- c_sqrtexp2bessel(c_sqrtexp, kappa)
  expect_equal(c_bessel2sqrtexp(c_bessel, kappa), c_sqrtexp)

  # Test vector inputs
  c_vec <- c(1, 2, 3)
  kappa_vec <- c(2, 3, 4)
  c_bessel_vec <- c_sqrtexp2bessel(c_vec, kappa_vec)
  expect_equal(c_bessel2sqrtexp(c_bessel_vec, kappa_vec), c_vec)

  # Test error handling
  expect_error(c_sqrtexp2bessel(-1, 2), "c must be non-negative")
  expect_error(c_bessel2sqrtexp(1, -2), "kappa must be non-negative")
})

test_that("identity is a no-op and round-trips", {
  x <- c(-2, -0.5, 0, 1.25, 3)
  eta <- link_transform(x, "identity", inverse = FALSE)
  back <- link_transform(eta, "identity", inverse = TRUE)
  expect_identical(eta, x)
  expect_identical(back, x)
})

test_that("log and inverse log round-trip on positive values", {
  x <- c(0.1, 0.5, 1, 2, 10)
  eta <- link_transform(x, "log", inverse = FALSE)
  back <- link_transform(eta, "log", inverse = TRUE)
  expect_equal(back, x, tolerance = 1e-12)
})

test_that("log1p and inverse expm1 round-trip for x > -1", {
  x <- c(-0.5, -0.1, 0, 0.3, 5)
  eta <- link_transform(x, "log1p", inverse = FALSE)
  back <- link_transform(eta, "log1p", inverse = TRUE)
  expect_equal(back, x, tolerance = 1e-12)
})

test_that("logm1 (brms::logm1) round-trip for x > 1", {
  testthat::skip_if_not_installed("brms")
  x <- c(1.1, 2, 5, 10)
  eta <- link_transform(x, "logm1", inverse = FALSE)
  back <- link_transform(eta, "logm1", inverse = TRUE)
  expect_equal(back, x, tolerance = 1e-12)
})

test_that("inverse link is its own inverse (x != 0)", {
  x <- c(-3, -0.5, 0.2, 4)
  eta <- link_transform(x, "inverse", inverse = FALSE)
  back <- link_transform(eta, "inverse", inverse = TRUE)
  expect_equal(back, x, tolerance = 1e-12)
})

test_that("sqrt link round-trip for x >= 0", {
  x <- c(0, 0.01, 1, 2, 9)
  eta <- link_transform(x, "sqrt", inverse = FALSE)
  back <- link_transform(eta, "sqrt", inverse = TRUE)
  expect_equal(back, x, tolerance = 1e-12)
})

test_that("logit round-trip for probabilities in (0,1)", {
  p <- c(0.1, 0.25, 0.5, 0.75, 0.9)
  eta <- link_transform(p, "logit", inverse = FALSE)
  back <- link_transform(eta, "logit", inverse = TRUE)
  expect_equal(back, p, tolerance = 1e-12)
})

test_that("probit round-trip for probabilities in (0,1)", {
  p <- c(0.01, 0.2, 0.5, 0.8, 0.99)
  eta <- link_transform(p, "probit", inverse = FALSE)
  back <- link_transform(eta, "probit", inverse = TRUE)
  expect_equal(back, p, tolerance = 1e-12)
})

test_that("tan_half round-trip on a safe interval (-pi, pi)", {
  x <- c(-2, -1, 0, 1, 2)  # well within (-pi, pi)
  eta <- link_transform(x, "tan_half", inverse = FALSE)   # tan(x/2)
  back <- link_transform(eta, "tan_half", inverse = TRUE) # 2*atan(eta)
  expect_equal(back, x, tolerance = 1e-12)
})


test_that("cloglog round-trip for probabilities in (0,1)", {
  p <- seq(0.1,0.9,by = 0.1)
  eta <- link_transform(p, "cloglog", inverse = FALSE)    # should be log(-log(1 - p))
  back <- link_transform(eta, "cloglog", inverse = TRUE)  # should be 1 - exp(-exp(eta))
  # This expectation is correct; it will fail with current implementation.
  expect_equal(back, p, tolerance = 1e-9)
})

test_that("loglog round-trip for probabilities in (0,1)", {
  p <- seq(0.1,0.9,by = 0.1)
  eta <- link_transform(p, "loglog", inverse = FALSE)     # log(-log(p))
  back <- link_transform(eta, "loglog", inverse = TRUE)   # exp(-exp(eta))
  expect_equal(back, p, tolerance = 1e-12)
})

test_that("loglog inverse maps reals to (0,1)", {
  eta <- c(-3, -2, 0, 1, 3)
  p <- link_transform(eta, "loglog", inverse = TRUE)
  expect_true(all(is.finite(p)))
  expect_true(all(p > 0 & p < 1))
})

test_that("loglog is monotone decreasing in p", {
  p <- c(0.1, 0.2, 0.4, 0.8)  # increasing p
  eta <- link_transform(p, "loglog", inverse = FALSE)
  # as p increases, eta decreases strictly
  expect_true(all(diff(eta) < 0))
})

test_that("loglog vectorization and NA propagation", {
  p <- c(0.2, NA_real_, 0.7)
  eta <- link_transform(p, "loglog", inverse = FALSE)
  back <- link_transform(eta, "loglog", inverse = TRUE)
  expect_true(is.na(eta[2]))
  expect_true(is.na(back[2]))
  expect_equal(back[c(1,3)], p[c(1,3)], tolerance = 1e-12)
})

test_that("vectorization works and NA positions are preserved", {
  p <- c(0.2, NA_real_, 0.8)
  eta <- link_transform(p, "logit", inverse = FALSE)
  back <- link_transform(eta, "logit", inverse = TRUE)
  expect_true(is.na(eta[2]))
  expect_true(is.na(back[2]))
  expect_equal(back[c(1,3)], p[c(1,3)], tolerance = 1e-12)
})

test_that("unknown link errors clearly", {
  expect_error(link_transform(1:3, "not_a_link"))
})

test_that("non-numeric values error", {
  # If you rely on glue(), also ensure glue is installed; otherwise any error is fine.
  expect_error(link_transform(c("a","b"), "log"))
})

test_that("NULL link is treated as identity", {
  x <- c(-2, -0.5, 0, 1.25, 3)
  eta <- link_transform(x, NULL, inverse = FALSE)
  back <- link_transform(eta, NULL, inverse = TRUE)
  expect_identical(eta, x)
  expect_identical(back, x)
})
