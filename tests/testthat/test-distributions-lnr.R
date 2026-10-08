test_that("dlnr validates inputs and log scale output", {
  rt <- c(0.5, 0.8)
  response <- c(1, 2)
  m <- c(-1, 0)
  s <- 1
  ndt <- 0.2

  dens <- dlnr(rt, response, m = m, s = s, ndt = ndt)
  log_dens <- dlnr(rt, response, m = m, s = s, ndt = ndt, log = TRUE)

  expect_equal(dens, exp(log_dens))
  expect_error(
    dlnr(rt, c(1, 3), m = m, s = s, ndt = ndt),
    "response must contain integers"
  )
  expect_error(
    dlnr(rt, response, m = m, s = 0, ndt = ndt),
    "finite and positive"
  )
  expect_error(
    dlnr(rt, response, m = m, s = s, ndt = -0.1),
    "finite and non-negative"
  )
})

test_that("rlnr returns valid response times and responses", {
  set.seed(123)
  res <- rlnr(200, m = c(-1, 0, 0.5), s = 1, ndt = 0.2)

  expect_s3_class(res, "data.frame")
  expect_equal(names(res), c("rt", "response"))
  expect_true(all(res$rt > 0.2))
  expect_true(all(res$response %in% 1:3))
})

test_that("plnr and qlnr are internally consistent", {
  m <- c(-1, 0)
  s <- 1
  ndt <- 0.2
  p <- c(0.1, 0.5, 0.9, 0.99)

  q <- qlnr(p, m = m, s = s, ndt = ndt)

  expect_equal(plnr(q, m = m, s = s, ndt = ndt), p, tolerance = 1e-6)
  expect_true(all(diff(plnr(seq(0.21, 3, length.out = 20),
                             m = m, s = s, ndt = ndt)) >= 0))
})

test_that("plnr(response = ) gives defective CDFs that sum to the marginal", {
  m <- c(-1, 0, 0.5)
  s <- c(0.4, 1.2, 0.8)
  ndt <- 0.2
  q <- c(0.3, 0.5, 1, 2, 10)

  defective <- vapply(seq_along(m), function(r) {
    plnr(q, response = r, m = m, s = s, ndt = ndt)
  }, numeric(length(q)))

  expect_equal(rowSums(defective), plnr(q, m = m, s = s, ndt = ndt),
               tolerance = 1e-8)
  # at a q far beyond the race each defective CDF is that response's probability
  far <- vapply(seq_along(m), function(r) {
    plnr(1e6, response = r, m = m, s = s, ndt = ndt)
  }, numeric(1))
  expect_equal(sum(far), 1, tolerance = 1e-8)
  expect_equal(far[1], 0.7531, tolerance = 1e-3)

  expect_equal(
    plnr(q, response = 1, m = m, s = s, ndt = ndt, lower.tail = FALSE),
    1 - defective[, 1], tolerance = 1e-8
  )
  expect_error(
    plnr(q, response = c(1, 2), m = m, s = s, ndt = ndt),
    "response has 2 entries but q has 5"
  )
})

test_that("dlnr integrates to one over the responses", {
  ndt <- 0.2
  for (par in list(list(m = c(-1, 0), s = c(1, 1)),
                   list(m = c(-1, 0, 0.5), s = 1),
                   list(m = c(-1, 0, 0.5), s = c(0.4, 1.2, 0.8)))) {
    total <- sum(vapply(seq_along(par$m), function(r) {
      stats::integrate(
        function(x) dlnr(x, rep(r, length(x)), m = par$m, s = par$s, ndt = ndt),
        ndt, Inf, rel.tol = 1e-10
      )$value
    }, numeric(1)))
    expect_equal(total, 1, tolerance = 1e-8)
  }
})

test_that("rlnr reproduces its own density", {
  withr::local_seed(99)
  m <- c(-1, 0, 0.5)
  s <- c(0.4, 1.2, 0.8)
  ndt <- 0.2
  n <- 2e5
  draws <- rlnr(n, m = m, s = s, ndt = ndt)

  for (r in seq_along(m)) {
    expected <- stats::integrate(
      function(x) dlnr(x, rep(r, length(x)), m = m, s = s, ndt = ndt),
      ndt, Inf, rel.tol = 1e-10
    )$value
    # 4 binomial standard errors
    expect_equal(mean(draws$response == r), expected,
                 tolerance = 4 * sqrt(expected * (1 - expected) / n) / expected)
  }
  probs <- c(0.1, 0.25, 0.5, 0.75, 0.9)
  expect_equal(unname(stats::quantile(draws$rt, probs)),
               qlnr(probs, m = m, s = s, ndt = ndt), tolerance = 0.01)
})

test_that("dlnr is zero outside its support and propagates NA", {
  m <- c(-1, 0)
  expect_equal(dlnr(c(0.1, 0.2, NA), c(1, 1, 1), m = m, s = 1, ndt = 0.2),
               c(0, 0, NA))
  expect_equal(dlnr(c(0.1, 0.2, NA), c(1, 1, 1), m = m, s = 1, ndt = 0.2,
                    log = TRUE),
               c(-Inf, -Inf, NA))
})

test_that("the lnr distribution functions refuse mis-sized arguments", {
  m <- c(-1, 0)
  expect_error(dlnr(c(0.5, 0.8), 1, m = m, s = 1, ndt = 0.2),
               "response has 1 entries but rt has 2")
  expect_error(dlnr(c(0.5, 0.8), c(1, 2), m = m, s = c(1, 1, 1), ndt = 0.2),
               "s has 3 entries but the race has 2 accumulators")
  expect_error(dlnr(c(0.5, 0.8), c(1, 2), m = m, s = 1, ndt = c(0.1, 0.2, 0.3)),
               "ndt has 3 entries but there are 2 response times")
})

test_that("rlnr(0) returns an empty data.frame with typed columns", {
  out <- rlnr(0, m = c(-1, 0), s = 1, ndt = 0.2)
  expect_equal(nrow(out), 0)
  expect_type(out$rt, "double")
  expect_type(out$response, "integer")
})

test_that("plnr keeps the upper tail where the lower tail has rounded to one", {
  m <- c(-1, 0)
  s <- c(1, 1)
  ndt <- 0.2
  q <- c(10, 100, 500, 5000)

  direct <- vapply(q, function(qi) {
    sum(stats::plnorm(qi - ndt, meanlog = m, sdlog = s,
                      lower.tail = FALSE, log.p = TRUE))
  }, numeric(1))

  expect_equal(
    plnr(q, m = m, s = s, ndt = ndt, lower.tail = FALSE, log.p = TRUE),
    direct
  )
})
