test_that(".sdm_log_z() matches integrate() over the #615 cells for both activations", {
  skip_on_cran() # about 1 s of integrate()

  lc <- c(-1, 0.5, 2, 3.5, 5)
  kappa <- c(0.5, 2, 10, 30, 100, 300)
  cells <- rbind(
    expand.grid(lc = lc, kappa = kappa),
    data.frame(lc = c(6, 7, 8, 10, 12, -1, 12), kappa = c(300, 30, 5, 1, 0.5, 3000, 3000)),
    expand.grid(lc = c(-1, 0, 2, 5), kappa = c(0.05, 0.1))
  )
  log_peak <- cells$lc + 0.5 * log(cells$kappa / (2 * pi))
  for (activation in c("vonmises", "laplace")) {
    reference <- mapply(sdm_log_z_integrate, log_peak, cells$kappa,
      MoreArgs = list(activation = activation)
    )
    error <- abs(.sdm_log_z(log_peak, cells$kappa, activation) - reference)
    expect_lt(max(error), 1e-7, label = glue("max |error| in log Z, {activation}"))
    expect_lt(max(error[cells$kappa <= 30 & cells$lc <= 5]), 1e-10,
      label = glue("max |error| in log Z at kappa <= 30, {activation}")
    )
  }
})

test_that(".sdm_log_z() keeps NA cells and the constant activation of c = 0", {
  expect_identical(.sdm_log_z(c(NA, 1), 2, "laplace")[1], NA_real_)
  expect_equal(.sdm_log_z(-Inf, c(0.5, 30), "vonmises"), rep(log(2 * pi), 2))
  expect_equal(.sdm_log_z(-Inf, c(0.5, 30), "laplace"), rep(log(2 * pi), 2))
})

test_that("the dsdm() density of every cell integrates to 1", {
  pars <- data.frame(c = c(3, 20, 0.5, 40), kappa = c(3, 30, 0.2, 8))
  for (i in seq_len(nrow(sdm_cells))) {
    cell <- sdm_cells[i, ]
    for (j in seq_len(nrow(pars))) {
      # split at the kink of the Laplace activation
      mass <- sum(vapply(list(c(-pi, 0), c(0, pi)), function(range) {
        stats::integrate(dsdm, range[1], range[2],
          c = pars$c[j], kappa = pars$kappa[j],
          activation = cell$activation, scaling = cell$scaling,
          rel.tol = 1e-12, subdivisions = 1000L
        )$value
      }, numeric(1)))
      expect_equal(mass, 1, tolerance = 1e-9,
        label = glue("mass of {cell$activation}, {cell$scaling}, c {pars$c[j]}, kappa {pars$kappa[j]}")
      )
    }
  }
})

test_that("two scalings give the same density once c is converted between them", {
  x <- c(-pi, -1, -0.1, 0, 0.05, 0.7, 2.5)
  for (activation in c("vonmises", "laplace")) {
    for (c in c(0.7, 4, 30)) {
      for (kappa in c(0.3, 3, 25)) {
        for (s1 in c("density", "peak", "range")) {
          cstar <- sdm_cstar_by_hand(c, kappa, activation, s1)
          for (s2 in setdiff(c("density", "peak", "range"), s1)) {
            c2 <- cstar / sdm_cstar_by_hand(1, kappa, activation, s2)
            d1 <- dsdm(x, 0, c, kappa, log = TRUE, activation = activation, scaling = s1)
            d2 <- dsdm(x, 0, c2, kappa, log = TRUE, activation = activation, scaling = s2)
            expect_lt(max(abs(d1 - d2)), 1e-12,
              label = glue("{activation}: {s1} c {c} vs {s2}, kappa {kappa}")
            )
          }
        }
      }
    }
  }
})

test_that("the activation of dsdm() is the one written out in ?sdm", {
  x <- c(-3, -0.4, 0, 0.2, 1.5, 3)
  # c* exp(log s(x)) - log Z, with Z from integrate()
  for (i in seq_len(nrow(sdm_cells))) {
    cell <- sdm_cells[i, ]
    cstar <- sdm_cstar_by_hand(5, 4, cell$activation, cell$scaling)
    s <- if (cell$activation == "vonmises") exp(4 * (cos(x) - 1)) else exp(-4 * abs(x))
    expected <- cstar * s - sdm_log_z_integrate(log(cstar), 4, cell$activation)
    expect_equal(dsdm(x, 0, 5, 4, log = TRUE, activation = cell$activation, scaling = cell$scaling),
      expected,
      tolerance = 1e-12
    )
  }
})

test_that("dsdm() gives the log density where log Z passes 709", {
  # stopped with "non-finite function value" before #558; log_lik() and so
  # loo() call dsdm(); log c 5 with kappa 200 is in the #615 grid
  expected <- exp(5) * sqrt(200 / (2 * pi)) * exp(200 * (cos(0.01) - 1)) -
    sdm_log_z_integrate(5 + 0.5 * log(200 / (2 * pi)), 200, "vonmises")
  expect_equal(dsdm(0.01, 0, exp(5), 200, log = TRUE), expected, tolerance = 1e-12)
  expect_true(is.finite(dsdm(0.01, 0, exp(12), 3000, log = TRUE)))
})

test_that("psdm(), qsdm() and rsdm() follow the density of the Laplace cells", {
  withr::local_seed(558)
  q <- c(-2, -0.3, -0.01, 0.02, 0.5, 2.8)
  for (scaling in c("density", "peak", "range")) {
    lower <- vapply(q, function(b) {
      sum(vapply(list(c(-pi, min(b, 0)), c(min(b, 0), b)), function(r) {
        if (r[2] <= r[1]) return(0)
        stats::integrate(dsdm, r[1], r[2], c = 6, kappa = 5, activation = "laplace",
          scaling = scaling, rel.tol = 1e-12)$value
      }, numeric(1)))
    }, numeric(1))
    p <- psdm(q, c = 6, kappa = 5, activation = "laplace", scaling = scaling)
    expect_equal(p, lower, tolerance = 1e-9)
    expect_equal(qsdm(p, c = 6, kappa = 5, activation = "laplace", scaling = scaling), q, tolerance = 1e-8)
  }
  # 20000 draws: the mean of cos(x) has an SD of about 0.005
  draws <- rsdm(20000, c = 6, kappa = 5, activation = "laplace")
  mean_cos <- stats::integrate(function(x) cos(x) * dsdm(x, c = 6, kappa = 5, activation = "laplace"),
    0, pi, rel.tol = 1e-10)$value * 2
  expect_lt(abs(mean(cos(draws)) - mean_cos), 0.02)
})

test_that("the d/p/q/r functions refuse an unknown activation or scaling", {
  expect_error(dsdm(0.5, activation = "gauss"), "should be one of")
  expect_error(psdm(0.5, scaling = "height"), "scaling must be one of 'density', 'peak' or 'range'", fixed = TRUE)
  expect_error(rsdm(2, activation = "laplace", parametrization = "sqrtexp"),
    "`parametrization` is deprecated and applies only to the von Mises activation", fixed = TRUE
  )
  expect_error(qsdm(0.5, scaling = "peak", parametrization = "bessel"),
    "`parametrization` is deprecated and applies only to the von Mises activation", fixed = TRUE
  )
})

test_that("parametrization is deprecated and keeps the values it gave", {
  expect_warning(
    old <- dsdm(0.3, 0, 4, 3, parametrization = "sqrtexp"),
    '`parametrization = "sqrtexp"` is deprecated since bmm 1.5.0', fixed = TRUE
  )
  expect_identical(old, dsdm(0.3, 0, 4, 3, scaling = "density"))
  expect_warning(
    bessel <- dsdm(0.3, 0, 4, 3, parametrization = "bessel"),
    "c_bessel2sqrtexp(c, kappa)", fixed = TRUE
  )
  expect_equal(bessel, dsdm(0.3, 0, c_bessel2sqrtexp(4, 3), 3))
})
