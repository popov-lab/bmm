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
  expect_identical(.sdm_log_z(NA_real_, 2, "laplace"), NA_real_)
  expect_identical(dsdm(numeric(0), 0, 3, 2), numeric(0))
  expect_identical(dsdm(0.3, 0, NA, 2), NA_real_)
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

test_that("kappa = 0 gives the uniform density where the cell has a limit there", {
  expect_equal(dsdm(c(-2, 0.5), 0, 2, 0, activation = "laplace", scaling = "density"), rep(1 / (2 * pi), 2))
  expect_equal(dsdm(c(-2, 0.5), 0, 2, 0, activation = "laplace"), rep(1 / (2 * pi), 2))
  expect_equal(dsdm(c(-2, 0.5), 0, 2, 0), rep(1 / (2 * pi), 2))
  expect_equal(
    suppressWarnings(dsdm(c(-2, 0.5), 0, 3, 0, parametrization = "bessel")),
    rep(1 / (2 * pi), 2)
  )
})

# Discrete response options (#558, K2) -----------------------------------------

test_that("dsdm() with n_options times the step sums to 1 over the grid for every cell", {
  for (K in c(4L, 90L, 360L)) {
    x <- 2 * pi * (seq_len(K) - K / 2) / K
    for (i in seq_len(nrow(sdm_cells))) {
      cell <- sdm_cells[i, ]
      for (par in list(c(0, 2, 3), c(0.3, 40, 0.5), c(-1.7, 0.2, 80))) {
        total <- sum(dsdm(x, par[1], par[2], par[3],
          n_options = K, activation = cell$activation, scaling = cell$scaling
        )) * 2 * pi / K
        expect_lt(abs(total - 1), 1e-12,
          label = glue("|sum - 1| at K {K}, {cell$activation}, {cell$scaling}, mu {par[1]}")
        )
      }
    }
  }
})

test_that("dsdm() with n_options is the softmax over the options on the density scale", {
  K <- 12
  x <- 2 * pi * (-5:6) / K
  a <- 7 * exp(-3 * (1 - cos(x - 0.2)))
  by_hand <- a - log(sum(exp(a))) - log(2 * pi / K)
  expect_equal(dsdm(x, 0.2, 7, 3, n_options = K, scaling = "peak", log = TRUE), by_hand, tolerance = 1e-13)
  # a response is taken at its nearest option
  expect_equal(dsdm(x + 0.2, 0.2, 7, 3, n_options = K, scaling = "peak"), exp(by_hand), tolerance = 1e-13)
})

test_that("dsdm() with options is the softmax over the offered options", {
  options <- deg2rad(c(0, 45, 180, -90))
  values <- dsdm(options, 0, 4, 2, n_options = 360, options = options, log = TRUE)
  expect_lt(abs(sum(exp(values)) * 2 * pi / 360 - 1), 1e-13)
  a <- .sdm_log_cstar(log(4), 2, "vonmises", "density")
  a <- exp(a - 2 * 2 * sin(options / 2)^2)
  expect_equal(values, a - log(sum(exp(a))) - log(2 * pi / 360), tolerance = 1e-13)
  # a response the trial did not offer has no probability
  expect_identical(dsdm(deg2rad(90), 0, 4, 2, n_options = 360, options = options), 0)
  expect_identical(dsdm(deg2rad(90), 0, 4, 2, n_options = 360, options = options, log = TRUE), -Inf)
})

test_that("dsdm() refuses options without n_options and a bad n_options", {
  expect_error(dsdm(0, options = c(0, 1)), "need `n_options`")
  expect_error(dsdm(0, n_options = 7), "even whole number of at least 4")
  expect_error(rsdm(3, n_options = 5.5), "even whole number of at least 4")
})

# The discrete von Mises SDM is the continuous one where the activation is
# wider than a few steps (k2_logs/discrete_vs_continuous.md): the error in log
# Z falls as exp(-2.2 r^2) in r = w / h, w = 1 / sqrt(max(c*, 1) kappa), and
# is below 1e-10 from r = 2.25 on. Laplace never converges: its kink at the
# target makes the step error of the grid sum O(h^2 c* kappa)
test_that("the discrete density approaches the continuous one with the measured curve", {
  # on both grids
  x <- 2 * pi * c(-44, -10, 0, 1, 7, 45) / 90
  cells <- expand.grid(lc = c(-1, 0, 1, 2, 3, 5), kappa = c(0.5, 2, 10, 50, 300), K = c(90, 360))
  for (i in seq_len(nrow(cells))) {
    cell <- cells[i, ]
    log_cstar <- .sdm_log_cstar(cell$lc, cell$kappa, "vonmises", "density")
    r <- 1 / sqrt(max(exp(log_cstar), 1) * cell$kappa) / (2 * pi / cell$K)
    discrete <- dsdm(x, 0, exp(cell$lc), cell$kappa, n_options = cell$K, log = TRUE)
    continuous <- dsdm(x, 0, exp(cell$lc), cell$kappa, log = TRUE)
    if (r < 0.75) next
    bound <- if (r >= 2.25) 1e-10 else max(1e-10, exp(-4 * r^2))
    expect_lt(max(abs(discrete - continuous)), bound,
      label = glue("|discrete - continuous| at log c {cell$lc}, kappa {cell$kappa}, K {cell$K}, r {signif(r, 3)}")
    )
  }
  laplace <- abs(dsdm(0, 0, 50, 20, n_options = 360, activation = "laplace", log = TRUE) -
    dsdm(0, 0, 50, 20, activation = "laplace", log = TRUE))
  expect_gt(laplace, 1e-3)
})

test_that("rsdm() with n_options or options draws option centres with their probabilities", {
  withr::local_seed(558)
  draws <- rsdm(40000, 0, 3, 2, n_options = 8)
  steps <- draws / (2 * pi / 8)
  expect_lt(max(abs(steps - round(steps))), 1e-12)
  expect_true(all(draws > -pi & draws <= pi))
  x <- 2 * pi * (-3:4) / 8
  expected <- dsdm(x, 0, 3, 2, n_options = 8) * 2 * pi / 8
  observed <- tabulate(round(steps) + 4, 8) / 40000
  expect_lt(max(abs(observed - expected)), 0.01)

  options <- deg2rad(c(0, 120, -120))
  draws <- rsdm(30000, 0, 3, 2, n_options = 360, options = options)
  expect_setequal(round(rad2deg(unique(draws))), c(0, 120, -120))
  expected <- dsdm(options, 0, 3, 2, n_options = 360, options = options) * 2 * pi / 360
  observed <- vapply(options, function(o) mean(abs(draws - o) < 1e-9), numeric(1))
  expect_lt(max(abs(observed - expected)), 0.01)
})

# The softmax of d f under Gumbel noise, with the von Mises similarity scaled
# to 1 at the target and 0 opposite it, f = (s - s(pi)) / (1 - s(pi)), is
# sdm(scaling = "range") with c = d: the constant cancels in the softmax.
# dtcc() is on the TCC branch until TCC-K3 brings both together
test_that("dsdm() with n_options is tcc(noise = 'gumbel') with a von Mises similarity", {
  skip_if_not(exists("dtcc", envir = asNamespace("bmm"), inherits = FALSE), "dtcc() not on this branch")
  dtcc <- get("dtcc", envir = asNamespace("bmm"))
  x <- 2 * pi * (-179:180) / 360
  for (kappa in c(0.5, 3, 30)) {
    s <- exp(-2 * kappa * sin(2 * pi * (0:180) / 360 / 2)^2)
    f <- (s - exp(-2 * kappa)) / -expm1(-2 * kappa)
    for (d in c(0.5, 4, 20)) {
      expect_equal(
        dtcc(x, d, similarity = f, noise = "gumbel", log = TRUE) - log(2 * pi / 360),
        dsdm(x, 0, d, kappa, n_options = 360, scaling = "range", log = TRUE),
        tolerance = 1e-10
      )
    }
  }
})

# K1's Laplace density scaling took its length from kappa alone
test_that("dsdm() gives one value per c with a single kappa in every cell", {
  for (i in seq_len(nrow(sdm_cells))) {
    cell <- sdm_cells[i, ]
    one_by_one <- vapply(c(5, 6, 7), function(c) {
      dsdm(0.3, 0, c, 2, activation = cell$activation, scaling = cell$scaling)
    }, numeric(1))
    expect_equal(dsdm(0.3, 0, c(5, 6, 7), 2, activation = cell$activation, scaling = cell$scaling),
                 one_by_one, label = glue("{cell$activation}, {cell$scaling}"))
  }
})
