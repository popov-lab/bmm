n_rows <- function(y) length(y)

ld_mixture <- function(y, mu_nt = numeric(0), kappa = 5, w_mem = 1, w_guess = 0) {
  n <- n_rows(y)
  cosd <- cbind(cos(y), vapply(mu_nt, function(m) cos(y - m), numeric(n)))
  logw <- matrix(log(w_mem), nrow = n, ncol = length(w_mem), byrow = TRUE)
  .circmix_ld(cosd, logw, rep(log(w_guess), n), rep(kappa, n))
}

# the table is built for a relative error near 2e-11, so a looser bound would
# miss a coarser table or a broken slope
test_that(".circmix_kappa() inverts .circmix_J() over the tabulated range", {
  kappa <- exp(seq(log(0.02), log(4e4), length.out = 400))
  expect_lt(max(abs(.circmix_kappa(.circmix_J(kappa)) / kappa - 1)), 5e-11)
})

test_that(".circmix_kappa() uses the exact asymptotics outside the table", {
  # the neglected term is O(J^2), 3e-12 just below the table
  J <- c(1e-12, 1e-9, 1e-7, 9.9e-6)
  expect_lt(max(abs(.circmix_J(.circmix_kappa(J)) / J - 1)), 1e-11)
  expect_equal(.circmix_kappa(1e7), 1e7 + 0.5, tolerance = 1e-12)
})

test_that(".circmix_kappa() is continuous at the edges of the table", {
  tab <- .circmix_kappa_table()
  edges <- exp(tab$logJ_min + c(0, (length(tab$logkappa) - 1) * tab$dlogJ))
  for (edge in edges) {
    across <- .circmix_kappa(edge * c(1 - 1e-12, 1 + 1e-12))
    expect_lt(abs(across[2] / across[1] - 1), 1e-10)
  }
})

test_that(".circmix_J() reproduces the limits of kappa * I1/I0", {
  expect_equal(.circmix_J(1e-4), 1e-8 / 2, tolerance = 1e-6)
  expect_equal(.circmix_J(1e4), 1e4 - 0.5, tolerance = 1e-6)
  expect_equal(.circmix_J(0), 0)
  # the series that takes over above kappa 1e3, as in Stan's circmix_J(), has
  # no step there that its gradient would see; three terms left 2e-13
  kappa <- 1e3 * (1 + 1e-12)
  expect_equal(
    .circmix_J(kappa),
    kappa * besselI(kappa, 1, expon.scaled = TRUE) / besselI(kappa, 0, expon.scaled = TRUE),
    tolerance = 1e-14
  )
})

test_that(".circmix_ld() is a normalised density on the circle", {
  total <- stats::integrate(
    function(y) exp(ld_mixture(y, mu_nt = 1.2, kappa = 8, w_mem = c(0.5, 0.3), w_guess = 0.2)),
    -pi, pi
  )
  expect_equal(total$value, 1, tolerance = 1e-6)
})

test_that(".circmix_ld() reduces to a von Mises with one component and no guessing", {
  y <- seq(-pi, pi, length.out = 25)
  expect_equal(
    ld_mixture(y, kappa = 12, w_mem = 1, w_guess = 0),
    brms::dvon_mises(y, mu = 0, kappa = 12, log = TRUE),
    tolerance = 1e-10
  )
})

test_that(".circmix_het_ld() agrees with .circmix_ld() when the components share kappa", {
  y <- seq(-pi, pi, length.out = 25)
  cosd <- cbind(cos(y), cos(y - 1.2))
  logw <- matrix(log(c(0.5, 0.3)), nrow = length(y), ncol = 2, byrow = TRUE)
  guess <- rep(log(0.2), length(y))
  expect_equal(
    .circmix_het_ld(cosd, logw, guess, matrix(9, nrow = length(y), ncol = 2)),
    .circmix_ld(cosd, logw, guess, rep(9, length(y)))
  )
})

test_that(".circmix_vp_ld() returns the constant-precision density when tau is zero", {
  y <- seq(-pi, pi, length.out = 25)
  cosd <- matrix(cos(y))
  logw <- matrix(log(0.8), nrow = length(y))
  guess <- rep(log(0.2), length(y))
  expect_equal(
    .circmix_vp_ld(cosd, logw, guess, rep(10, length(y)), rep(0, length(y))),
    .circmix_ld(cosd, logw, guess, rep(10, length(y)))
  )
  # kappa = 0 gives J = 0: no gamma to integrate over
  expect_equal(
    .circmix_vp_ld(cosd, logw, guess, rep(0, length(y)), rep(0.5, length(y))),
    .circmix_ld(cosd, logw, guess, rep(0, length(y)))
  )
})

# tau -> Inf at fixed kappa leaves the mean of J at J(kappa) but sends the
# shape to 0, so the mass of J collapses onto 0 and memory becomes uniform.
# vM(kappa) is the opposite limit, tau -> 0.
test_that(".circmix_vp_ld() has uniform memory components when tau is infinite", {
  y <- c(0, 1, pi)
  n <- length(y)
  for (p_mem in c(1, 0.8)) {
    expect_equal(
      .circmix_vp_ld(matrix(cos(y)), matrix(log(p_mem), n), rep(log1p(-p_mem), n), rep(5, n), rep(Inf, n)),
      rep(-log(2 * pi), n),
      label = paste("p_mem", p_mem)
    )
  }
  # J underflows to 0 at kappa 1e-170, a von Mises that is uniform to 1e-170
  expect_equal(
    .circmix_vp_ld(matrix(cos(y)), matrix(0, n), rep(-Inf, n), rep(1e-170, n), rep(1, n)),
    rep(-log(2 * pi), n)
  )
  # rows with and without tau = Inf in one call keep their own densities
  mixed <- .circmix_vp_ld(matrix(cos(y)), matrix(0, n), rep(-Inf, n), rep(5, n), c(Inf, 2, 0))
  expect_equal(mixed[1], -log(2 * pi))
  expect_equal(mixed[2], .circmix_vp_ld(matrix(cos(1)), matrix(0), -Inf, 5, 2))
  expect_equal(mixed[3], .circmix_ld(matrix(cos(pi)), matrix(0), -Inf, 5))
  # a scalar kappa recycles over the tau vector, with or without an Inf in it
  expect_equal(.circmix_vp_ld(matrix(cos(y)), matrix(0, n), rep(-Inf, n), 5, c(Inf, 2, 0)), mixed)
  expect_equal(
    .circmix_vp_ld(matrix(cos(y)), matrix(0, n), rep(-Inf, n), 5, c(2, 2, 0)),
    .circmix_vp_ld(matrix(cos(y)), matrix(0, n), rep(-Inf, n), rep(5, n), c(2, 2, 0))
  )
})

# Near the tau = Inf limit the memory density exceeds 1 / (2 pi) by
# pi sqrt(2) J(kappa) / sqrt(tau) at most, so below 1e-17 of that ratio it is
# uniform to double precision. Integrating there failed: past tau 7e305 the
# tail grid's 2 pi kappa overflowed (NA), and past 4.5e306 its end, 40 tau, did.
test_that(".circmix_vp_ld() is uniform, not NA, as tau nears the limit of a double", {
  y <- c(0, 1, pi)
  cells <- expand.grid(y = y, tau = c(1e306, 4e306, 5e306, 1e308))
  n <- nrow(cells)
  expect_identical(
    .circmix_vp_ld(matrix(cos(cells$y)), matrix(log(0.8), n), rep(log(0.2), n), 5, cells$tau),
    rep(-log(2 * pi), n)
  )
  # J(kappa) / sqrt(tau) = 1e-10 is integrated, and the excess is the bound;
  # just below 1e-17 the fallback takes over
  J <- .circmix_J(5)
  ratio <- c(1e-10, 0.99e-17)
  ld <- .circmix_vp_ld(matrix(1, 2), matrix(0, 2), rep(-Inf, 2), 5, (J / ratio)^2)
  expect_equal((ld[1] + log(2 * pi)) / ratio[1], pi * sqrt(2), tolerance = 1e-2)
  expect_identical(ld[2], -log(2 * pi))
  # kappa 1e140 is integrated: at tau 1e306 the grid reaches kappa 4e307, and
  # at tau 1e308 it would end past the largest double
  cells <- expand.grid(y = y, tau = c(1e306, 1e308))
  expect_equal(
    .circmix_vp_ld(matrix(cos(cells$y)), matrix(0, 6), rep(-Inf, 6), 1e140, cells$tau),
    rep(-log(2 * pi), 6),
    tolerance = 1e-10
  )
  expect_equal(
    .circmix_log_besselI0_scaled(c(1e300, 1e308)),
    -0.5 * (log(2 * pi) + log(c(1e300, 1e308))),
    tolerance = 1e-14
  )
})

test_that(".circmix_vp_ld() is NA where tau is NA, as it is where kappa is", {
  y <- c(0, 1, pi)
  ld <- .circmix_vp_ld(matrix(cos(y)), matrix(log(0.8), 3), rep(log(0.2), 3), 5, c(NA, 2, NaN))
  expect_true(all(is.na(ld[c(1, 3)])))
  expect_equal(ld[2], .circmix_vp_ld(matrix(cos(1)), matrix(log(0.8)), log(0.2), 5, 2))
  expect_true(is.na(.circmix_vp_ld(matrix(1), matrix(0), -Inf, NA, 2)))
})

test_that(".circmix_vp_ld() is a normalised density on the circle", {
  total <- stats::integrate(
    function(y) {
      n <- length(y)
      exp(.circmix_vp_ld(
        matrix(cos(y)), matrix(log(0.8), nrow = n), rep(log(0.2), n),
        rep(10, n), rep(2, n)
      ))
    },
    -pi, pi
  )
  expect_equal(total$value, 1, tolerance = 1e-6)
})

# The grid follows each component's integrand rather than the prior, so the
# density integrates to one only up to the quadrature error. The part of that
# error that depends on the shape alone acts as a log prior on tau: an offset of
# 2.5e-4 moved the maximum-likelihood estimate of log tau by 0.13 standard
# errors at N = 5000.
test_that(".circmix_vp_ld() is a normalised density at every shape", {
  total <- function(shape) {
    tau <- .circmix_J(10) / shape
    stats::integrate(
      function(y) {
        n <- length(y)
        exp(.circmix_vp_ld(
          matrix(cos(y)), matrix(log(0.8), nrow = n), rep(log(0.2), n),
          rep(10, n), rep(tau, n)
        ))
      },
      -pi, pi,
      rel.tol = 1e-12
    )$value
  }
  # below shape 2 the totals miss 1 by 2.5e-11 to 1.7e-9
  for (shape in c(0.05, 0.3, 0.75, 1)) {
    expect_lt(abs(total(shape) - 1), 1e-8, label = paste("shape", shape))
  }
})

# Integrates in log J from J = 1e-30 rather than from a gamma quantile: for a
# response far from the target the integrand peaks far below the prior's bulk,
# where a quantile bound cuts it off (and returns -Inf at large kappa). Below
# 1e-30 the von Mises is uniform to sqrt(2e-30) = 1.4e-15, so that mass enters
# as uniform; at kappa 10 it is 92% of the prior at shape 0.001 and 22% at
# shape 0.02. With 600 pieces instead of 300, or the cut at 1e-60, the value
# moves by at most 2e-13 at shapes 0.001 to 1.9.
vp_reference <- function(y, kappa, tau, p_mem, mu = 0, w = 1) {
  shape <- .circmix_J(kappa) / tau
  integrand <- function(t) {
    J <- exp(t)
    k <- .circmix_kappa(J)
    memory <- Reduce(`+`, Map(function(m, wm) wm * exp(k * cos(y - m) - .circmix_log_besselI0(k)), mu, w))
    density <- (p_mem * memory + (1 - p_mem)) / (2 * pi)
    density * stats::dgamma(J, shape = shape, scale = tau) * J
  }
  uniform <- stats::pgamma(1e-30, shape, scale = tau) * (p_mem * sum(w) + 1 - p_mem) / (2 * pi)
  upper <- log(stats::qgamma(1e-20, shape, scale = tau, lower.tail = FALSE))
  breaks <- seq(log(1e-30), upper, length.out = 300)
  log(uniform + sum(vapply(seq_len(length(breaks) - 1), function(i) {
    stats::integrate(integrand, breaks[i], breaks[i + 1],
      rel.tol = 1e-12, abs.tol = 0, subdivisions = 2000L
    )$value
  }, numeric(1))))
}

vp_errors <- function(y, kappa, shape, p_mem, ...) {
  tau <- .circmix_J(kappa) / shape
  n <- length(y)
  .circmix_vp_ld(
    matrix(cos(y)), matrix(log(p_mem), nrow = n), rep(log1p(-p_mem), n),
    rep(kappa, n), rep(tau, n), ...
  ) - vapply(y, vp_reference, numeric(1), kappa = kappa, tau = tau, p_mem = p_mem)
}

# Below shape 2 the errors at these cells are 2e-10 to 4e-9, so 1e-7 leaves a
# margin of 25 and still catches a slip of 0.1% in one term of the analytic
# tail. At kappa 0.001 and tau 1e-6 the gamma quantile that ends the tail grid
# lies below J0 = 1e-4, so the grid takes its fallback span; that cell is held
# to 1e-9 (error 5e-11).
test_that(".circmix_vp_ld() matches adaptive integration at every shape", {
  skip_on_cran()
  y <- c(0, 0.3, 1.5, pi)
  shapes <- c(0.001, 0.02, 0.1, 0.3, 0.74, 1, 1.9, 2, 2.1, 10)
  for (shape in shapes) {
    expect_lt(
      max(abs(vp_errors(y, kappa = 10, shape = shape, p_mem = 0.8))),
      1e-7,
      label = paste("shape", shape)
    )
  }
  expect_lt(
    max(abs(vp_errors(y, kappa = 0.001, shape = .circmix_J(0.001) / 1e-6, p_mem = 0.8))),
    1e-9
  )
})

# The tail grid's step depends on the node count, which 81 nodes alone cannot
# pin. Errors at kappa 10: 41 nodes 2.8e-5 (shape 0.02) and 5.7e-5 (shape 1.5),
# 121 nodes 5.5e-10 and 1e-12; each bound leaves a margin of at least 5.
test_that(".circmix_vp_ld() holds its accuracy below shape 2 at 41 and 121 nodes", {
  skip_on_cran()
  y <- c(0, 0.3, 1.5, pi)
  cases <- list(
    list(nodes = 41L, shape = 0.02, bound = 1.5e-4),
    list(nodes = 41L, shape = 1.5, bound = 3e-4),
    list(nodes = 121L, shape = 0.02, bound = 5e-9),
    list(nodes = 121L, shape = 1.5, bound = 5e-9)
  )
  for (case in cases) {
    expect_lt(
      max(abs(vp_errors(y, kappa = 10, shape = case$shape, p_mem = 0.8, nodes = case$nodes))),
      case$bound,
      label = paste("nodes", case$nodes, "shape", case$shape)
    )
  }
})

# Without a guess weight nothing masks the memory density of a response far
# from the target, whose integrand peaks many prior SDs below E[log J]: the
# capacity versions once an item holds a slot.
test_that(".circmix_vp_ld() is accurate for far responses without guessing", {
  skip_on_cran()
  y <- c(0, 1, 2, pi)
  expect_lt(max(abs(vp_errors(y, kappa = 200, shape = 10, p_mem = 1))), 1e-8)
  expect_lt(max(abs(vp_errors(y, kappa = 50, shape = 100, p_mem = 1))), 1e-8)
  # errors 1.1e-8, 1.5e-8 and 4.2e-9
  expect_lt(max(abs(vp_errors(y, kappa = 50, shape = 0.75, p_mem = 1))), 1e-7)
  expect_lt(max(abs(vp_errors(y, kappa = 50, shape = 0.3, p_mem = 1))), 1e-7)
  expect_lt(max(abs(vp_errors(y, kappa = 50, shape = 0.05, p_mem = 1))), 1e-7)

  # each component needs its own window when the components sit apart
  tau <- .circmix_J(100) / 20
  y <- c(0.3, 1.25, 2.4)
  n <- length(y)
  expect_lt(
    max(abs(
      .circmix_vp_ld(
        cbind(cos(y), cos(y - 2.5)), matrix(log(c(0.7, 0.3)), n, 2, byrow = TRUE),
        rep(-Inf, n), rep(100, n), rep(tau, n)
      ) - vapply(y, vp_reference, numeric(1),
        kappa = 100, tau = tau, p_mem = 1, mu = c(0, 2.5), w = c(0.7, 0.3)
      )
    )),
    1e-8
  )
})

# Summing kappa cos_d - log I0(kappa) cancels the gamma weight of a grid node
# away once kappa there nears 1e17 (kappa 1e12 already costs 2e-5), which moved the log
# density of a response on the item up by 2.4 nats at tau 1e20 and 5.5 at 1e200.
# At a huge tau nearly all precision mass sits at J ~ 0, so that response is
# uniform to within 1e-7 here. At a huge kappa, kappa(J) = J + 1/2 and the
# von Mises density at its mean is sqrt(kappa / (2 pi)) (1 + O(1 / kappa)), so
# the density is E[sqrt(J / (2 pi))] to within 1e-9 from kappa 1e12 on.
test_that(".circmix_vp_ld() is accurate at a huge kappa or tau for a response on the item", {
  tau <- c(1e20, 1e200, 1e20, 1e200)
  kappa <- c(5, 5, 50, 50)
  expect_equal(
    .circmix_vp_ld(matrix(1, 4), matrix(0, 4), rep(-Inf, 4), kappa, tau),
    rep(-log(2 * pi), 4),
    tolerance = 1e-7
  )

  # past kappa ~1e16, 1 / tau + 1 - cos_d rounds to 0 at cos_d = 1 unless the
  # mode scheme's start groups (1 - cos_d) first; past kappa ~1e154 the
  # slopes were 0/0
  cells <- expand.grid(kappa = c(1e12, 1e15, 1e17, 1e100, 1e200, 1e300), shape = c(2, 5, 100))
  tau <- .circmix_J(cells$kappa) / cells$shape
  expect_lt(
    max(abs(
      .circmix_vp_ld(matrix(1, 18), matrix(0, 18), rep(-Inf, 18), cells$kappa, tau) -
        (0.5 * log(tau / (2 * pi)) + lgamma(cells$shape + 0.5) - lgamma(cells$shape))
    )),
    1e-9
  )
})

# For a response on the item and tau = J, the first slope is
# shape - 1 + (1 - A) d kappa / d log J = shape - 1 + A / (1 + A) and the second
# tends to -1. A rounds to 1 from kappa ~1e16, which took the drift to 0 when it
# was computed as (cos_d - A), and past kappa ~1e154 dA / d kappa and (1 - A)^2
# both underflow, which made the second slope 0/0.
test_that(".circmix_vp_slopes() holds at a huge kappa for a response on the item", {
  J <- 10^c(12, 16, 50, 155, 300)
  slopes <- .circmix_vp_slopes(log(J), 1, 2, J, .circmix_kappa_table())
  expect_equal(slopes$d1, rep(1.5, 5), tolerance = 1e-12)
  expect_equal(slopes$d2, rep(-1, 5), tolerance = 1e-12)
})

# Below a gamma shape of 2 the J < 1e-4 part of the integral is analytic and
# the grid starts at 1e-4; above it the grid follows the integrand. With
# intercept-only kappa and tau every row crosses at once, so the step between
# the two schemes at shape 2 is a step in the whole log likelihood. Measured at
# y = 0, 1, pi: 1.2e-8 to 1.8e-8 at 81 nodes for kappa 1 to 100 (2.7e-7 at
# kappa 1000), and at 41 nodes 3.6e-5, 1.8e-4 and 5.8e-4 at kappa 1, 10 and
# 100. The bounds hold each step to about five times its size, so a change to
# either scheme shows up.
test_that("the two variable-precision schemes meet at shape 2", {
  tab <- .circmix_kappa_table()
  cos_d <- cos(c(0, 1, pi))
  cases <- list(
    list(nodes = 81L, kappa = 1, bound = 1e-7),
    list(nodes = 81L, kappa = 10, bound = 1e-7),
    list(nodes = 81L, kappa = 100, bound = 1e-7),
    list(nodes = 41L, kappa = 1, bound = 2e-4),
    list(nodes = 41L, kappa = 10, bound = 1e-3),
    list(nodes = 41L, kappa = 100, bound = 3e-3)
  )
  for (case in cases) {
    shape <- rep(2, 3)
    tau <- rep(.circmix_J(case$kappa) / 2, 3)
    jump <- .circmix_vp_component_tail(cos_d, shape, tau, case$nodes, tab) -
      .circmix_vp_component_mode(cos_d, shape, tau, case$nodes, tab)
    expect_lt(max(abs(jump)), case$bound, label = paste("nodes", case$nodes, "kappa", case$kappa))
  }
})

# From shape 2 up the likelihood stays the mode-centred one, so that existing
# fits keep their values; shapes 3 and 20 catch the switch moving up, 2 - 1e-9
# and 1.5 catch it moving down.
test_that(".circmix_vp_component() takes the tail scheme below shape 2 and the mode scheme from 2", {
  tab <- .circmix_kappa_table()
  shape <- c(0.3, 1.5, 2 - 1e-9, 2, 2 + 1e-9, 3, 20)
  cos_d <- cos(seq(0, pi, length.out = length(shape)))
  tau <- .circmix_J(10) / shape
  below <- shape < 2
  expect_identical(
    .circmix_vp_component(cos_d, shape, tau, 81L, tab),
    c(
      .circmix_vp_component_tail(cos_d[below], shape[below], tau[below], 81L, tab),
      .circmix_vp_component_mode(cos_d[!below], shape[!below], tau[!below], 81L, tab)
    )
  )
})

test_that(".circmix_slots() gives a slot distribution that is continuous in K", {
  slot_weights <- function(K, set_size, max_slots) {
    allocation <- .circmix_slots(K, set_size)
    weights <- numeric(max_slots + 1)
    weights[allocation$slots + 1] <- 1 - allocation$extra
    weights[allocation$slots + 2] <- weights[allocation$slots + 2] + allocation$extra
    weights
  }
  expect_equal(slot_weights(6 - 1e-9, 3, 4), slot_weights(6, 3, 4), tolerance = 1e-8)
  expect_equal(slot_weights(3 - 1e-9, 3, 4), slot_weights(3, 3, 4), tolerance = 1e-8)
  expect_equal(.circmix_slots(3.5, 1), list(slots = 3, extra = 0.5))
})

slot_averaging_density <- function(y, K, set_size, kappa = 6, tau = 0) {
  n <- length(y)
  .circmix_slot_averaging_ld(
    matrix(cos(y)), matrix(0, nrow = n), rep(K, n), rep(set_size, n),
    rep(kappa, n), rep(tau, n), 81L
  )
}

test_that(".circmix_slot_averaging_ld() is a normalised density on the circle", {
  for (case in list(c(3, 5, 0), c(7, 2, 0), c(7, 2, 1))) {
    total <- stats::integrate(
      function(y) exp(slot_averaging_density(y, case[1], case[2], tau = case[3])),
      -pi, pi
    )
    expect_equal(total$value, 1, tolerance = 1e-6)
  }
})

test_that(".circmix_slot_averaging_ld() is accurate for far responses under variable precision", {
  y <- c(1.5, 2, pi)
  three_slots <- .circmix_kappa(3 * .circmix_J(20))
  expect_lt(
    max(abs(
      slot_averaging_density(y, K = 3, set_size = 1, kappa = 20, tau = 2) -
        vapply(y, vp_reference, numeric(1), kappa = three_slots, tau = 2, p_mem = 1)
    )),
    1e-8
  )
  # one slot at gamma shape 0.2, so the three the item holds sit at 0.6, below
  # the switch at 2 and below 0.74, where the mode-centred grid alone loses
  # 1e-5; the error is 9.7e-9
  tau <- .circmix_J(5) / 0.2
  three_slots <- .circmix_kappa(3 * .circmix_J(5))
  expect_lt(
    max(abs(
      slot_averaging_density(y, K = 3, set_size = 1, kappa = 5, tau = tau) -
        vapply(y, vp_reference, numeric(1), kappa = three_slots, tau = tau, p_mem = 1)
    )),
    1e-7
  )
})

test_that(".circmix_slot_averaging_ld() sums the information of the slots an item holds", {
  y <- seq(-pi, pi, length.out = 30)
  two_slots <- .circmix_kappa(2 * .circmix_J(6))
  expect_equal(
    slot_averaging_density(y, K = 8, set_size = 4),
    brms::dvon_mises(y, mu = 0, kappa = two_slots, log = TRUE),
    tolerance = 1e-8
  )
  # with fewer slots than items, the items that received none are guessed
  expect_equal(
    slot_averaging_density(y, K = 2, set_size = 4),
    ld_mixture(y, kappa = 6, w_mem = 0.5, w_guess = 0.5),
    tolerance = 1e-8
  )
})

test_that("the slot-averaging gradient in K is the right-hand slope at a whole K / set_size", {
  skip_unless_stan_tests()
  skip_on_cran()
  tab <- .circmix_kappa_table()
  code <- paste0(
    "functions {\n",
    read_lines2(system.file("stan_chunks", "circmix_functions.stan", package = "bmm")),
    "\n}\n",
    "data { int NK; vector[NK] logk; vector[NK] dlogk; real logJ_min; real dlogJ; }
    parameters { real<lower=0> K; }
    model {
      target += circmix_slot_averaging_ld([cos(0.3)]', [0.0]', K, 2, 6, 0.5, 81,
                                          logk, dlogk, logJ_min, dlogJ);
    }"
  )
  model <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(code))
  utils::capture.output(diagnosis <- model$diagnose(
    data = list(
      NK = length(tab$logkappa), logk = tab$logkappa, dlogk = tab$dlogkappa,
      logJ_min = tab$logJ_min, dlogJ = tab$dlogJ
    ),
    init = list(list(K = 4)),
    error = 1e6
  ))
  # the density has a kink in K there, so diagnose()'s central difference
  # straddles two slopes; the gradient is held to the slope on the side the
  # allocation belongs to (K = 4 is two slots each, the start of [4, 6)), taken
  # from the R density and moved to the log K scale Stan samples on
  ld <- function(K) {
    .circmix_slot_averaging_ld(matrix(cos(0.3)), matrix(0), K, 2, 6, 0.5, 81L)
  }
  h <- 1e-7
  right_slope <- (ld(4 + h) - ld(4)) / h
  expect_equal(diagnosis$gradients()$model, 4 * right_slope + 1, tolerance = 1e-5)
})

test_that(".circmix_slot_averaging_ld() mixes the two slot-count branches", {
  y <- seq(-pi, pi, length.out = 9)
  n <- length(y)
  cosd <- cbind(cos(y), cos(y - 1.4))
  logw <- matrix(log(c(0.7, 0.3)), nrow = n, ncol = 2, byrow = TRUE)
  single <- .circmix_J(6)
  branch_ld <- function(slots) {
    .circmix_vp_ld(cosd, logw, rep(-Inf, n), rep(.circmix_kappa(slots * single), n),
                   rep(0.4, n), 81L)
  }
  # K = 3.5 over 2 items: 1 slot with probability 0.25, 2 with 0.75
  expect_equal(
    .circmix_slot_averaging_ld(cosd, logw, rep(3.5, n), rep(2, n), rep(6, n), rep(0.4, n), 81L),
    log(0.25 * exp(branch_ld(1)) + 0.75 * exp(branch_ld(2))),
    tolerance = 1e-12
  )
  # K = 3.5 over 4 items: no slot (guessed) with probability 0.125
  expect_equal(
    .circmix_slot_averaging_ld(cosd, logw, rep(3.5, n), rep(4, n), rep(6, n), rep(0.4, n), 81L),
    log(0.125 / (2 * pi) + 0.875 * exp(branch_ld(1))),
    tolerance = 1e-12
  )
  # a missing K leaves its own row NA and the others untouched
  expect_identical(
    is.na(.circmix_slot_averaging_ld(cosd[1:2, ], logw[1:2, ], c(NA, 3.5), c(2, 2), c(6, 6), c(0.4, 0.4), 81L)),
    c(TRUE, FALSE)
  )
})

# Compares the sampler against the density it is meant to draw from. n is large
# enough that a bin proportion has a standard error near 5e-4, so the tolerance
# below is about ten standard errors and does not need a seed to be stable.
expect_samples_match_density <- function(draws, density, n_bins = 20) {
  breaks <- seq(-pi, pi, length.out = n_bins + 1)
  observed <- as.numeric(table(cut(draws, breaks))) / length(draws)
  expected <- vapply(seq_len(n_bins), function(i) {
    stats::integrate(density, breaks[i], breaks[i + 1])$value
  }, numeric(1))
  expect_lt(max(abs(observed - expected)), 0.005)
}

test_that(".rcircmix() draws from the constant-precision mixture density", {
  n <- 2e5
  mu_nt <- 1.2
  draws <- .rcircmix(
    mu = matrix(c(0, mu_nt), nrow = n, ncol = 2, byrow = TRUE),
    logw = matrix(log(c(0.5, 0.3)), nrow = n, ncol = 2, byrow = TRUE),
    logw_guess = rep(log(0.2), n), kappa = rep(8, n), tau = rep(0, n)
  )
  expect_samples_match_density(draws, function(y) {
    exp(ld_mixture(y, mu_nt = mu_nt, kappa = 8, w_mem = c(0.5, 0.3), w_guess = 0.2))
  })
})

test_that(".rcircmix() draws from the variable-precision mixture density", {
  n <- 2e5
  draws <- .rcircmix(
    mu = matrix(0, nrow = n, ncol = 1),
    logw = matrix(log(0.8), nrow = n, ncol = 1),
    logw_guess = rep(log(0.2), n), kappa = rep(10, n), tau = rep(2, n)
  )
  expect_samples_match_density(draws, function(y) {
    exp(.circmix_vp_ld(
      matrix(cos(y)), matrix(log(0.8), nrow = length(y)),
      rep(log(0.2), length(y)), rep(10, length(y)), rep(2, length(y))
    ))
  })
})

test_that(".rcircmix_slot_averaging() draws from the slot-averaging density", {
  n <- 2e5
  draws <- .rcircmix_slot_averaging(
    mu = matrix(0, nrow = n), logw = matrix(0, nrow = n),
    K = rep(7, n), set_size = rep(2, n), kappa = rep(6, n), tau = rep(1, n)
  )
  expect_samples_match_density(draws, function(y) {
    exp(slot_averaging_density(y, K = 7, set_size = 2, tau = 1))
  })
})

test_that(".circmix_add_variable_precision() inserts tau next to kappa", {
  spec <- list(
    parameters = list(mu = "location", kappa = "precision", thetat = "weight"),
    links = list(mu = "tan_half", kappa = "log", thetat = "logit"),
    priors = list(kappa = list(main = "normal(2, 1)")),
    init_ranges = list(kappa = c(3, 8))
  )
  out <- .circmix_add_variable_precision(spec)
  expect_equal(names(out$parameters), c("mu", "kappa", "tau", "thetat"))
  expect_equal(names(out$links), c("mu", "kappa", "tau", "thetat"))
  expect_equal(out$links$tau, "log")
  expect_named(out$priors$tau, c("main", "effects", "sd"))
  expect_false(is.null(out$init_ranges$tau))
})

test_that(".circmix_recycle() extends every argument to the longest", {
  out <- .circmix_recycle(x = 1:3, mu = 0, kappa = c(5, 6, 7))
  expect_equal(out$x, 1:3)
  expect_equal(out$mu, c(0, 0, 0))
  expect_equal(out$kappa, c(5, 6, 7))
})

test_that(".circmix_cos() reads only as many non-targets as the set size has", {
  nt <- c(1, 2, 3)
  expect_equal(dim(.circmix_cos(0.5, 0, nt, set_size = 1)), c(1L, 1L))
  expect_equal(dim(.circmix_cos(0.5, 0, nt, set_size = 3)), c(1L, 3L))
  expect_equal(
    .circmix_cos(0.5, 0, nt, set_size = 3)[1, ],
    cos(0.5 - c(0, 1, 2))
  )
  # padding beyond the set size cannot change the result
  expect_equal(
    .circmix_cos(0.5, 0, c(1, 2, 999), set_size = 3),
    .circmix_cos(0.5, 0, c(1, 2, -7), set_size = 3)
  )
})

test_that(".circmix_bounds() derives natural-scale bounds from the links", {
  bounds <- .circmix_bounds(list(mu = "tan_half", kappa = "log", thetat = "logit"))
  expect_equal(bounds$lb, c(NA, 0, 0))
  expect_equal(bounds$ub, c(NA, NA, 1))

  bounds <- .circmix_bounds(list(kappa = "softplus", thetat = "probit", x = "identity"))
  expect_equal(bounds$lb, c(0, 0, NA))
  expect_equal(bounds$ub, c(NA, 1, NA))
})

test_that(".circmix_prep_tau() is zero when the model has no tau parameter", {
  expect_equal(.circmix_prep_tau(list(dpars = list(kappa = 1), ndraws = 4), 1), rep(0, 4))
})

test_that(".circmix_prep_nodes() falls back to the default node count", {
  expect_equal(.circmix_prep_nodes(list(family = list())), 81L)
  expect_equal(.circmix_prep_nodes(list(family = list(vp_nodes = 41L))), 41L)
})

test_that(".circmix_check_variable_precision() rejects unusable node counts", {
  expect_error(.circmix_check_variable_precision(NA, 41L), "TRUE or FALSE")
  expect_error(.circmix_check_variable_precision(FALSE, 21L), "at least 41")
  expect_error(.circmix_check_variable_precision(FALSE, 42L), "at least 41")
  # the model stores as.integer(vp_nodes), so 42.5 would become an even 42
  expect_error(.circmix_check_variable_precision(FALSE, 42.5), "at least 41")
  expect_error(.circmix_check_variable_precision(FALSE, 41.5), "at least 41")
  expect_error(.circmix_check_variable_precision(FALSE, NA), "at least 41")
  expect_error(.circmix_check_variable_precision(FALSE, c(41, 43)), "at least 41")
  expect_silent(.circmix_check_variable_precision(TRUE, 41L))
  expect_silent(.circmix_check_variable_precision(TRUE, 81L))
  expect_silent(.circmix_check_variable_precision(TRUE, 43))
})

test_that(".circmix_prep_nt() returns the covariates in numeric order", {
  prep <- list(data = list(
    Y = 1:3, vreal10 = c(10, 0, 0), vreal2 = c(2, 0, 0), vreal1 = c(1, 0, 0)
  ))
  expect_equal(unname(.circmix_prep_nt(prep, 1)), c(1, 2, 10))
})

test_that(".circmix_custom_family() adds tau only when the model estimates it", {
  model <- list(
    links = list(mu = "tan_half", kappa = "log", tau = "log", thetat = "logit"),
    vp_nodes = 41L, variable_precision = FALSE
  )
  stub_lik <- function(i, prep) NULL
  stub_pred <- function(i, prep, ...) NULL
  stub_epred <- posterior_epred_undefined("demo")
  constant <- .circmix_custom_family(
    model, "demo", "thetat",
    log_lik = stub_lik, posterior_predict = stub_pred,
    posterior_epred = stub_epred
  )
  expect_identical(constant$posterior_epred, stub_epred)
  expect_equal(constant$dpars, c("mu", "kappa", "thetat"))
  expect_equal(constant$core_dpars, c("mu", "kappa", "0.0", "thetat"))

  model$variable_precision <- TRUE
  varying <- .circmix_custom_family(
    model, "demo", "thetat",
    log_lik = stub_lik, posterior_predict = stub_pred,
    posterior_epred = stub_epred
  )
  expect_equal(varying$dpars, c("mu", "kappa", "tau", "thetat"))
  expect_equal(varying$core_dpars, varying$dpars)
  # brms normalises the bounds into a list of Stan bound expressions
  expect_equal(as.numeric(unlist(varying$lb)), c(NA, 0, 0, 0))
  expect_equal(as.numeric(unlist(varying$ub)), c(NA, NA, NA, 1))
})

test_that(".circmix_family_vars() emits only bare indices", {
  vars <- .circmix_family_vars(vint = TRUE, n_vreal = 3)
  indexed <- grep("\\[", vars, value = TRUE)
  expect_equal(indexed, c("vint1[n]", "vreal1[n]", "vreal2[n]", "vreal3[n]"))
  expect_true(all(grepl("^(vint|vreal)[0-9]+\\[n\\]$", indexed)))
  expect_true(all(.circmix_table_vars() %in% vars))
})

test_that(".circmix_aterm() builds the response formula for each covariate combination", {
  expect_equal(as.character(.circmix_aterm("dev_rad")), "dev_rad ~ 1")
  expect_equal(
    as.character(.circmix_aterm("dev_rad", vint = "ss_numeric", vreal = c("nt1", "nt2"))),
    "dev_rad | vint(ss_numeric) + vreal(nt1, nt2) ~ 1"
  )
})

test_that(".circmix_stan_wrapper() packs the vreal arguments and passes the table through", {
  wrapper <- .circmix_stan_wrapper(
    "imm_full",
    dpars = c("mu", "kappa", "tau", "c", "a", "s"),
    vint = "ss", vreal = list(nt = 2, dist = 2), nodes = 41L
  )
  expect_match(wrapper, "real imm_full_lpdf\\(real y, real mu, real kappa, real tau")
  expect_match(wrapper, "int ss, real nt1, real nt2, real dist1, real dist2")
  expect_match(wrapper, "to_vector\\(\\{nt1, nt2\\}\\), to_vector\\(\\{dist1, dist2\\}\\)", fixed = FALSE)
  expect_match(wrapper, "imm_full_core\\(")
  expect_match(wrapper, "41, circmix_logk, circmix_dlogk, circmix_logJ_min, circmix_dlogJ")
})

test_that(".circmix_stan_wrapper() can pass the core a value its family omits", {
  wrapper <- .circmix_stan_wrapper(
    "mixture2p_simple",
    dpars = c("mu", "kappa", "thetat"),
    core_dpars = c("mu", "kappa", "0.0", "thetat")
  )
  # the whole signature, because an empty argument left by a zero-length
  # covariate group still matches any prefix assertion
  expect_match(
    wrapper,
    paste0(
      "real mixture2p_simple_lpdf\\(real y, real mu, real kappa, real thetat, ",
      "data vector circmix_logk, data vector circmix_dlogk, ",
      "data real circmix_logJ_min, data real circmix_dlogJ\\) \\{"
    )
  )
  expect_match(wrapper, "mixture2p_simple_core\\(y, mu, kappa, 0.0, thetat, 81")
})

test_that(".circmix_stan_wrapper() emits Stan that parses", {
  skip_unless_stan_tests()
  skip_on_cran()
  wrappers <- c(
    .circmix_stan_wrapper("no_covariates", c("mu", "kappa")),
    .circmix_stan_wrapper("with_set_size", c("mu", "kappa"), vint = "ss"),
    .circmix_stan_wrapper(
      "with_covariates", c("mu", "kappa", "tau"),
      vint = "ss", vreal = list(nt = 2, dist = 2)
    )
  )
  stubs <- c(
    "real no_covariates_core(real y, real mu, real kappa, int n, data vector a, data vector b, data real c, data real d) { return 0; }",
    "real with_set_size_core(real y, real mu, real kappa, int ss, int n, data vector a, data vector b, data real c, data real d) { return 0; }",
    "real with_covariates_core(real y, real mu, real kappa, real tau, int ss, vector nt, vector dist, int n, data vector a, data vector b, data real c, data real d) { return 0; }"
  )
  code <- paste0(
    "functions {\n", paste(stubs, collapse = "\n"), "\n",
    paste(wrappers, collapse = "\n"), "\n}\nmodel {}\n"
  )
  model <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(code), compile = FALSE)
  expect_silent(model$check_syntax(quiet = TRUE))
})

# The functions are evaluated through a fixed_param run rather than through
# cmdstanr::expose_functions(), which links against RcppParallel's libtbb and so
# fails on toolchains where that disagrees with the one cmdstan was built with.
# sig_figs is raised because cmdstan writes only 6 significant digits by default,
# which alone accounts for differences of ~1e-7. The cells *_e are one-component
# variable-precision rows with their own kappa and tau; a call without them
# passes none.
circmix_stan_values <- function(data) {
  data <- utils::modifyList(list(
    NE = 0L, y_e = numeric(0), logw_e = numeric(0), logw_guess_e = numeric(0),
    kappa_e = numeric(0), tau_e = numeric(0)
  ), data)
  chunks <- c("circmix_functions.stan", "mixture2p_functions.stan", "mixture3p_functions.stan", "imm_functions.stan")
  code <- paste0(
    "functions {\n",
    paste(vapply(chunks, function(chunk) {
      read_lines2(system.file("stan_chunks", chunk, package = "bmm"))
    }, character(1)), collapse = "\n"),
    "\n}\n",
    "data {
      int NK; vector[NK] logk; vector[NK] dlogk; real logJ_min; real dlogJ;
      int NC; vector[NC] cosd; vector[NC] logw; real logw_guess; vector[NC] kappa_het;
      int NKAP; vector[NKAP] kappa;
      int NJ; vector[NJ] Jtest;
      int NTAU; vector[NTAU] tau; real kappa_vp; int nodes;
      int NI0; vector[NI0] kappa_i0;
      real K_slots; int ss_slots; int ss_empty;
      vector[2] cosd_far; vector[2] logw_far; real kappa_far; int NFAR; vector[NFAR] tau_far;
      vector[3] nt_loc; int ss_loc; vector[2] dist_w;
      int NE; vector[NE] y_e; vector[NE] logw_e; vector[NE] logw_guess_e;
      vector[NE] kappa_e; vector[NE] tau_e;
    }
    generated quantities {
      vector[NKAP] out_J; vector[NKAP] out_ld; vector[NJ] out_kappa;
      real out_het; vector[NTAU] out_vp; vector[2] out_slots;
      vector[NTAU] out_sa; vector[NFAR] out_far; vector[NI0] out_logi0; real out_vp0;
      vector[4] out_branches; vector[4] out_branches_empty; vector[ss_loc] out_loc;
      vector[NTAU] out_shape;
      vector[2] out_w2_simple; vector[2] out_w2_slot; vector[4] out_w3_simple;
      vector[4] out_w3_slot; vector[3] out_w3_sa; vector[4] out_imm_abc;
      vector[4] out_imm_bsc; vector[4] out_imm_full; vector[2] out_w2_full;
      vector[2] out_w3_slot1; vector[1] out_w3_sa1;
      vector[NE] out_e_shape; vector[NE] out_e_vp; vector[NE] out_e_sa;
      for (i in 1:NKAP) {
        out_J[i] = circmix_J(kappa[i]);
        out_ld[i] = circmix_ld(cosd, logw, logw_guess, kappa[i]);
      }
      for (i in 1:NJ) {
        out_kappa[i] = circmix_kappa(Jtest[i], logk, dlogk, logJ_min, dlogJ);
      }
      out_het = circmix_het_ld(cosd, logw, logw_guess, kappa_het);
      for (i in 1:NI0) {
        out_logi0[i] = log_modified_bessel_first_kind(0, kappa_i0[i]);
      }
      for (i in 1:NTAU) {
        out_vp[i] = circmix_vp_ld(cosd, logw, logw_guess, kappa_vp, tau[i], nodes,
                                  logk, dlogk, logJ_min, dlogJ);
        out_sa[i] = circmix_slot_averaging_ld(cosd, logw, K_slots, ss_slots,
                                              kappa_vp, tau[i], nodes,
                                              logk, dlogk, logJ_min, dlogJ);
      }
      out_slots = circmix_slots(K_slots, ss_slots);
      out_vp0 = circmix_vp_ld(cosd, logw, logw_guess, 0, 0.5, nodes,
                              logk, dlogk, logJ_min, dlogJ);
      for (i in 1:NFAR) {
        out_far[i] = circmix_vp_ld(cosd_far, logw_far, negative_infinity(), kappa_far,
                                   tau_far[i], nodes, logk, dlogk, logJ_min, dlogJ);
      }
      for (i in 1:NTAU) {
        out_shape[i] = circmix_vp_shape(kappa_vp, tau[i]);
      }
      out_branches = circmix_slot_averaging_branches(K_slots, ss_slots, kappa_vp,
                                                     logk, dlogk, logJ_min, dlogJ);
      out_branches_empty = circmix_slot_averaging_branches(K_slots, ss_empty,
                                                           kappa_vp, logk, dlogk,
                                                           logJ_min, dlogJ);
      out_loc = circmix_locations(0.3, nt_loc, ss_loc);
      out_w2_simple = mixture2p_simple_logw(0.7);
      out_w2_slot = mixture2p_slot_logw(2.5, 3);
      out_w3_simple = mixture3p_simple_logw(1.2, 0.3, 3);
      out_w3_slot = mixture3p_slot_logw(2.5, 0.8, 3);
      out_w3_sa = mixture3p_slot_averaging_logw(0.8, 3);
      out_imm_abc = imm_abc_logw(2, 0.5, 1, 3);
      out_imm_bsc = imm_bsc_logw(2, 3, 1, 3, dist_w);
      out_imm_full = imm_full_logw(2, 0.5, 3, 1, 3, dist_w);
      out_w2_full = mixture2p_slot_logw(3.5, 3);
      out_w3_slot1 = mixture3p_slot_logw(0.6, 0.8, 1);
      out_w3_sa1 = mixture3p_slot_averaging_logw(0.8, 1);
      for (i in 1:NE) {
        out_e_shape[i] = circmix_vp_shape(kappa_e[i], tau_e[i]);
        out_e_vp[i] = circmix_vp_ld([cos(y_e[i])]', [logw_e[i]]', logw_guess_e[i], kappa_e[i],
                                    tau_e[i], nodes, logk, dlogk, logJ_min, dlogJ);
        out_e_sa[i] = circmix_slot_averaging_ld([cos(y_e[i])]', [0.0]', K_slots, ss_slots,
                                                kappa_e[i], tau_e[i], nodes,
                                                logk, dlogk, logJ_min, dlogJ);
      }
    }"
  )
  model <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(code))
  fit <- model$sample(
    data = data, fixed_param = TRUE, iter_sampling = 1, chains = 1,
    refresh = 0, show_messages = FALSE, sig_figs = 18
  )
  function(variable) as.numeric(fit$draws(variable))
}

test_that("the Stan and R implementations of the shared core agree", {
  skip_unless_stan_tests()
  skip_on_cran()

  tab <- .circmix_kappa_table()
  y <- 0.4
  cosd <- cos(y - c(0, 1.2, -2.0))
  logw <- log(c(0.4, 0.2, 0.15))
  logw_guess <- log(0.25)
  kappa <- c(0.05, 0.5, 2, 10, 100, 1000, 6e4)
  # at a node Stan's binary search and R's floor() pick different cells, which
  # must give the same value; only the point inside the last cell would see an
  # off-by-one that loses that cell
  n_tab <- length(tab$logkappa)
  node <- function(i) exp(tab$logJ_min + (i - 1) * tab$dlogJ)
  Jtest <- c(
    1e-7, 1e-3, 0.5, 9.486, 500, 5e4, 1e7,
    node(c(2, 1000, n_tab - 1, n_tab)), sqrt(node(n_tab - 1) * node(n_tab))
  )
  # either side of the point where besselI() underflows and R switches to the
  # asymptotic series; Stan's log_modified_bessel_first_kind() is independent
  kappa_i0 <- c(1, 1e3, 9e4, 9.99e4, 1.0005e5, 1.001e5, 1e6)
  kappa_het <- c(3, 12, 40)
  tau <- c(0, 0.5, 2, .circmix_J(10) / 0.75)
  # far from both components and no guessing: each component's own window
  cosd_far <- cos(2.9 - c(0, 1.2))
  logw_far <- log(c(0.7, 0.3))
  tau_far <- .circmix_J(200) / c(0.75, 10, 1000)

  stan <- circmix_stan_values(list(
    NK = length(tab$logkappa), logk = tab$logkappa, dlogk = tab$dlogkappa,
    logJ_min = tab$logJ_min, dlogJ = tab$dlogJ,
    NC = length(cosd), cosd = cosd, logw = logw, logw_guess = logw_guess,
    kappa_het = kappa_het, NKAP = length(kappa), kappa = kappa,
    NJ = length(Jtest), Jtest = Jtest, NTAU = length(tau), tau = tau,
    NI0 = length(kappa_i0), kappa_i0 = kappa_i0,
    kappa_vp = 10, nodes = 81L, K_slots = 3.5, ss_slots = 2L, ss_empty = 4L,
    cosd_far = cosd_far, logw_far = logw_far, kappa_far = 200,
    NFAR = length(tau_far), tau_far = tau_far,
    nt_loc = c(1, 2, 3), ss_loc = 3L, dist_w = c(0.5, 2)
  ))

  n <- length(kappa)
  as_rows <- function(x) matrix(x, nrow = n, ncol = length(x), byrow = TRUE)

  expect_equal(stan("out_J"), .circmix_J(kappa), tolerance = 1e-10)
  expect_lt(max(abs(stan("out_kappa") / .circmix_kappa(Jtest) - 1)), 1e-13)
  expect_lt(max(abs(stan("out_logi0") / .circmix_log_besselI0(kappa_i0) - 1)), 1e-14)
  expect_equal(
    stan("out_ld"),
    .circmix_ld(as_rows(cosd), as_rows(logw), rep(logw_guess, n), kappa),
    tolerance = 1e-12
  )
  expect_equal(
    stan("out_het"),
    .circmix_het_ld(
      matrix(cosd, nrow = 1), matrix(logw, nrow = 1), logw_guess,
      matrix(kappa_het, nrow = 1)
    ),
    tolerance = 1e-12
  )
  expect_equal(
    stan("out_vp"),
    vapply(tau, function(c_tau) {
      .circmix_vp_ld(matrix(cosd, nrow = 1), matrix(logw, nrow = 1), logw_guess, 10, c_tau)
    }, numeric(1)),
    tolerance = 1e-12
  )
  expect_equal(stan("out_slots"), unname(unlist(.circmix_slots(3.5, 2))))
  # J(0) = 0, so the gamma on J is a point mass at zero and memory is uniform
  expect_equal(
    stan("out_vp0"),
    .circmix_vp_ld(matrix(cosd, nrow = 1), matrix(logw, nrow = 1), logw_guess, 0, 0.5),
    tolerance = 1e-12
  )
  expect_equal(
    stan("out_far"),
    vapply(tau_far, function(c_tau) {
      .circmix_vp_ld(matrix(cosd_far, nrow = 1), matrix(logw_far, nrow = 1), -Inf, 200, c_tau)
    }, numeric(1)),
    tolerance = 1e-12
  )
  expect_equal(
    stan("out_sa"),
    vapply(tau, function(c_tau) {
      .circmix_slot_averaging_ld(
        matrix(cosd, nrow = 1), matrix(logw, nrow = 1), 3.5, 2L, 10, c_tau, 81L
      )
    }, numeric(1)),
    tolerance = 1e-12
  )

  shape <- .circmix_J(10) / tau
  expect_equal(stan("out_shape"), c(Inf, shape[-1]), tolerance = 1e-12)
  for (case in list(list("out_branches", 2L), list("out_branches_empty", 4L))) {
    branch <- .circmix_slot_averaging_branches(3.5, case[[2]], 10)
    expect_equal(stan(case[[1]]), as.numeric(unlist(branch)), tolerance = 1e-12)
  }
  expect_equal(stan("out_loc"), c(0.3, 1, 2))

  # each model's weight builder against the R twin its d*() function uses
  weights <- list(
    out_w2_simple = .mixture2p_weights(0.7),
    out_w2_slot = .mixture2p_weights(2.5 / 3),
    out_w3_simple = .mixture3p_softmax_weights(1.2, 0.3, 3),
    out_w3_slot = .mixture3p_nested_weights(2.5 / 3, 0.8, 3),
    out_imm_abc = .imm_log_weights(2, 0.5, 0, 1, 3, c(0.5, 2), "abc"),
    out_imm_bsc = .imm_log_weights(2, 1, 3, 1, 3, c(0.5, 2), "bsc"),
    out_imm_full = .imm_log_weights(2, 0.5, 3, 1, 3, c(0.5, 2), "full"),
    # K at or above the set size stores everything; at set size 1 there is no swap
    out_w2_full = .mixture2p_weights(1),
    out_w3_slot1 = .mixture3p_nested_weights(0.6, 0.8, 1)
  )
  for (name in names(weights)) {
    expect_equal(
      stan(name), c(weights[[name]]$logw, weights[[name]]$logw_guess),
      tolerance = 1e-12, label = name
    )
  }
  expect_equal(stan("out_w3_sa"), c(.mixture3p_swap_weights(0.8, 3)), tolerance = 1e-12)
  expect_equal(stan("out_w3_sa1"), c(.mixture3p_swap_weights(0.8, 1)))
})

test_that("variable precision says what kappa means, in every model and version", {
  models <- list(
    mixture2p("y", variable_precision = TRUE),
    mixture3p("y", "nt", "ss", version = "slot", variable_precision = TRUE),
    imm("y", "nt", "d", "ss", variable_precision = TRUE)
  )
  for (model in models) {
    expect_match(model$parameters$kappa, "Fisher information J(kappa) is the mean of J", fixed = TRUE)
    expect_no_match(model$parameters$kappa, "..", fixed = TRUE)
  }
  # under slot averaging an item held in j slots has mean J j * J(kappa), so
  # J(kappa) is the mean for one slot, not for a trial
  slot_averaging <- list(
    mixture2p("y", set_size = "ss", version = "slot_averaging", variable_precision = TRUE),
    mixture3p("y", "nt", "ss", version = "slot_averaging", variable_precision = TRUE)
  )
  for (model in slot_averaging) {
    expect_match(model$parameters$kappa, "mean Fisher information of one slot", fixed = TRUE)
    expect_no_match(model$parameters$kappa, "mean of J across trials", fixed = TRUE)
    expect_match(model$parameters$tau, "sum of j such draws", fixed = TRUE)
    expect_no_match(model$parameters$kappa, "..", fixed = TRUE)
  }
  expect_no_match(mixture2p("y")$parameters$kappa, "variable precision")
})

test_that("dmixture2p() is accurate at a small gamma shape at the default nodes", {
  # gamma shape J(1) / 20 = 0.0223, where most of the prior mass of J lies
  # below the grid; the error is 3.9e-10
  expect_lt(
    abs(dmixture2p(0.3, kappa = 1, tau = 20, log = TRUE) -
      vp_reference(0.3, kappa = 1, tau = 20, p_mem = 0.6)),
    1e-8
  )
})

test_that("Stan matches R at small gamma shapes and either side of the tail switch", {
  skip_unless_stan_tests()
  skip_on_cran()

  # 0.73 sits below 0.74, where the mode-centred grid alone loses 1e-5; 1.999
  # and 2.001 sit either side of the switch to the analytic tail, 1.9 and 2.1
  # further out; at 0.001 the floor of 40 tau ends the tail grid
  tab <- .circmix_kappa_table()
  cosd <- cos(0.4 - c(0, 1.2))
  logw <- log(c(0.5, 0.3))
  shapes <- c(0.001, 0.05, 0.3, 0.73, 1.9, 1.999, 2.001, 2.1)
  tau <- .circmix_J(10) / shapes
  tau_far <- .circmix_J(200) / c(0.05, 0.3, 1.9)
  cosd_far <- cos(2.9 - c(0, 1.2))
  logw_far <- log(c(0.7, 0.3))
  data <- list(
    NK = length(tab$logkappa), logk = tab$logkappa, dlogk = tab$dlogkappa,
    logJ_min = tab$logJ_min, dlogJ = tab$dlogJ,
    NC = 2L, cosd = cosd, logw = logw, logw_guess = log(0.2),
    kappa_het = c(3, 12), NKAP = 1L, kappa = array(2), NJ = 1L, Jtest = array(0.5),
    NTAU = length(tau), tau = tau, NI0 = 1L, kappa_i0 = array(1),
    kappa_vp = 10, nodes = 81L, K_slots = 1.5, ss_slots = 2L, ss_empty = 4L,
    cosd_far = cosd_far, logw_far = logw_far, kappa_far = 200,
    NFAR = length(tau_far), tau_far = tau_far,
    nt_loc = c(1, 2, 3), ss_loc = 3L, dist_w = c(0.5, 2)
  )
  stan <- circmix_stan_values(data)

  expect_equal(stan("out_shape"), shapes, tolerance = 1e-12)
  expect_equal(
    stan("out_vp"),
    vapply(tau, function(c_tau) {
      .circmix_vp_ld(matrix(cosd, nrow = 1), matrix(logw, nrow = 1), log(0.2), 10, c_tau)
    }, numeric(1)),
    tolerance = 1e-10
  )
  expect_equal(
    stan("out_far"),
    vapply(tau_far, function(c_tau) {
      .circmix_vp_ld(matrix(cosd_far, nrow = 1), matrix(logw_far, nrow = 1), -Inf, 200, c_tau)
    }, numeric(1)),
    tolerance = 1e-10
  )
  # two edges of the tail grid's span. With tau below 2.5e-6 the gamma quantile
  # that ends the grid lies below J0, so the grid takes its fallback span. At
  # shape 0.001 the floor of 40 tau ends it; a response at the target sees the
  # floor (a floor of 10 tau moves its log density by 2e-6), one at 0.4 barely
  # does (1e-10)
  tau_floor <- .circmix_J(10) / 0.001
  stan_edges <- circmix_stan_values(utils::modifyList(data, list(
    kappa_vp = 0.001, NTAU = 1L, tau = array(1e-6),
    cosd_far = cos(c(0, -1.2)), kappa_far = 10, NFAR = 1L, tau_far = array(tau_floor)
  )))
  expect_equal(
    stan_edges("out_vp"),
    .circmix_vp_ld(matrix(cosd, nrow = 1), matrix(logw, nrow = 1), log(0.2), 0.001, 1e-6),
    tolerance = 1e-10
  )
  expect_equal(
    stan_edges("out_far"),
    .circmix_vp_ld(matrix(cos(c(0, -1.2)), nrow = 1), matrix(logw_far, nrow = 1), -Inf, 10, tau_floor),
    tolerance = 1e-10
  )
  expect_equal(
    stan("out_sa"),
    vapply(tau, function(c_tau) {
      .circmix_slot_averaging_ld(
        matrix(cosd, nrow = 1), matrix(logw, nrow = 1), 1.5, 2L, 10, c_tau, 81L
      )
    }, numeric(1)),
    tolerance = 1e-10
  )
})

# Memory is uniform at tau = Inf, where J that underflows to 0, and wherever
# J(kappa) / sqrt(tau) < 1e-17 (.circmix_vp_uniform()). Integrating there gave
# NaN at a subnormal gamma shape (kappa 1e-155 at tau 1, 0.01 at 1e305) and NaN
# or +Inf once the tail grid overflowed (tau 1e306 on). kappa 1e140 is still
# integrated: at tau 1e306 on a grid that reaches kappa 4e307, and at 1e308 on
# one capped below the largest double.
test_that("Stan matches R where memory is uniform and near that limit", {
  skip_unless_stan_tests()
  skip_on_cran()

  tab <- .circmix_kappa_table()
  limits <- data.frame(
    kappa = c(5, 1e-170, 1e-155, 0.01, 5, 5, 5),
    tau = c(Inf, 1, 1, 1e305, 1e306, 5e306, 1e308)
  )
  cells <- merge(expand.grid(y = c(0, 1, pi), p_mem = c(1, 0.8)), limits)
  cells <- rbind(cells, expand.grid(y = c(0, pi), p_mem = 1, kappa = 1e140, tau = c(1e306, 1e308)))
  n <- nrow(cells)
  uniform <- seq_len(n - 4)
  stan <- circmix_stan_values(list(
    NK = length(tab$logkappa), logk = tab$logkappa, dlogk = tab$dlogkappa,
    logJ_min = tab$logJ_min, dlogJ = tab$dlogJ,
    NC = 1L, cosd = array(1), logw = array(0), logw_guess = -Inf,
    kappa_het = array(3), NKAP = 1L, kappa = array(2), NJ = 1L, Jtest = array(0.5),
    NTAU = 1L, tau = array(1), NI0 = 1L, kappa_i0 = array(1),
    kappa_vp = 5, nodes = 81L, K_slots = 1.5, ss_slots = 2L, ss_empty = 4L,
    cosd_far = c(1, 1), logw_far = log(c(0.7, 0.3)), kappa_far = 5,
    NFAR = 1L, tau_far = array(1),
    nt_loc = c(1, 2, 3), ss_loc = 3L, dist_w = c(0.5, 2),
    NE = n, y_e = cells$y, logw_e = log(cells$p_mem), logw_guess_e = log1p(-cells$p_mem),
    kappa_e = cells$kappa, tau_e = cells$tau
  ))

  expect_equal(stan("out_e_shape")[uniform], rep(Inf, length(uniform)))
  expect_equal(stan("out_e_vp"), rep(-log(2 * pi), n), tolerance = 1e-10)
  expect_equal(stan("out_e_vp")[uniform], rep(-log(2 * pi), length(uniform)), tolerance = 1e-15)
  expect_equal(
    stan("out_e_vp"),
    .circmix_vp_ld(
      matrix(cos(cells$y)), matrix(log(cells$p_mem)), log1p(-cells$p_mem),
      cells$kappa, cells$tau
    ),
    tolerance = 1e-12
  )
  expect_equal(
    stan("out_e_sa"),
    .circmix_slot_averaging_ld(
      matrix(cos(cells$y)), matrix(0, n), rep(1.5, n), rep(2L, n),
      cells$kappa, cells$tau, 81L
    ),
    tolerance = 1e-12
  )
})

# the R twin of this check is ".circmix_vp_ld() is accurate at a huge kappa or
# tau for a response on the item"; at tau 1e17 the tail grid reaches kappa 1e17
test_that("Stan matches R at a huge kappa or tau for a response on the item", {
  skip_unless_stan_tests()
  skip_on_cran()

  tab <- .circmix_kappa_table()
  cosd <- cos(c(0, 1.2))
  logw <- log(c(0.5, 0.3))
  tau <- c(1e17, 1e20, 1e200)
  data <- list(
    NK = length(tab$logkappa), logk = tab$logkappa, dlogk = tab$dlogkappa,
    logJ_min = tab$logJ_min, dlogJ = tab$dlogJ,
    NC = 2L, cosd = cosd, logw = logw, logw_guess = log(0.2),
    kappa_het = c(3, 12), NKAP = 1L, kappa = array(2), NJ = 1L, Jtest = array(0.5),
    NTAU = length(tau), tau = tau, NI0 = 1L, kappa_i0 = array(1),
    kappa_vp = 50, nodes = 81L, K_slots = 1.5, ss_slots = 2L, ss_empty = 4L,
    cosd_far = c(1, 1), logw_far = log(c(0.7, 0.3)), kappa_far = 5,
    NFAR = length(tau), tau_far = tau,
    nt_loc = c(1, 2, 3), ss_loc = 3L, dist_w = c(0.5, 2)
  )
  stan <- circmix_stan_values(data)

  expect_equal(
    stan("out_vp"),
    vapply(tau, function(c_tau) {
      .circmix_vp_ld(matrix(cosd, nrow = 1), matrix(logw, nrow = 1), log(0.2), 50, c_tau)
    }, numeric(1)),
    tolerance = 1e-10
  )
  expect_equal(
    stan("out_far"),
    .circmix_vp_ld(matrix(1, 3), matrix(0, 3), rep(-Inf, 3), 5, tau),
    tolerance = 1e-10
  )
  expect_equal(stan("out_far"), rep(-log(2 * pi), 3), tolerance = 1e-7)

  # the mode-centred scheme at a huge kappa, past the kappa ~1e16 where the
  # start's rate + 1 - c rounds to 0 unless (1 - c) is grouped first
  for (kappa_mode in c(1e17, 1e100)) {
    tau_mode <- .circmix_J(kappa_mode) / c(2, 5, 100)
    stan_mode <- circmix_stan_values(utils::modifyList(data, list(kappa_far = kappa_mode, tau_far = tau_mode)))
    expect_equal(
      stan_mode("out_far"),
      .circmix_vp_ld(matrix(1, 3), matrix(0, 3), rep(-Inf, 3), kappa_mode, tau_mode),
      tolerance = 1e-10
    )
  }
})

test_that("the Stan gradient of the variable-precision density is right below shape 2 and at a huge kappa", {
  skip_unless_stan_tests()
  skip_on_cran()

  # below kappa = 0.00447, J(kappa) < J0 / 10 and Stan Math's reverse-mode
  # gamma_p() drops its derivative as flat, which it is not at small shapes;
  # kappa 10 at shape 0.5 is an ordinary cell of the same branch. At kappa 1e15
  # the slope in log kappa was off by 4% to 11% and at 1e50 it was 0, in both
  # schemes, because the reverse-mode derivative of a Bessel-ratio J(kappa) is
  # erratic from kappa ~1e15 and 0 from 1e17.
  tab <- .circmix_kappa_table()
  cells <- rbind(
    expand.grid(y = c(0, 1, pi), shape = c(0.1, 1.5), kappa = c(0.0044, 0.001)),
    data.frame(y = c(0, 1), shape = 0.5, kappa = 10),
    expand.grid(y = c(0, 1), shape = c(0.5, 5), kappa = c(1e15, 1e50))
  )
  cells$tau <- .circmix_J(cells$kappa) / cells$shape
  n <- nrow(cells)
  code <- paste0(
    "functions {\n",
    read_lines2(system.file("stan_chunks", "circmix_functions.stan", package = "bmm")),
    "\n}\n",
    "data { int NK; vector[NK] logk; vector[NK] dlogk; real logJ_min; real dlogJ; int N; vector[N] y; }
    parameters { vector[N] log_kappa; vector[N] log_tau; }
    model {
      for (i in 1:N) {
        target += circmix_vp_ld([cos(y[i])]', [0.0]', negative_infinity(), exp(log_kappa[i]),
                                exp(log_tau[i]), 81, logk, dlogk, logJ_min, dlogJ);
      }
    }"
  )
  model <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(code))
  utils::capture.output(diagnosis <- model$diagnose(
    data = list(
      NK = length(tab$logkappa), logk = tab$logkappa, dlogk = tab$dlogkappa,
      logJ_min = tab$logJ_min, dlogJ = tab$dlogJ, N = n, y = cells$y
    ),
    init = list(list(log_kappa = log(cells$kappa), log_tau = log(cells$tau))),
    error = 1e6
  ))

  ld <- function(kappa, tau) {
    .circmix_vp_ld(matrix(cos(cells$y)), matrix(0, n), rep(-Inf, n), kappa, tau, 81L)
  }
  h <- 1e-5
  up <- exp(h)
  down <- exp(-h)
  slope <- c(
    ld(cells$kappa * up, cells$tau) - ld(cells$kappa * down, cells$tau),
    ld(cells$kappa, cells$tau * up) - ld(cells$kappa, cells$tau * down)
  ) / (2 * h)
  # relative, because diagnose() reports six significant digits
  expect_lt(max(abs(diagnosis$gradients()$model - slope) / pmax(abs(slope), 1e-3)), 1e-4)
})

test_that(".circmix_vp_ld() integrates below the shape where trigamma() is NaN, as Stan does", {
  # at shape 1e-200 nearly all precision mass sits at J ~ 0, so the response is uniform
  tau <- .circmix_J(10) / 1e-200
  expect_equal(.circmix_vp_ld(matrix(cos(1)), matrix(0), -Inf, 10, tau), -log(2 * pi), tolerance = 1e-10)
})
