# =============================================================================
# Numerics tests for the RDM start-point kernels (Stage 4 of #353): continuity
# of .dwald_full()/.pwald_full() across the midpoint switch, and far-tail
# finiteness of drdm()/prdm(). See R/distributions.R for the derivation.
# =============================================================================

# .dwald_full()/.pwald_full() fall back to the plain-Wald midpoint
# approximation below A* = .rdm_midpoint_ratio * s * sqrt(t); a fit's gradient
# would see a kink at that boundary if the two branches disagreed in value or
# slope for a start-point range that is small but not exactly zero. Central
# differences (eps = 1e-6) are taken just inside and just outside the switch
# and compared to each other, not to an analytic gradient, so the test is
# blind to a shared error in both branches but not to a discontinuity at the
# seam itself.
test_that("the midpoint switch in .dwald_full()/.pwald_full() is continuous in value and gradient", {
  s <- 1.1
  t <- 0.6
  drift <- 2.3
  gap <- 0.9
  a_star <- .rdm_midpoint_ratio * s * sqrt(t)
  a_lo <- a_star * (1 - 1e-6)
  a_hi <- a_star * (1 + 1e-6)

  f_d <- function(t, drift, gap, A, s) {
    .dwald_full(t, drift = drift, bound = gap + A, A = A, s = s, log = TRUE)
  }
  f_p <- function(t, drift, gap, A, s) {
    .pwald_full(t, drift = drift, bound = gap + A, A = A, s = s,
               lower.tail = FALSE, log.p = TRUE)
  }

  expect_lt(abs(f_d(t, drift, gap, a_hi, s) - f_d(t, drift, gap, a_lo, s)), 1e-8)
  expect_lt(abs(f_p(t, drift, gap, a_hi, s) - f_p(t, drift, gap, a_lo, s)), 1e-8)

  eps <- 1e-6
  central_diff <- function(fn, par, A) {
    args <- list(t = t, drift = drift, gap = gap, A = A, s = s)
    plus <- args
    plus[[par]] <- plus[[par]] + eps
    minus <- args
    minus[[par]] <- minus[[par]] - eps
    (do.call(fn, plus) - do.call(fn, minus)) / (2 * eps)
  }
  for (par in c("t", "drift", "gap", "s")) {
    expect_equal(
      central_diff(f_d, par, a_hi), central_diff(f_d, par, a_lo),
      tolerance = 1e-4
    )
    expect_equal(
      central_diff(f_p, par, a_hi), central_diff(f_p, par, a_lo),
      tolerance = 1e-4
    )
  }
})

# The raw (non-log-space) Wald forms underflow to -Inf/NaN once the density or
# survivor gets this far into the tail; drdm()/prdm() must stay finite at all
# three sp regimes (fixed at 0, numerically tiny but nonzero, and a magnitude
# a fit would actually estimate).
test_that("drdm() and prdm() stay finite far into the tail", {
  for (sp in c(0, 1e-8, 0.3)) {
    ld <- drdm(20, 1, drift = 8, gap = 0.5, ndt = 0, s = 0.5, sp = sp, log = TRUE)
    ls <- prdm(20, drift = 8, gap = 0.5, ndt = 0, s = 0.5, sp = sp,
              lower.tail = FALSE, log.p = TRUE)
    expect_true(is.finite(ld))
    expect_true(is.finite(ls))
    expect_lt(ld, -700)
    expect_lt(ls, -700)
  }
})
