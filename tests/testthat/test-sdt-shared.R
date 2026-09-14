# Test SDT shared infrastructure

############################################################################# !
# CDF HELPER TESTS                                                       ####
############################################################################# !

test_that(".sdt_cdf returns correct values for known inputs", {
  # Normal: Phi(0) = 0.5
  expect_equal(bmm:::.sdt_cdf(0, "normal"), 0.5)
  # Logistic: inv_logit(0) = 0.5
  expect_equal(bmm:::.sdt_cdf(0, "logistic"), 0.5)
  # Gumbel_min (smallest extreme value) at 0: 1 - exp(-exp(0)) = 1 - exp(-1)
  expect_equal(bmm:::.sdt_cdf(0, "gumbel_min"), 1 - exp(-1), tolerance = 1e-10)
  # Gumbel_max (largest extreme value) at 0: exp(-exp(0)) = exp(-1)
  expect_equal(bmm:::.sdt_cdf(0, "gumbel_max"), exp(-1), tolerance = 1e-10)
})

test_that("gumbel labels follow the extreme-value convention", {
  # the min-type is the mirror image of the max-type about zero
  eta <- c(-2, -0.5, 0, 0.5, 2)
  expect_equal(bmm:::.sdt_cdf(eta, "gumbel_min"),
               1 - bmm:::.sdt_cdf(-eta, "gumbel_max"))
  # gumbel_max is the log-log / evd::pgumbel parameterisation
  expect_equal(bmm:::.sdt_cdf(eta, "gumbel_max"), exp(-exp(-eta)))
})

test_that("every dist argument offers exactly the registry's distributions", {
  # the registry order defines the dist_type integer passed to Stan, so a
  # signature that drifts out of step with it becomes an off-by-one
  fns <- list(sdt_yn, dsdt_yn, rsdt_yn, sdt_d, sdt_criterion,
              sdt_mafc, dsdt_mafc, rsdt_mafc,
              sdt_rating, dsdt_rating, rsdt_rating)
  for (f in fns) {
    expect_equal(eval(formals(f)$dist), names(bmm:::.sdt_dists))
  }
})

test_that("a two-category rating reduces to the binary likelihood", {
  # both models must implement the same decision rule -- evidence above the
  # criterion is an "old" response -- which for the asymmetric extreme-value
  # distributions is only true if neither evaluates its cdf on the mirrored axis
  for (di in names(bmm:::.sdt_dists)) {
    for (stim in c(0L, 1L)) {
      expect_equal(
        bmm:::.sdt_category_probs(-0.2, 1.5, 1.3, stim, di)[2],
        dsdt_yn(1, 1, stim, d = 1.5, criterion = -0.2,
                sdratio = 1.3, dist = di),
        info = paste(di, "stimulus", stim)
      )
    }
  }
})

test_that("the m-AFC closed forms sit on the right extreme-value branch", {
  # taking the max of largest-extreme-value variates is what yields the softmax;
  # the smallest-extreme-value case is the Gamma ratio. Swapping these fits the
  # mirror model, and the two agree at m = 2, so check m > 2.
  for (m in c(4L, 8L)) {
    expect_equal(bmm:::.mafc_pc_r(1.2, m, "gumbel_max"),
                 1 / (1 + (m - 1) * exp(-1.2)), info = paste("m =", m))
    expect_equal(
      bmm:::.mafc_pc_r(1.2, m, "gumbel_min"),
      exp(lgamma(1 + exp(-1.2)) + lgamma(m) - lgamma(m + exp(-1.2))),
      info = paste("m =", m)
    )
    expect_false(isTRUE(all.equal(bmm:::.mafc_pc_r(1.2, m, "gumbel_min"),
                                  bmm:::.mafc_pc_r(1.2, m, "gumbel_max"))))
  }
})

test_that("the sdratio prior stays inside the calibrated quadrature range", {
  # .ranking_gh_n() is calibrated over sdratio in [0.5, 2.0]. If the default
  # prior is ever widened past that, the sampler will visit ratios the
  # quadrature was never verified for, so the two must be changed together.
  models <- list(sdt_yn("n_old", "stimulus", "n_trials"),
                 sdt_ranking(paste0("rank", 1:4), m = 4, dist = "normal"))
  for (model in models) {
    main <- model$default_priors$sdratio$main
    sd <- as.numeric(sub("normal\\(0, ([0-9.]+)\\)", "\\1", main))
    expect_true(is.finite(sd), info = main)
    covered <- diff(pnorm(log(c(0.5, 2)), 0, sd))
    expect_gt(covered, 0.95)
  }
})

test_that("sdratio priors agree across the SDT models", {
  yn <- sdt_yn("n_old", "stimulus", "n_trials")$default_priors$sdratio
  rk <- sdt_ranking(paste0("rank", 1:4), m = 4,
                    dist = "normal")$default_priors$sdratio
  expect_identical(yn, rk)
})

test_that("quantile functions invert their cdfs", {
  p <- c(0.05, 0.25, 0.5, 0.75, 0.95)
  for (d in names(bmm:::.sdt_dists)) {
    expect_equal(bmm:::.sdt_cdf(bmm:::.sdt_dists[[d]]$qf(p), d), p,
                 tolerance = 1e-10, info = d)
  }
})


############################################################################# !
# LOG-SCALE CDF TESTS                                                    ####
############################################################################# !

test_that("lcdf and lccdf agree with the probability scale where it is exact", {
  # only the central range: further out the naive log(cdf) is the inaccurate
  # side of the comparison, which is the whole reason the log-scale pair exists
  eta <- seq(-2, 2, by = 0.5)
  for (d in names(bmm:::.sdt_dists)) {
    p <- bmm:::.sdt_cdf(eta, d)
    expect_equal(bmm:::.sdt_log_cdf(eta, d), log(p), tolerance = 1e-10, info = d)
    expect_equal(bmm:::.sdt_log_ccdf(eta, d), log(1 - p),
                 tolerance = 1e-10, info = d)
  }
})

test_that("lcdf and lccdf stay finite where the probability scale underflows", {
  # the naive path returns log(0) = -Inf here; Stan and the log-scale path do not
  for (d in names(bmm:::.sdt_dists)) {
    expect_true(is.finite(bmm:::.sdt_log_cdf(-40, d)), info = d)
    expect_true(is.finite(bmm:::.sdt_log_ccdf(40, d)), info = d)
  }
  expect_equal(log(bmm:::.sdt_cdf(-40, "normal")), -Inf)
})

test_that(".sdt_cdf is vectorized over eta", {
  eta <- c(-1, 0, 1)
  result <- bmm:::.sdt_cdf(eta, "normal")
  expect_length(result, 3)
  expect_equal(result, pnorm(eta))
})


############################################################################# !
# DPRIME AND CRITERION TESTS                                              ####
############################################################################# !

test_that("sdt_d computes correct values for normal distribution", {
  hr <- pnorm(1)
  far <- pnorm(-1)
  expect_equal(sdt_d(hr, far, dist = "normal"), 2, tolerance = 1e-10)
})

test_that("sdt_criterion computes correct values for normal distribution", {
  hr <- pnorm(1)
  far <- pnorm(-1)
  expect_equal(sdt_criterion(hr, far, dist = "normal"), 0, tolerance = 1e-10)
})

test_that("sdt_d validates input", {
  expect_error(sdt_d(0, 0.5), "between 0 and 1")
  expect_error(sdt_d(1, 0.5), "between 0 and 1")
  expect_error(sdt_d(0.5, 0), "between 0 and 1")
})


############################################################################# !
# VECTORIZATION TESTS                                                      ####
############################################################################# !

test_that(".sdt_eta is vectorized over all arguments", {
  eta <- bmm:::.sdt_eta(1.5, 0.2, c(0, 1))
  expect_length(eta, 2)
  expect_equal(eta[1], -1.5/2 - 0.2)
  expect_equal(eta[2], 1.5/2 - 0.2)

  eta <- bmm:::.sdt_eta(c(1.0, 2.0), 0.2, c(0, 1))
  expect_length(eta, 2)
  expect_equal(eta[1], -1.0/2 - 0.2)
  expect_equal(eta[2], 2.0/2 - 0.2)

  # `d` is d_a, so the separation in noise units is d * sqrt((1 + r^2) / 2)
  s <- sqrt((1 + c(1.0, 1.3)^2) / 2)
  eta <- bmm:::.sdt_eta(1.5, 0.2, c(0, 1), sdratio = c(1.0, 1.3))
  expect_length(eta, 2)
  expect_equal(eta[1], -1.5 * s[1] / 2 - 0.2)
  expect_equal(eta[2], (1.5 * s[2] / 2 - 0.2) / 1.3)
})

test_that(".sdt_eta divides only the signal trials by sdratio", {
  # the noise trial carries no scale division, so its eta must match an
  # equal-variance call with the separation widened to d * s
  s <- sqrt((1 + 1.5^2) / 2)
  eta_uv <- bmm:::.sdt_eta(1.5, 0.2, 0, sdratio = 1.5)
  eta_ev <- bmm:::.sdt_eta(1.5 * s, 0.2, 0, sdratio = 1)
  expect_equal(eta_uv, eta_ev)
})
