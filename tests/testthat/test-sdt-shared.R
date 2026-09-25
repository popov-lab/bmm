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
              sdt_mafc, dsdt_mafc, rsdt_mafc)
  for (f in fns) {
    expect_equal(eval(formals(f)$dist), names(bmm:::.sdt_dists))
  }
})

test_that("the m-AFC probabilities agree with the registry's own cdf", {
  # P(correct) is the probability the signal variate beats m - 1 distractors,
  # int f(x - d) F(x)^(m - 1) dx, so a closed form is only right if it agrees
  # with the cdf the rest of the model uses. Taking the density off that same
  # cdf by central difference keeps the check independent of every branch of
  # .mafc_pc_r(): a swapped gumbel label and a self-consistent but mis-scaled
  # exponent both survive a comparison against the formula itself. gumbel_min
  # and gumbel_max coincide at m = 2, which is why the grid goes past it.
  # 1e-7 is the implementation's side: its 40-point Gauss-Hermite normal branch
  # sits 4.5e-8 from the integral at m = 8, d = 0, every other cell below 2e-10.
  pc_integral <- function(d, m, dist) {
    cdf <- bmm:::.sdt_dists[[dist]]$cdf
    dens <- function(x) (cdf(x + 1e-5) - cdf(x - 1e-5)) / 2e-5
    stats::integrate(function(x) dens(x - d) * cdf(x)^(m - 1),
                     -Inf, Inf, rel.tol = 1e-10)$value
  }
  for (dist in names(bmm:::.sdt_dists)) {
    for (m in c(2L, 3L, 5L, 8L)) {
      for (d in c(0, 0.8, 2)) {
        expect_equal(bmm:::.mafc_pc_r(d, m, dist), pc_integral(d, m, dist),
                     tolerance = 1e-7,
                     info = paste(dist, "m =", m, "d =", d))
      }
    }
  }
})

test_that("the default sdratio prior brackets the group-level estimates", {
  # two group-level SD ratios are on record: Mickes et al. (2007) Table 1
  # averages sd(lure)/sd(target) = 0.79 over their 13 retained subjects, i.e.
  # 1.26 signal over noise, and the broeder_schuetz_2009_e3 posterior gives
  # 1.46 [1.25, 1.71]. The prior has to reach past the widest ratio those
  # support and still keep its mass off ratios above 2, which recognition does
  # not produce -- one-sided, "covers the empirical range" cannot fail.
  main <- sdt_yn("n_old", "stimulus", "n_trials")$default_priors$sdratio$main
  sd <- as.numeric(sub("normal\\(0, ([0-9.]+)\\)", "\\1", main))
  expect_true(is.finite(sd), info = main)
  expect_gt(exp(qnorm(0.975, 0, sd)), 1.71)
  expect_gt(diff(pnorm(log(c(0.5, 2)), 0, sd)), 0.95)
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
# SHARED FORMULA CHECKS                                                    ####
############################################################################# !

sdt_stim_data <- function() {
  data.frame(
    n_old = c(10, 40, 15, 35),
    stimulus = c(0L, 1L, 0L, 1L),
    n_trials = c(50, 50, 50, 50),
    condition = c("A", "A", "B", "B"),
    id = c(1, 1, 2, 2)
  )
}

# check_model() must run first, or add_missing_parameters() overwrites a formula
# the user gave a fixed parameter: on a raw model `sdratio ~ stimulus` comes back
# as `sdratio ~ 1` and a direct check_formula() call tests nothing
sdt_check_formula <- function(formula, data = sdt_stim_data(),
                              model = sdt_yn("n_old", "stimulus",
                                             "n_trials")) {
  model <- check_model(model, data, formula)
  check_formula(model, data, formula)
}

test_that("a stimulus term on a parameter formula is refused", {
  # the likelihood already consumes stimulus to set each row's role, so the term
  # is confounded with the parameter's own intercept -- profiled flat to 6e-11
  expect_error(
    sdt_check_formula(bmf(d ~ stimulus, criterion ~ 1)),
    "'d'.*uses the stimulus variable 'stimulus'"
  )
  expect_error(
    sdt_check_formula(bmf(d ~ 1, criterion ~ stimulus)),
    "'criterion'.*uses the stimulus variable"
  )
  expect_error(
    sdt_check_formula(bmf(d ~ 1, criterion ~ 1, sdratio ~ stimulus)),
    "'sdratio'.*uses the stimulus variable"
  )
  # both offenders are named, not just the first
  expect_error(
    sdt_check_formula(bmf(d ~ stimulus, criterion ~ stimulus)),
    "'d', 'criterion'"
  )
})

test_that("a stimulus grouping factor is refused with the same message", {
  # rhs_vars() reports random-effect grouping variables, and a random intercept
  # per stimulus level is the same confound as a fixed effect
  expect_error(
    sdt_check_formula(bmf(d ~ 1 + (1 | stimulus), criterion ~ 1)),
    "as a predictor or as a grouping factor"
  )
})

test_that("the stimulus check leaves legitimate predictors alone", {
  expect_no_error(
    sdt_check_formula(bmf(d ~ condition, criterion ~ 1 + (1 | id)))
  )
  # neither the confound check nor the all-intercept sdratio heuristic applies
  # to a design with a real predictor
  expect_no_warning(
    sdt_check_formula(bmf(d ~ condition, criterion ~ 1 + (1 | id)))
  )
  # and it fires through the full pipeline, not only on a direct call
  expect_error(
    bmm(bmf(d ~ stimulus, criterion ~ 1), sdt_stim_data(),
        sdt_yn("n_old", "stimulus", "n_trials"),
        backend = "mock", mock_fit = 1, rename = FALSE),
    "uses the stimulus variable"
  )
})

test_that("a format without a stimulus other-var skips the confound check", {
  # sdt_mafc and sdt_ranking have no stimulus other-var, and check_formula.sdt
  # must fall through to NextMethod() instead of failing on a NULL stim_var
  model <- structure(
    list(name = "stub", parameters = list(par1 = ""),
         resp_vars = list(response = "n_old"),
         other_vars = list(stimulus = NULL)),
    class = c("bmmodel", "sdt", "stub")
  )
  formula <- bmf(par1 ~ 1)
  expect_identical(check_formula(model, sdt_stim_data(), formula), formula)
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
