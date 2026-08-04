cd_data <- function(n = 100, seed = 1) {
  withr::with_seed(seed, {
    dat <- data.frame(target = runif(n, -pi, pi))
    dat$probe <- wrap(dat$target + sample(c(rep(0, n / 2), runif(n / 2, -pi, pi))))
    dat$resp <- rmixture2p_cd(n, wrap(dat$probe - dat$target), kappa = 6, p_mem = 0.7)
    dat
  })
}

test_that("mixture2p_cd is independent of the continuous reproduction model", {
  m <- mixture2p_cd(response = "resp", probe = "probe", target = "target")
  expect_s3_class(m, "bmmodel")
  expect_s3_class(m, "change_detection")
  expect_s3_class(m, "mixture2p_cd")
  expect_false(inherits(m, "mixture2p"))
  expect_false(inherits(m, "circular"))
  expect_equal(m$resp_vars, list(response = "resp", probe = "probe", target = "target"))
})

test_that("mixture2p_cd fixes mu and criterion but leaves them estimable", {
  m <- mixture2p_cd(response = "resp", probe = "probe", target = "target")
  expect_equal(m$fixed_parameters, list(mu = 0, criterion = 0))

  # a formula for the criterion removes it from the fixed parameters
  freed <- check_model(m, cd_data(), bmf(kappa ~ 1, thetat ~ 1, criterion ~ 1))
  expect_false("criterion" %in% names(freed$fixed_parameters))
  still_fixed <- check_model(m, cd_data(), bmf(kappa ~ 1, thetat ~ 1))
  expect_true("criterion" %in% names(still_fixed$fixed_parameters))
})

test_that("mixture2p() dispatches on the response arguments", {
  expect_s3_class(mixture2p(resp_error = "y"), "mixture2p")
  expect_false(inherits(mixture2p(resp_error = "y"), "change_detection"))

  cd <- mixture2p(response = "resp", probe = "probe", target = "target")
  expect_s3_class(cd, "mixture2p_cd")

  expect_error(mixture2p(), "Provide either")
  expect_error(mixture2p(resp_error = "y", response = "r"), "not both")
  expect_error(mixture2p(response = "r", probe = "p"), "all of")
})

test_that("mixture2p_de and mixture2p agree for continuous reproduction", {
  expect_equal(class(mixture2p_de("y")), class(mixture2p(resp_error = "y")))
  expect_error(mixture2p_de(), "resp_error")
})

test_that("check_data.change_detection centres the probe on the target", {
  m <- mixture2p_cd(response = "resp", probe = "probe", target = "target")
  dat <- data.frame(resp = c(0, 1, 1), probe = c(0.5, 1, -0.5), target = 0.5)
  checked <- check_data(m, dat, bmf(kappa ~ 1, thetat ~ 1))
  expect_equal(checked$probe_centered, wrap(dat$probe - dat$target))
})

test_that("check_data.change_detection rejects a non-binary response", {
  m <- mixture2p_cd(response = "resp", probe = "probe", target = "target")
  dat <- data.frame(resp = c(0, 2, 1), probe = c(0, 1, 2), target = 0)
  expect_error(check_data(m, dat, bmf(kappa ~ 1, thetat ~ 1)), "binary")
})

test_that("dmixture2p_cd returns a proper probability mass function", {
  p <- dmixture2p_cd(c(0, 1), probe = 0.7, kappa = 5, p_mem = 0.8)
  expect_length(p, 2)
  expect_true(all(p >= 0 & p <= 1))
  expect_equal(sum(p), 1)
})

test_that("dmixture2p_cd rejects invalid arguments", {
  expect_error(dmixture2p_cd(1, 0, kappa = -1), "kappa")
  expect_error(dmixture2p_cd(1, 0, p_mem = 1.2), "p_mem")
  expect_error(dmixture2p_cd(2, 0), "binary")
})

test_that("an identical probe is rejected far above chance", {
  # the implemented decision rule must not pin same-trial accuracy at 0.5
  p_change <- dmixture2p_cd(1, probe = 0, kappa = 8, p_mem = 0.8)
  expect_lt(p_change, 0.35)
  expect_gt(dmixture2p_cd(1, probe = pi, kappa = 8, p_mem = 0.8), p_change)
})

test_that("P(change) increases monotonically with probe distance", {
  probes <- seq(0, pi, length.out = 25)
  p <- dmixture2p_cd(rep(1, 25), probes, kappa = 6, p_mem = 0.75)
  expect_true(all(diff(p) >= -1e-12))
})

test_that("better memory improves change detection", {
  args <- list(response = 1, probe = pi / 2, p_mem = 0.8)
  expect_gt(
    do.call(dmixture2p_cd, c(args, kappa = 20)),
    do.call(dmixture2p_cd, c(args, kappa = 3))
  )
  # with no memory at all, responding is at the guessing rate set by the criterion
  expect_equal(
    dmixture2p_cd(1, probe = 0, kappa = 6, p_mem = 0),
    1 - .cd_crit_angle(6) / pi
  )
})

test_that("the criterion shifts responding in the expected direction", {
  p <- vapply(c(-1, 0, 1), function(b) {
    dmixture2p_cd(1, probe = 1, kappa = 6, p_mem = 0.7, criterion = b)
  }, numeric(1))
  expect_true(all(diff(p) < 0))
})

test_that("the decision boundary is free of p_mem at an unbiased criterion", {
  # Lin & Oberauer (2022), Appendix B
  expect_equal(.cd_crit_angle(7, criterion = 0, p_mem = 0.2),
               .cd_crit_angle(7, criterion = 0, p_mem = 0.9))
  expect_equal(.cd_crit_angle(7, criterion = 0, p_mem = 1),
               .cd_crit_angle(7, criterion = 1e-12, p_mem = 1), tolerance = 1e-6)
})

test_that("the arc integral matches a dense reference rule", {
  reference <- function(probe, kappa, p_mem) {
    d <- .cd_crit_angle(kappa)
    x <- seq(probe - d, probe + d, length.out = 200001)
    dens <- exp(kappa * cos(x) - .cd_log_i0(kappa)) / (2 * pi)
    1 - (p_mem * sum((dens[-1] + dens[-length(dens)]) / 2) * (x[2] - x[1]) +
      (1 - p_mem) * d / pi)
  }
  grid <- expand.grid(probe = c(0, 0.8, 2.2), kappa = c(1, 5, 25), p_mem = c(0.3, 0.9))
  got <- dmixture2p_cd(rep(1, nrow(grid)), grid$probe, grid$kappa, grid$p_mem)
  want <- mapply(reference, grid$probe, grid$kappa, grid$p_mem)
  expect_equal(got, want, tolerance = 1e-8)
})

test_that("rmixture2p_cd recovers the density it was generated from", {
  p <- withr::with_seed(3, mean(rmixture2p_cd(2e5, probe = 1.2, kappa = 6, p_mem = 0.7)))
  expect_equal(p, dmixture2p_cd(1, 1.2, kappa = 6, p_mem = 0.7), tolerance = 0.01)
})

test_that("the mixture2p_cd pipeline runs and emits compact Stan data", {
  dat <- cd_data()
  m <- mixture2p_cd(response = "resp", probe = "probe", target = "target")
  fit <- bmm(bmf(kappa ~ 1, thetat ~ 1), dat, m,
             backend = "mock", mock_fit = 1, rename = FALSE)
  expect_equal(fit$fit, 1)

  sc <- stancode(bmf(kappa ~ 1, thetat ~ 1), dat, model = m)
  expect_match(sc, "mixture2p_cd_lpmf")
  expect_match(sc, "array\\[N\\] real vreal1;")
})

test_that("the Stan criterion code is only emitted when it is estimated", {
  dat <- cd_data()
  m <- mixture2p_cd(response = "resp", probe = "probe", target = "target")
  fixed <- standata(bmf(kappa ~ 1, thetat ~ 1), dat, model = m)
  free <- standata(bmf(kappa ~ 1, thetat ~ 1, criterion ~ 1), dat, model = m)
  expect_equal(as.integer(fixed$cd_free_criterion), 0L)
  expect_equal(as.integer(free$cd_free_criterion), 1L)
})
