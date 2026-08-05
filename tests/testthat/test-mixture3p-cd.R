cd3p_data <- function(n = 100, seed = 1) {
  withr::with_seed(seed, {
    dat <- data.frame(
      target = 0,
      probe = runif(n, -pi, pi),
      nt1 = runif(n, -pi, pi),
      nt2 = runif(n, -pi, pi),
      ss = 3
    )
    dat$resp <- rmixture3p_cd(
      n, dat$probe, nt_features = cbind(dat$nt1, dat$nt2),
      kappa = 6, thetat = 1.2, thetant = 0.1
    )
    dat
  })
}

cd3p_model <- function() {
  mixture3p_cd(
    response = "resp", probe = "probe", target = "target",
    nt_features = c("nt1", "nt2"), set_size = "ss"
  )
}

test_that("mixture3p_cd is independent of the continuous reproduction model", {
  m <- cd3p_model()
  expect_s3_class(m, "change_detection")
  expect_s3_class(m, "non_targets")
  expect_s3_class(m, "mixture3p_cd")
  expect_false(inherits(m, "mixture3p"))
  expect_false(inherits(m, "circular"))
  expect_equal(m$fixed_parameters, list(mu = 0, criterion = 0))
})

test_that("mixture3p() dispatches on the response arguments", {
  de <- mixture3p(resp_error = "y", nt_features = c("nt1", "nt2"), set_size = "ss")
  expect_s3_class(de, "mixture3p")
  expect_false(inherits(de, "change_detection"))

  cd <- mixture3p(response = "r", probe = "p", target = "t",
                  nt_features = c("nt1", "nt2"), set_size = "ss")
  expect_s3_class(cd, "mixture3p_cd")

  expect_error(mixture3p(nt_features = "nt1", set_size = 2), "Provide either")
  expect_error(
    mixture3p(resp_error = "y", response = "r", nt_features = "nt1", set_size = 2),
    "not both"
  )
})

test_that("dmixture3p_cd returns a proper probability mass function", {
  p <- dmixture3p_cd(c(0, 1), probe = 0.7, nt_features = c(1.5, -2),
                     kappa = 5, thetat = 1, thetant = 0.2)
  expect_equal(sum(p), 1)
  expect_true(all(p >= 0 & p <= 1))
})

test_that("an identical probe is rejected far above chance", {
  p <- dmixture3p_cd(1, probe = 0, nt_features = c(2, -2.5), kappa = 8,
                     thetat = 1.5, thetant = -0.5)
  expect_lt(p, 0.35)
})

test_that("mixture3p_cd predicts an intrusion cost", {
  # a probe matching a non-target is harder to reject than a new probe at the
  # same distance from the target (Lin & Oberauer, 2022)
  nt <- 2.0
  intrusion <- dmixture3p_cd(1, probe = nt, nt_features = c(nt, -0.9),
                             kappa = 8, thetat = 1, thetant = 0.5)
  new_probe <- dmixture3p_cd(1, probe = -nt, nt_features = c(nt, -0.9),
                             kappa = 8, thetat = 1, thetant = 0.5)
  expect_lt(intrusion, new_probe)
})

test_that("inactive lures drop out of the retrieval distribution", {
  args <- list(response = 1, probe = 1.3, kappa = 6, thetat = 1, thetant = 0.4)
  with_lure <- do.call(dmixture3p_cd, c(args, list(
    nt_features = c(1.3, 2.0), lure_idx = c(1, 0)
  )))
  only_lure <- do.call(dmixture3p_cd, c(args, list(
    nt_features = 1.3, lure_idx = 1
  )))
  expect_equal(with_lure, only_lure)
})

test_that("mixture3p_cd matches mixture2p_cd when non-targets carry no weight", {
  p3 <- dmixture3p_cd(1, probe = 1.1, nt_features = c(2, -2), lure_idx = c(0, 0),
                      kappa = 6, thetat = 0, thetant = -100)
  p2 <- dmixture2p_cd(1, probe = 1.1, kappa = 6, p_mem = 0.5)
  expect_equal(p3, p2, tolerance = 1e-8)
})

test_that("the mixture3p_cd pipeline runs and wires its Stan data", {
  dat <- cd3p_data()
  m <- cd3p_model()
  f <- bmf(kappa ~ 1, thetat ~ 1, thetant ~ 1)
  expect_equal(bmm(f, dat, m, backend = "mock", mock_fit = 1, rename = FALSE)$fit, 1)

  sd_ <- standata(f, dat, model = m)
  expect_equal(dim(sd_$cd_nt_features), c(nrow(dat), 2L))
  expect_equal(dim(sd_$cd_lure_idx), c(nrow(dat), 2L))
  expect_equal(as.integer(sd_$cd_free_criterion), 0L)

  free <- standata(bmf(kappa ~ 1, thetat ~ 1, thetant ~ 1, criterion ~ 1),
                   dat, model = m)
  expect_equal(as.integer(free$cd_free_criterion), 1L)
})
