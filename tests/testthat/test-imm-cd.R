immcd_data <- function(n = 100, seed = 1) {
  withr::with_seed(seed, {
    dat <- data.frame(
      target = 0,
      probe = runif(n, -pi, pi),
      nt1 = runif(n, -pi, pi),
      nt2 = runif(n, -pi, pi),
      d1 = runif(n, 0.2, 2),
      d2 = runif(n, 0.2, 2),
      ss = 3
    )
    dat$resp <- rimm_cd(
      n, dat$probe,
      nt_features = cbind(dat$nt1, dat$nt2),
      nt_distances = cbind(dat$d1, dat$d2),
      kappa = 6, c = 5, a = 0.6, s = 1.35
    )
    dat
  })
}

immcd_model <- function() {
  imm_cd(
    response = "resp", probe = "probe", target = "target",
    nt_features = c("nt1", "nt2"), nt_distances = c("d1", "d2"), set_size = "ss"
  )
}

test_that("imm_cd is independent of the continuous reproduction model", {
  m <- immcd_model()
  expect_s3_class(m, "change_detection")
  expect_s3_class(m, "non_targets")
  expect_s3_class(m, "imm_cd")
  expect_s3_class(m, "imm_cd_full")
  expect_false(inherits(m, "imm"))
  expect_false(inherits(m, "circular"))
  expect_equal(m$fixed_parameters, list(mu = 0, criterion = 0))
})

test_that("imm_cd only supports the full version for now", {
  expect_error(
    imm_cd(response = "r", probe = "p", target = "t", nt_features = "nt1",
           nt_distances = "d1", set_size = 2, version = "abc"),
    "Only version"
  )
})

test_that("imm() dispatches on the response arguments", {
  de <- imm(resp_error = "y", nt_features = "nt1", nt_distances = "d1", set_size = 2)
  expect_s3_class(de, "imm")
  expect_s3_class(de, "imm_full")
  expect_false(inherits(de, "change_detection"))

  cd <- imm(response = "r", probe = "p", target = "t", nt_features = "nt1",
            nt_distances = "d1", set_size = 2)
  expect_s3_class(cd, "imm_cd_full")

  expect_error(imm(nt_features = "nt1", set_size = 2), "Provide either")
})

test_that("dimm_cd returns a proper probability mass function", {
  p <- dimm_cd(c(0, 1), probe = 0.7, nt_features = c(1.5, -2),
               nt_distances = c(0.5, 1.2), kappa = 6, c = 5, a = 0.6, s = 1.35)
  expect_equal(sum(p), 1)
  expect_true(all(p >= 0 & p <= 1))
})

test_that("an identical probe is rejected far above chance", {
  args <- list(nt_features = c(2, -2.5), nt_distances = c(1.5, 2),
               kappa = 8, c = 20, a = 0.3, s = 1.5)
  same <- do.call(dimm_cd, c(list(response = 1, probe = 0), args))
  far <- do.call(dimm_cd, c(list(response = 1, probe = pi), args))
  expect_lt(same, 0.3)
  expect_gt(far, 0.9)
})

test_that("nearer non-targets interfere more", {
  args <- list(response = 1, probe = 2.0, nt_features = c(2.0, -2.5),
               kappa = 8, c = 6, a = 0.4, s = 1.2)
  near <- do.call(dimm_cd, c(args, list(nt_distances = c(0.1, 2))))
  far <- do.call(dimm_cd, c(args, list(nt_distances = c(3, 2))))
  # a nearby non-target matching the probe is harder to call a change
  expect_lt(near, far)
})

test_that("stronger context activation improves change detection", {
  args <- list(response = 1, probe = pi / 2, nt_features = c(2, -2.5),
               nt_distances = c(1, 1.5), kappa = 8, a = 0.4, s = 1.2)
  expect_gt(do.call(dimm_cd, c(args, c = 20)), do.call(dimm_cd, c(args, c = 0.5)))
})

test_that("inactive lures drop out of the retrieval distribution", {
  args <- list(response = 1, probe = 1.3, kappa = 6, c = 5, a = 0.6, s = 1.35)
  with_lure <- do.call(dimm_cd, c(args, list(
    nt_features = c(1.3, 2.0), nt_distances = c(0.5, 1), lure_idx = c(1, 0)
  )))
  only_lure <- do.call(dimm_cd, c(args, list(
    nt_features = 1.3, nt_distances = 0.5, lure_idx = 1
  )))
  expect_equal(with_lure, only_lure)
})

test_that("the imm_cd pipeline runs and wires its Stan data", {
  dat <- immcd_data()
  m <- immcd_model()
  f <- bmf(kappa ~ 1, c ~ 1, a ~ 1, s ~ 1)
  expect_equal(bmm(f, dat, m, backend = "mock", mock_fit = 1, rename = FALSE)$fit, 1)

  sd_ <- standata(f, dat, model = m)
  expect_equal(dim(sd_$cd_nt_features), c(nrow(dat), 2L))
  expect_equal(dim(sd_$cd_nt_distances), c(nrow(dat), 2L))
  expect_equal(as.integer(sd_$cd_free_criterion), 0L)

  free <- standata(bmf(kappa ~ 1, c ~ 1, a ~ 1, s ~ 1, criterion ~ 1),
                   dat, model = m)
  expect_equal(as.integer(free$cd_free_criterion), 1L)
})
