sdmcd_data <- function(n = 100, seed = 1) {
  withr::with_seed(seed, {
    dat <- data.frame(target = 0, probe = runif(n, -pi, pi))
    dat$resp <- rsdm_cd(n, dat$probe, c = 5, kappa = 4)
    dat
  })
}

test_that("sdm_cd is independent of the continuous reproduction model", {
  m <- sdm_cd(response = "resp", probe = "probe", target = "target")
  expect_s3_class(m, "change_detection")
  expect_s3_class(m, "sdm_cd")
  expect_false(inherits(m, "sdm"))
  expect_false(inherits(m, "circular"))
  expect_equal(m$fixed_parameters, list(mu = 0, criterion = 0))
})

test_that("the continuous reproduction sdm keeps its own classes and methods", {
  de <- sdm(resp_error = "y")
  expect_s3_class(de, "sdm")
  expect_s3_class(de, "sdm_simple")
  expect_false(inherits(de, "change_detection"))
  # the methods develop dispatches on must still exist
  expect_true(!is.null(getS3method("configure_model", "sdm", optional = TRUE)))
  expect_true(!is.null(getS3method("postprocess_brm", "sdm", optional = TRUE)))
})

test_that("sdm() dispatches on the response arguments", {
  expect_s3_class(sdm(response = "r", probe = "p", target = "t"), "sdm_cd")
  expect_error(sdm(), "Provide either")
  expect_error(sdm(resp_error = "y", response = "r"), "not both")
  expect_error(sdm(response = "r", probe = "p"), "all of")
})

test_that("dsdm_cd returns a proper probability mass function", {
  p <- dsdm_cd(c(0, 1), probe = 0.7, c = 5, kappa = 4)
  expect_equal(sum(p), 1)
  expect_true(all(p >= 0 & p <= 1))
})

test_that("dsdm_cd rejects invalid arguments", {
  expect_error(dsdm_cd(1, 0, kappa = -1), "kappa")
  expect_error(dsdm_cd(1, 0, c = -1), "c must")
  expect_error(dsdm_cd(2, 0), "binary")
})

test_that("an identical probe is rejected far above chance", {
  expect_lt(dsdm_cd(1, probe = 0, c = 5, kappa = 4), 0.3)
  expect_gt(dsdm_cd(1, probe = pi, c = 5, kappa = 4), 0.9)
})

test_that("P(change) increases monotonically with probe distance", {
  p <- dsdm_cd(rep(1, 25), seq(0, pi, length.out = 25), c = 5, kappa = 4)
  expect_true(all(diff(p) >= -1e-9))
})

test_that("stronger memory improves change detection", {
  expect_gt(dsdm_cd(1, probe = pi / 2, c = 15, kappa = 8),
            dsdm_cd(1, probe = pi / 2, c = 1, kappa = 8))
  expect_lt(dsdm_cd(1, probe = 0, c = 15, kappa = 8),
            dsdm_cd(1, probe = 0, c = 1, kappa = 8))
})

test_that("the adaptive quadrature holds up where a fixed rule fails", {
  # the sdm density is far more peaked than a von Mises; check the integral
  # against a peak-split adaptive reference across the sampled range
  reference <- function(lo, hi, c, kappa) {
    eta <- function(u) c * exp(kappa * (cos(u) - 1)) * sqrt(kappa / (2 * pi))
    peak <- eta(0)
    pieces <- unique(sort(c(lo, 0, hi)))
    tot <- sum(vapply(seq_len(length(pieces) - 1L), function(i) {
      stats::integrate(function(x) exp(eta(x) - peak), pieces[i], pieces[i + 1L],
                       rel.tol = 1e-12, subdivisions = 5000)$value
    }, numeric(1)))
    peak + log(tot)
  }
  grid <- expand.grid(c = c(1, 10, 40), kappa = c(1, 20, 100), hw = c(0.2, 1.5))
  err <- mapply(function(cc, kk, hh) {
    abs(.sdm_log_int(-hh, hh, 0, cc, kk) - reference(-hh, hh, cc, kk))
  }, grid$c, grid$kappa, grid$hw)
  expect_lt(max(err), 1e-8)
})

test_that("the normalising constant is rotation invariant", {
  expect_equal(.sdm_log_int(-pi, pi, 0, 6, 12),
               .sdm_log_int(-pi + 1.3, pi + 1.3, 1.3, 6, 12))
})

test_that("the criterion shifts responding in the expected direction", {
  p <- vapply(c(-0.5, 0, 0.5), function(b) {
    dsdm_cd(1, probe = 0.5, c = 4, kappa = 6, criterion = b)
  }, numeric(1))
  expect_true(all(diff(p) <= 0))
})

test_that("rsdm_cd recovers the density it was generated from", {
  p <- withr::with_seed(4, mean(rsdm_cd(1e5, probe = 1.0, c = 5, kappa = 4)))
  expect_equal(p, dsdm_cd(1, 1.0, c = 5, kappa = 4), tolerance = 0.02)
})

test_that("the sdm_cd pipeline runs and wires its Stan data", {
  dat <- sdmcd_data()
  m <- sdm_cd(response = "resp", probe = "probe", target = "target")
  f <- bmf(c ~ 1, kappa ~ 1)
  expect_equal(bmm(f, dat, m, backend = "mock", mock_fit = 1, rename = FALSE)$fit, 1)

  sc <- stancode(f, dat, model = m)
  expect_match(sc, "sdm_cd_lpmf")
  expect_match(sc, "array\\[N\\] real vreal1;")
})
