test_that("sdm simple CD constructor works", {
  m <- sdm(response = "resp", probe = "probe", target = "target",
           task = "cd")
  expect_s3_class(m, "bmmodel")
  expect_s3_class(m, "change_detection")
  expect_s3_class(m, "sdm")
  expect_s3_class(m, "sdm_simple")
  expect_s3_class(m, "sdm_simple_cd")
  expect_false(inherits(m, "circular"))
  expect_equal(m$resp_vars$response, "resp")
  expect_true(m$void_mu)
  expect_equal(m$fixed_parameters$beta, 0)
})

test_that("sdm DE constructor still works", {
  m <- sdm(resp_error = "y")
  expect_s3_class(m, "circular")
  expect_s3_class(m, "sdm_simple_de")
  expect_false(inherits(m, "change_detection"))
})

test_that("sdm CD validates required arguments", {
  expect_error(sdm(task = "cd"), "response")
  expect_error(sdm(response = "r", task = "cd"), "probe")
  expect_error(sdm(response = "r", probe = "p", task = "cd"), "target")
})

test_that("dsdm_cd returns valid probabilities", {
  p <- dsdm_cd(1, probe = 0, c = 4, kappa = 3)
  expect_true(p >= 0 && p <= 1)
})

test_that("dsdm_cd: same probe has lower P(change)", {
  p_same <- dsdm_cd(1, probe = 0, c = 4, kappa = 3)
  p_diff <- dsdm_cd(1, probe = pi, c = 4, kappa = 3)
  expect_lt(p_same, p_diff)
})

test_that("rsdm_cd returns binary values", {
  set.seed(42)
  r <- rsdm_cd(100, probe = 1.0, c = 4, kappa = 3)
  expect_true(all(r %in% c(0, 1)))
})

test_that("sdm simple CD pipeline runs with mock backend", {
  dat <- data.frame(
    resp = c(rep(0, 50), rep(1, 50)),
    probe = runif(100, -pi, pi),
    target = runif(100, -pi, pi)
  )

  m <- sdm(response = "resp", probe = "probe", target = "target",
           task = "cd")
  f <- bmf(c ~ 1, kappa ~ 1)

  mock_fit <- bmm(f, dat, m, backend = "mock", mock_fit = 1, rename = FALSE)
  expect_equal(mock_fit$fit, 1)
  expect_type(mock_fit$bmm, "list")
})
