test_that("mixture2p CD constructor works", {
  m <- mixture2p_cd(response = "resp", probe = "probe", target = "target")
  expect_s3_class(m, "bmmodel")
  expect_s3_class(m, "change_detection")
  expect_s3_class(m, "mixture2p")
  expect_s3_class(m, "mixture2p_cd")
  expect_false(inherits(m, "circular"))
  expect_equal(m$resp_vars$response, "resp")
  expect_equal(m$resp_vars$probe, "probe")
  expect_equal(m$resp_vars$target, "target")
  expect_true(m$void_mu)
  expect_equal(m$fixed_parameters$beta, 0)
})

test_that("mixture2p(task = 'cd') still works", {
  expect_warning(
    m <- mixture2p(response = "resp", probe = "probe", target = "target",
                   task = "cd"),
    "deprecated"
  )
  expect_s3_class(m, "mixture2p_cd")
})

test_that("mixture2p DE constructor still works", {
  m <- mixture2p(resp_error = "y")
  expect_s3_class(m, "circular")
  expect_s3_class(m, "mixture2p_de")
  expect_false(inherits(m, "change_detection"))
})

test_that("mixture2p CD validates required arguments", {
  expect_error(mixture2p_cd(probe = "p", target = "t"), "response")
  expect_error(mixture2p_cd(response = "r", target = "t"), "probe")
  expect_error(mixture2p_cd(response = "r", probe = "p"), "target")
})

test_that("check_data.change_detection validates data", {
  m <- mixture2p_cd(response = "resp", probe = "probe", target = "target")
  f <- bmf(kappa ~ 1, thetat ~ 1)

  dat <- data.frame(
    resp = c(0, 1, 1, 0, 1),
    probe = c(0.5, 1.0, -0.5, 2.0, -1.5),
    target = c(0.5, 0.5, 0.5, 0.5, 0.5)
  )

  checked <- check_data(m, dat, f)
  expect_true("probe_centered" %in% colnames(checked))
  expect_equal(checked$probe_centered, wrap(dat$probe - dat$target))
})

test_that("check_data.change_detection rejects invalid response", {
  m <- mixture2p_cd(response = "resp", probe = "probe", target = "target")
  f <- bmf(kappa ~ 1, thetat ~ 1)
  dat <- data.frame(resp = c(0, 2, 1), probe = c(0, 1, 2), target = c(0, 0, 0))
  expect_error(check_data(m, dat, f), "binary")
})

test_that("dmixture2p_cd returns valid probabilities", {
  probs <- dmixture2p_cd(c(0, 1), probe = 0, kappa = 5, thetat = 0.8)
  expect_length(probs, 2)
  expect_true(all(probs >= 0 & probs <= 1))
  expect_equal(sum(probs), 1, tolerance = 1e-4)
})

test_that("dmixture2p_cd accepts p_target alias", {
  probs <- dmixture2p_cd(c(0, 1), probe = 0, kappa = 5, p_target = 0.8)
  expect_length(probs, 2)
  expect_equal(sum(probs), 1, tolerance = 1e-4)
})

test_that("dmixture2p_cd: same probe has lower P(change)", {
  p_change_same <- dmixture2p_cd(1, probe = 0, kappa = 5, thetat = 0.8)
  p_change_diff <- dmixture2p_cd(1, probe = pi, kappa = 5, thetat = 0.8)
  expect_lt(p_change_same, p_change_diff)
})

test_that("rmixture2p_cd returns binary values", {
  set.seed(42)
  r <- rmixture2p_cd(100, probe = 1.0, kappa = 5, thetat = 0.8)
  expect_true(all(r %in% c(0, 1)))
})

test_that("mixture2p CD pipeline runs with mock backend", {
  dat <- data.frame(
    resp = c(rep(0, 50), rep(1, 50)),
    probe = runif(100, -pi, pi),
    target = runif(100, -pi, pi)
  )

  m <- mixture2p_cd(response = "resp", probe = "probe", target = "target")
  f <- bmf(kappa ~ 1, thetat ~ 1)

  mock_fit <- bmm(f, dat, m, backend = "mock", mock_fit = 1, rename = FALSE)
  expect_equal(mock_fit$fit, 1)
  expect_type(mock_fit$bmm, "list")
})

test_that("mixture2p_cd stancode uses compact probe data", {
  dat <- data.frame(resp = c(0, 1), probe = c(0.2, -0.1), target = c(0, 0))
  sc <- stancode(
    bmf(kappa ~ 1, thetat ~ 1),
    dat,
    model = mixture2p_cd(response = "resp", probe = "probe", target = "target")
  )

  expect_match(sc, "vector\\[2\\] probe_cd;", perl = TRUE)
  expect_no_match(sc, "vreal1")
})
