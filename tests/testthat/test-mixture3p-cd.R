test_that("mixture3p CD constructor works", {
  m <- mixture3p(response = "resp", probe = "probe", target = "target",
                 nt_features = paste0("nt", 1:2, "_loc"), set_size = 3,
                 task = "cd")
  expect_s3_class(m, "bmmodel")
  expect_s3_class(m, "change_detection")
  expect_s3_class(m, "non_targets")
  expect_s3_class(m, "mixture3p")
  expect_s3_class(m, "mixture3p_cd")
  expect_false(inherits(m, "circular"))
  expect_equal(m$resp_vars$response, "resp")
  expect_true(m$void_mu)
  expect_equal(m$fixed_parameters$beta, 0)
})

test_that("mixture3p DE constructor still works", {
  m <- mixture3p(resp_error = "y", nt_features = paste0("nt", 1:2), set_size = 3)
  expect_s3_class(m, "circular")
  expect_s3_class(m, "mixture3p_de")
  expect_false(inherits(m, "change_detection"))
})

test_that("mixture3p CD validates required arguments", {
  expect_error(mixture3p(task = "cd"), "response")
  expect_error(
    mixture3p(response = "r", probe = "p", target = "t", task = "cd"),
    "nt_features"
  )
})

test_that("dmixture3p_cd returns valid probabilities", {
  p <- dmixture3p_cd(1, probe = 0, nt_features = c(1, -1),
                     lure_idx = c(1, 1), kappa = 5, thetat = log(6),
                     thetant = log(3))
  expect_true(p >= 0 && p <= 1)
})

test_that("mixture3p CD pipeline runs with mock backend", {
  set.seed(42)
  dat <- data.frame(
    resp = c(rep(0, 50), rep(1, 50)),
    probe = runif(100, -pi, pi),
    target = runif(100, -pi, pi),
    nt1_loc = runif(100, -pi, pi),
    nt2_loc = runif(100, -pi, pi),
    ss = 3
  )

  m <- mixture3p(response = "resp", probe = "probe", target = "target",
                 nt_features = c("nt1_loc", "nt2_loc"), set_size = "ss",
                 task = "cd")
  f <- bmf(kappa ~ 1, thetat ~ 1, thetant ~ 1)

  mock_fit <- bmm(f, dat, m, backend = "mock", mock_fit = 1, rename = FALSE)
  expect_equal(mock_fit$fit, 1)
  expect_type(mock_fit$bmm, "list")
})
