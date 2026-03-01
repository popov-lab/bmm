test_that("imm full CD constructor works", {
  m <- imm(response = "resp", probe = "probe", target = "target",
           nt_features = paste0("nt", 1:2, "_loc"),
           nt_distances = paste0("nt", 1:2, "_dist"),
           set_size = 3, task = "cd")
  expect_s3_class(m, "bmmodel")
  expect_s3_class(m, "change_detection")
  expect_s3_class(m, "non_targets")
  expect_s3_class(m, "imm")
  expect_s3_class(m, "imm_full")
  expect_s3_class(m, "imm_full_cd")
  expect_false(inherits(m, "circular"))
  expect_equal(m$resp_vars$response, "resp")
  expect_true(m$void_mu)
  expect_equal(m$fixed_parameters$beta, 0)
})

test_that("imm DE constructor still works", {
  m <- imm(resp_error = "y", nt_features = paste0("nt", 1:2),
           nt_distances = paste0("d", 1:2), set_size = 3)
  expect_s3_class(m, "circular")
  expect_s3_class(m, "imm_full_de")
  expect_false(inherits(m, "change_detection"))
})

test_that("imm CD validates required arguments", {
  expect_error(imm(task = "cd"), "response")
  expect_error(
    imm(response = "r", probe = "p", target = "t",
        nt_features = paste0("nt", 1:2), set_size = 3, task = "cd"),
    "nt_distances"
  )
})

test_that("dimm_cd returns valid probabilities", {
  p <- dimm_cd(1, probe = 0, nt_features = c(1, -1),
               nt_distances = c(0.5, 1.5), lure_idx = c(1, 1),
               kappa = 5, c = 2, a = 0.5, s = 2)
  expect_true(p >= 0 && p <= 1)
})

test_that("dimm_cd: same probe has lower P(change)", {
  p_same <- dimm_cd(1, probe = 0, nt_features = c(1), nt_distances = c(1),
                    lure_idx = c(1), kappa = 5, c = 2, a = 0.5, s = 2)
  p_diff <- dimm_cd(1, probe = pi, nt_features = c(1), nt_distances = c(1),
                    lure_idx = c(1), kappa = 5, c = 2, a = 0.5, s = 2)
  expect_lt(p_same, p_diff)
})

test_that("rimm_cd returns binary values", {
  set.seed(42)
  r <- rimm_cd(100, probe = 1.0, nt_features = c(1, -1),
               nt_distances = c(0.5, 1.5), lure_idx = c(1, 1),
               kappa = 5, c = 2, a = 0.5, s = 2)
  expect_true(all(r %in% c(0, 1)))
})

test_that("imm full CD pipeline runs with mock backend", {
  set.seed(42)
  dat <- data.frame(
    resp = c(rep(0, 50), rep(1, 50)),
    probe = runif(100, -pi, pi),
    target = runif(100, -pi, pi),
    nt1_loc = runif(100, -pi, pi),
    nt2_loc = runif(100, -pi, pi),
    nt1_dist = runif(100, 0, 3),
    nt2_dist = runif(100, 0, 3),
    ss = 3
  )

  m <- imm(response = "resp", probe = "probe", target = "target",
           nt_features = c("nt1_loc", "nt2_loc"),
           nt_distances = c("nt1_dist", "nt2_dist"),
           set_size = "ss", task = "cd")
  f <- bmf(kappa ~ 1, c ~ 1, a ~ 1, s ~ 1)

  mock_fit <- bmm(f, dat, m, backend = "mock", mock_fit = 1, rename = FALSE)
  expect_equal(mock_fit$fit, 1)
  expect_type(mock_fit$bmm, "list")
})
