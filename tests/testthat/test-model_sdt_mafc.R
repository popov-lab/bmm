# Test m-AFC SDT model specification and integration

############################################################################# !
# MODEL CONSTRUCTOR TESTS                                                ####
############################################################################# !

test_that("sdt_mafc model can be created", {
  model <- sdt_mafc("n_correct", n_trials = "n_trials", m = 4)
  expect_s3_class(model, "bmmodel")
  expect_s3_class(model, "sdt")
  expect_s3_class(model, "sdt_mafc")
})

test_that("sdt_mafc has correct parameters (only dprime, no criterion)", {
  model <- sdt_mafc("n_correct", n_trials = "n_trials", m = 4)
  expect_true("dprime" %in% names(model$parameters))
  expect_true("mu" %in% names(model$parameters))
  expect_false("criterion" %in% names(model$parameters))
  expect_equal(model$fixed_parameters$mu, 0)
})

test_that("sdt_mafc stores m correctly", {
  model <- sdt_mafc("n_correct", n_trials = "n_trials", m = 4)
  expect_equal(model$other_vars$m, 4L)

  model2 <- sdt_mafc("n_correct", n_trials = "n_trials", m = 8)
  expect_equal(model2$other_vars$m, 8L)
})

test_that("sdt_mafc has correct links", {
  model <- sdt_mafc("n_correct", n_trials = "n_trials", m = 4)
  expect_equal(model$links$mu, "identity")
  expect_equal(model$links$dprime, "identity")
})

test_that("sdt_mafc accepts custom links", {
  model <- sdt_mafc("n_correct", n_trials = "n_trials", m = 4,
                    links = list(dprime = "log"))
  expect_equal(model$links$dprime, "log")
})

test_that("sdt_mafc requires m >= 2", {
  expect_error(sdt_mafc("n_correct", n_trials = "n_trials", m = 1))
})

test_that("sdt_mafc requires dist = 'normal'", {
  expect_error(
    sdt_mafc("n_correct", n_trials = "n_trials",
             m = 4, dist = "logistic")
  )
})

test_that("sdt_mafc does not have a stimulus argument", {
  expect_silent(
    sdt_mafc("n_correct", n_trials = "n_trials", m = 4)
  )
})

test_that("sdt_mafc has init_ranges with mu and dprime", {
  model <- sdt_mafc("n_correct", n_trials = "n_trials", m = 4)
  expect_false(is.null(model$init_ranges))
  expect_true(all(c("mu", "dprime") %in% names(model$init_ranges)))
  expect_equal(model$init_ranges$dprime, c(0.5, 1.5))
})


############################################################################# !
# DISTRIBUTION FUNCTION TESTS                                            ####
############################################################################# !

test_that(".mafc_pc_r gives correct 2-AFC result", {
  d <- 1.5
  expected <- pnorm(d / sqrt(2))
  result <- bmm:::.mafc_pc_r(d, 2L)
  expect_equal(result, expected, tolerance = 1e-10)
})

test_that(".mafc_pc_r gives correct results for m >= 3", {
  d <- 1.5
  exact <- integrate(function(x) dnorm(x) * pnorm(x + d)^3, -Inf, Inf)$value
  result <- bmm:::.mafc_pc_r(d, 4L)
  expect_equal(result, exact, tolerance = 1e-6)
})

test_that(".mafc_pc_r works for various m values", {
  d <- 1.0
  pc2 <- bmm:::.mafc_pc_r(d, 2L)
  pc4 <- bmm:::.mafc_pc_r(d, 4L)
  pc8 <- bmm:::.mafc_pc_r(d, 8L)
  expect_true(pc2 > pc4)
  expect_true(pc4 > pc8)

  expect_true(pc2 > 0.5)
  expect_true(pc4 > 0.25)
  expect_true(pc8 > 0.125)
})

test_that(".mafc_pc_r at d'=0 gives chance performance", {
  expect_equal(bmm:::.mafc_pc_r(0, 2L), 0.5, tolerance = 1e-10)
  expect_equal(bmm:::.mafc_pc_r(0, 4L), 0.25, tolerance = 1e-6)
  expect_equal(bmm:::.mafc_pc_r(0, 8L), 0.125, tolerance = 1e-3)
})

test_that("rsdt_mafc generates valid mafc data", {
  set.seed(42)
  dat <- rsdt_mafc(n_per_cell = 100, n_subjects = 5, dprime = 1.5, m = 4)

  expect_true("id" %in% colnames(dat))
  expect_true("n_trials" %in% colnames(dat))
  expect_true("n_correct" %in% colnames(dat))
  expect_equal(nrow(dat), 5)
  expect_true(all(dat$n_correct >= 0))
  expect_true(all(dat$n_correct <= dat$n_trials))
  expect_true(all(dat$n_trials == 100))
})

test_that("dsdt_mafc computes valid mafc densities", {
  dens <- dsdt_mafc(n_correct = 80, n_trials = 100, dprime = 1.5, m = 4)
  expect_true(dens > 0 && dens <= 1)

  log_dens <- dsdt_mafc(n_correct = 80, n_trials = 100, dprime = 1.5, m = 4,
                        log = TRUE)
  expect_equal(log(dens), log_dens, tolerance = 1e-10)
})

test_that("dsdt_mafc validates inputs", {
  expect_error(dsdt_mafc(n_correct = NULL, n_trials = 100, dprime = 1.5, m = 4))
  expect_error(dsdt_mafc(n_correct = 80, n_trials = NULL, dprime = 1.5, m = 4))
  expect_error(dsdt_mafc(n_correct = 80, n_trials = 100, dprime = 1.5, m = NULL))
  expect_error(dsdt_mafc(n_correct = 110, n_trials = 100, dprime = 1.5, m = 4))
})

test_that("rsdt_mafc validates inputs", {
  expect_error(rsdt_mafc(n_per_cell = 50, n_subjects = 1, dprime = 1, m = NULL))
  expect_error(rsdt_mafc(n_per_cell = 50, n_subjects = 1, dprime = 1, m = 1))
})


############################################################################# !
# CHECK_DATA TESTS                                                       ####
############################################################################# !

test_that("check_data validates mafc data correctly", {
  model <- sdt_mafc("n_correct", n_trials = "n_trials", m = 4)
  dat <- rsdt_mafc(n_per_cell = 50, n_subjects = 3, dprime = 1.5, m = 4)
  formula <- bmf(dprime ~ 1)

  result <- check_data(model, dat, formula)
  expect_true("m_afc" %in% colnames(result))
  expect_equal(unique(result$m_afc), 4L)
})

test_that("check_data errors on missing columns for mafc", {
  model <- sdt_mafc("n_correct", n_trials = "n_trials", m = 4)
  dat <- data.frame(id = 1:5, wrong_col = 50)
  formula <- bmf(dprime ~ 1)
  expect_error(check_data(model, dat, formula))
})

test_that("check_data errors on invalid response counts for mafc", {
  model <- sdt_mafc("n_correct", n_trials = "n_trials", m = 4)
  dat <- data.frame(n_correct = c(-1, 10), n_trials = c(50, 50))
  formula <- bmf(dprime ~ 1)
  expect_error(check_data(model, dat, formula))
})

test_that("check_data errors when n_correct > n_trials for mafc", {
  model <- sdt_mafc("n_correct", n_trials = "n_trials", m = 4)
  dat <- data.frame(n_correct = c(60, 10), n_trials = c(50, 50))
  formula <- bmf(dprime ~ 1)
  expect_error(check_data(model, dat, formula))
})


############################################################################# !
# FORMULA CONSTRUCTION TESTS                                              ####
############################################################################# !

test_that("bmf2bf.sdt_mafc produces correct formula for mafc", {
  model <- sdt_mafc("n_correct", n_trials = "n_trials", m = 4)
  dat <- rsdt_mafc(n_per_cell = 50, n_subjects = 3, dprime = 1.5, m = 4)
  dat$m_afc <- 4L
  bf <- bmf2bf(model, bmf(dprime ~ 1))

  expect_s3_class(bf, "brmsformula")
  bf_str <- deparse(bf$formula, width.cutoff = 500)
  expect_true(grepl("n_correct", bf_str))
  expect_true(grepl("trials", bf_str))
  expect_true(grepl("vint.*m_afc", bf_str))
})


############################################################################# !
# STAN CODE GENERATION TESTS                                              ####
############################################################################# !

test_that("sdt_mafc produces valid stancode (2-AFC)", {
  dat <- rsdt_mafc(n_per_cell = 50, n_subjects = 3, dprime = 1.5, m = 2)
  model <- sdt_mafc("n_correct", n_trials = "n_trials", m = 2)
  formula <- bmf(dprime ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
  expect_true(grepl("mafc_pc", code))
  expect_true(grepl("sdt_mafc_lpmf", code))
})

test_that("sdt_mafc produces valid stancode (4-AFC)", {
  dat <- rsdt_mafc(n_per_cell = 50, n_subjects = 3, dprime = 1.5, m = 4)
  model <- sdt_mafc("n_correct", n_trials = "n_trials", m = 4)
  formula <- bmf(dprime ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
  expect_true(grepl("gh_nodes", code))
  expect_true(grepl("gh_weights", code))
})

test_that("sdt_mafc produces valid stancode (8-AFC)", {
  dat <- rsdt_mafc(n_per_cell = 50, n_subjects = 3, dprime = 1.5, m = 8)
  model <- sdt_mafc("n_correct", n_trials = "n_trials", m = 8)
  formula <- bmf(dprime ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
})

test_that("sdt_mafc stancode works with predictors on dprime", {
  set.seed(42)
  dat <- rsdt_mafc(n_per_cell = 50, n_subjects = 5, dprime = 1.5, m = 4)
  dat$condition <- rep(c("A", "B"), length.out = nrow(dat))
  model <- sdt_mafc("n_correct", n_trials = "n_trials", m = 4)
  formula <- bmf(dprime ~ condition)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
})

test_that("sdt_mafc stancode works with random effects on dprime", {
  dat <- rsdt_mafc(n_per_cell = 50, n_subjects = 5, dprime = 1.5, m = 4)
  model <- sdt_mafc("n_correct", n_trials = "n_trials", m = 4)
  formula <- bmf(dprime ~ 1 + (1 | id))
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
})
