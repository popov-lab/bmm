# Test Meta-d' SDT model specification and integration

############################################################################# !
# MODEL CONSTRUCTOR TESTS                                                ####
############################################################################# !

test_that("sdt_metad model can be created", {
  model <- sdt_metad(c("r1", "r2", "r3", "r4"), "stimulus")
  expect_s3_class(model, "bmmodel")
  expect_s3_class(model, "sdt")
  expect_s3_class(model, "sdt_metad")
})

test_that("sdt_metad has metad parameter", {
  model <- sdt_metad(c("r1", "r2", "r3", "r4"), "stimulus")
  expect_true("metad" %in% names(model$parameters))
  expect_true("metad" %in% names(model$links))
  expect_true("metad" %in% names(model$default_priors))
  expect_equal(model$links$metad, "identity")
})

test_that("sdt_metad has standard SDT parameters too", {
  model <- sdt_metad(c("r1", "r2", "r3", "r4"), "stimulus")
  expect_true("dprime" %in% names(model$parameters))
  expect_true("criterion" %in% names(model$parameters))
  expect_true("spacing" %in% names(model$parameters))
})

test_that("sdt_metad has sdratio as overridable fixed parameter", {
  model <- sdt_metad(c("r1", "r2", "r3", "r4"), "stimulus")
  expect_true("sdratio" %in% names(model$parameters))
  expect_equal(model$fixed_parameters$sdratio, 0)
})

test_that("sdt_metad works with different distributions", {
  for (d in c("normal", "logistic", "gumbel_min", "gumbel_max")) {
    model <- sdt_metad(c("r1", "r2", "r3", "r4"), "stimulus", dist = d)
    expect_s3_class(model, "sdt_metad")
  }
})


############################################################################# !
# DISTRIBUTION FUNCTION TESTS                                            ####
############################################################################# !

test_that("rsdt_metad generates meta-d' data", {
  dat <- rsdt_metad(n_per_cell = 100, n_subjects = 5, dprime = 1.5, criterion = 0,
                    metad = 1.0, n_ratings = 4, spacing = 0.5)
  expect_true(is.data.frame(dat))
  expect_true(all(c("r1", "r2", "r3", "r4", "stimulus", "id") %in% names(dat)))
  expect_equal(nrow(dat), 10)
  expect_true(all(rowSums(dat[, paste0("r", 1:4)]) == 100))
})

test_that("rsdt_metad with metad=dprime matches standard rating SDT", {
  set.seed(42)
  dat_metad <- rsdt_metad(n_per_cell = 1000, n_subjects = 1, dprime = 1.5,
                          criterion = 0, metad = 1.5, n_ratings = 4, spacing = 0.5)
  set.seed(42)
  dat_std <- rsdt_rating(n_per_cell = 1000, n_subjects = 1, dprime = 1.5,
                         criterion = 0, n_ratings = 4, spacing = 0.5)
  expect_equal(dat_metad[, paste0("r", 1:4)], dat_std[, paste0("r", 1:4)])
})

test_that("rsdt_metad requires metad argument", {
  expect_error(
    rsdt_metad(n_per_cell = 50, n_subjects = 1, dprime = 1.5, criterion = 0,
               n_ratings = 4, spacing = 0.5),
    "metad.*required"
  )
})

test_that("dsdt_metad computes density", {
  thresholds <- c(-0.5, 0.0, 0.5)

  d1 <- dsdt_metad(counts = c(5, 15, 25, 55), stimulus = 1,
                   dprime = 1.5, metad = 1.0, thresholds = thresholds)
  expect_true(d1 > 0)

  d2 <- dsdt_metad(counts = c(55, 25, 15, 5), stimulus = 0,
                   dprime = 1.5, metad = 1.0, thresholds = thresholds)
  expect_true(d2 > 0)
})

test_that("dsdt_metad with metad=dprime matches standard rating", {
  thresholds <- c(-0.5, 0.0, 0.5)
  counts <- c(10, 20, 30, 40)

  d_metad <- dsdt_metad(counts = counts, stimulus = 1,
                        dprime = 1.5, metad = 1.5, thresholds = thresholds,
                        log = TRUE)
  d_std <- dsdt_rating(counts = counts, stimulus = 1,
                       dprime = 1.5, thresholds = thresholds, log = TRUE)
  expect_equal(d_metad, d_std, tolerance = 1e-10)
})


############################################################################# !
# META-D' CATEGORY PROBABILITY TESTS                                     ####
############################################################################# !

test_that("meta-d' category probs sum to 1", {
  thresholds <- c(-0.5, 0.0, 0.5)

  probs_old <- bmm:::.sdt_metad_category_probs(
    thresholds, dprime = 1.5, metad = 1.0, stimulus = 1,
    sdratio = 1, dist = "normal"
  )
  expect_equal(sum(probs_old), 1, tolerance = 1e-10)

  probs_new <- bmm:::.sdt_metad_category_probs(
    thresholds, dprime = 1.5, metad = 1.0, stimulus = 0,
    sdratio = 1, dist = "normal"
  )
  expect_equal(sum(probs_new), 1, tolerance = 1e-10)
})

test_that("meta-d' with metad=dprime gives same probs as standard SDT", {
  thresholds <- c(-0.5, 0.0, 0.5)

  probs_metad <- bmm:::.sdt_metad_category_probs(
    thresholds, dprime = 1.5, metad = 1.5, stimulus = 1,
    sdratio = 1, dist = "normal"
  )
  probs_std <- bmm:::.sdt_category_probs(
    thresholds, shift = 1.5 / 2, sdratio = 1, stimulus = 1, dist = "normal"
  )
  expect_equal(probs_metad, probs_std, tolerance = 1e-10)
})

test_that("meta-d' lower than dprime shifts confidence toward criterion", {
  thresholds <- c(-0.5, 0.0, 0.5)

  probs_ideal <- bmm:::.sdt_metad_category_probs(
    thresholds, dprime = 2.0, metad = 2.0, stimulus = 1,
    sdratio = 1, dist = "normal"
  )
  probs_poor <- bmm:::.sdt_metad_category_probs(
    thresholds, dprime = 2.0, metad = 0.5, stimulus = 1,
    sdratio = 1, dist = "normal"
  )

  K <- length(probs_ideal)
  expect_true(probs_poor[K] < probs_ideal[K])
})

test_that("meta-d' works with all 4 distributions", {
  thresholds <- c(-0.5, 0.0, 0.5)

  for (d in c("normal", "logistic", "gumbel_min", "gumbel_max")) {
    probs <- bmm:::.sdt_metad_category_probs(
      thresholds, dprime = 1.5, metad = 1.0, stimulus = 1,
      sdratio = 1, dist = d
    )
    expect_equal(sum(probs), 1, tolerance = 1e-10)
    expect_true(all(probs > 0))
  }
})


############################################################################# !
# FORMULA CONSTRUCTION TESTS                                              ####
############################################################################# !

test_that("sdt_metad produces valid stancode with equidistant thresholds", {
  dat <- rsdt_metad(n_per_cell = 50, n_subjects = 3, dprime = 1.5, criterion = 0,
                    metad = 1.0, n_ratings = 4, spacing = 0.5)
  model <- sdt_metad(c("r1", "r2", "r3", "r4"), "stimulus")
  formula <- bmf(dprime ~ 1, metad ~ 1, criterion ~ 1, spacing ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
  expect_true(grepl("sdt_metad_lpmf", code, fixed = TRUE))
  expect_true(grepl("metad", code))
})

test_that("sdt_metad produces valid stancode with K=6", {
  dat <- rsdt_metad(n_per_cell = 50, n_subjects = 2, dprime = 1.5, criterion = 0,
                    metad = 1.0, n_ratings = 6, spacing = 0.3)
  model <- sdt_metad(c("r1", "r2", "r3", "r4", "r5", "r6"), "stimulus")
  formula <- bmf(dprime ~ 1, metad ~ 1, criterion ~ 1, spacing ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
})

test_that("sdt_metad produces valid stancode with predictors on metad", {
  dat <- rsdt_metad(n_per_cell = 50, n_subjects = 3, dprime = 1.5, criterion = 0,
                    metad = 1.0, n_ratings = 4, spacing = 0.5)
  dat$condition <- rep(c("A", "B"), length.out = nrow(dat))
  model <- sdt_metad(c("r1", "r2", "r3", "r4"), "stimulus")
  formula <- bmf(dprime ~ 1, metad ~ condition, criterion ~ 1, spacing ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
})

test_that("sdt_metad with all 4 distributions produces valid stancode", {
  dat <- rsdt_metad(n_per_cell = 50, n_subjects = 2, dprime = 1.5, criterion = 0,
                    metad = 1.0, n_ratings = 4, spacing = 0.5)
  for (d in c("normal", "logistic", "gumbel_min", "gumbel_max")) {
    model <- sdt_metad(c("r1", "r2", "r3", "r4"), "stimulus", dist = d)
    formula <- bmf(dprime ~ 1, metad ~ 1, criterion ~ 1, spacing ~ 1)
    code <- stancode(formula, data = dat, model = model)
    expect_true(nchar(code) > 0)
  }
})

test_that("sdt_metad with log_distance thresholds produces valid stancode", {
  dat <- rsdt_metad(n_per_cell = 50, n_subjects = 3, dprime = 1.5, criterion = 0,
                    metad = 1.0, n_ratings = 4, spacing = 0.5)
  model <- sdt_metad(c("r1", "r2", "r3", "r4"), "stimulus",
                     threshold_type = "log_distance")
  formula <- bmf(dprime ~ 1, metad ~ 1, criterion ~ 1,
                 delta1 ~ 1, delta3 ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
})
