# Test DPSDT (Dual Process SDT) model specification and integration

############################################################################# !
# MODEL CONSTRUCTOR TESTS                                                ####
############################################################################# !

test_that("sdt_dp model can be created", {
  model <- sdt_dp(c("r1", "r2", "r3", "r4"), "stimulus")
  expect_s3_class(model, "bmmodel")
  expect_s3_class(model, "sdt")
  expect_s3_class(model, "sdt_dp")
})

test_that("sdt_dp has Ro and Rn parameters", {
  model <- sdt_dp(c("r1", "r2", "r3", "r4"), "stimulus")
  expect_true("Ro" %in% names(model$parameters))
  expect_true("Ro" %in% names(model$links))
  expect_true("Ro" %in% names(model$default_priors))
  expect_equal(model$links$Ro, "identity")
  expect_true("Rn" %in% names(model$parameters))
  expect_true("Rn" %in% names(model$links))
  expect_true("Rn" %in% names(model$default_priors))
  expect_equal(model$links$Rn, "identity")
})

test_that("sdt_dp has standard SDT parameters too", {
  model <- sdt_dp(c("r1", "r2", "r3", "r4"), "stimulus")
  expect_true("dprime" %in% names(model$parameters))
  expect_true("criterion" %in% names(model$parameters))
  expect_true("spacing" %in% names(model$parameters))
})

test_that("sdt_dp has sdratio as overridable fixed parameter", {
  model <- sdt_dp(c("r1", "r2", "r3", "r4"), "stimulus")
  expect_true("sdratio" %in% names(model$parameters))
  expect_equal(model$fixed_parameters$sdratio, 0)
})

test_that("sdt_dp works with different distributions", {
  for (d in c("normal", "logistic", "gumbel_min", "gumbel_max")) {
    model <- sdt_dp(c("r1", "r2", "r3", "r4"), "stimulus", dist = d)
    expect_s3_class(model, "sdt_dp")
    expect_equal(model$other_vars$dist, d)
  }
})

test_that("sdt_dp works with log_distance thresholds", {
  model <- sdt_dp(c("r1", "r2", "r3", "r4"), "stimulus",
                  threshold_type = "log_distance")
  expect_true("delta1" %in% names(model$parameters))
  expect_true("delta3" %in% names(model$parameters))
  expect_false("spacing" %in% names(model$parameters))
})

test_that("sdt_dp model has NULL init_ranges", {
  model <- sdt_dp(paste0("r", 1:6), "stimulus")
  expect_null(model$init_ranges)
})


############################################################################# !
# DISTRIBUTION FUNCTION TESTS                                            ####
############################################################################# !

test_that("rsdt_dp generates DPSDT data", {
  dat <- rsdt_dp(n_per_cell = 100, n_subjects = 5, dprime = 1.0, criterion = 0,
                 Ro = 0.3, Rn = 0, n_ratings = 4, spacing = 0.5)
  expect_true(is.data.frame(dat))
  expect_true(all(c("r1", "r2", "r3", "r4", "stimulus", "id") %in% names(dat)))
  expect_equal(nrow(dat), 10)
  expect_true(all(rowSums(dat[, paste0("r", 1:4)]) == 100))
})

test_that("rsdt_dp Ro=0 Rn=0 matches standard rating SDT", {
  set.seed(42)
  dat_dp <- rsdt_dp(n_per_cell = 1000, n_subjects = 1, dprime = 1.5, criterion = 0,
                    Ro = 0, Rn = 0, n_ratings = 4, spacing = 0.5)
  set.seed(42)
  dat_std <- rsdt_rating(n_per_cell = 1000, n_subjects = 1, dprime = 1.5,
                         criterion = 0, n_ratings = 4, spacing = 0.5)
  expect_equal(dat_dp[, paste0("r", 1:4)], dat_std[, paste0("r", 1:4)])
})

test_that("rsdt_dp Ro shifts old items to highest confidence", {
  set.seed(42)
  dat_no_r <- rsdt_dp(n_per_cell = 1000, n_subjects = 1, dprime = 1.0, criterion = 0,
                      Ro = 0, Rn = 0, n_ratings = 4, spacing = 0.5)
  set.seed(42)
  dat_hi_r <- rsdt_dp(n_per_cell = 1000, n_subjects = 1, dprime = 1.0, criterion = 0,
                      Ro = 0.5, Rn = 0, n_ratings = 4, spacing = 0.5)

  old_no_r <- dat_no_r[dat_no_r$stimulus == 1, ]
  old_hi_r <- dat_hi_r[dat_hi_r$stimulus == 1, ]
  expect_true(old_hi_r$r4 > old_no_r$r4)
})

test_that("rsdt_dp requires Ro and Rn arguments", {
  expect_error(
    rsdt_dp(n_per_cell = 50, n_subjects = 1, dprime = 1.5, criterion = 0,
            Rn = 0, n_ratings = 4, spacing = 0.5),
    "Ro.*required"
  )
  expect_error(
    rsdt_dp(n_per_cell = 50, n_subjects = 1, dprime = 1.5, criterion = 0,
            Ro = 0, n_ratings = 4, spacing = 0.5),
    "Rn.*required"
  )
})

test_that("rsdt_dp validates Ro and Rn bounds", {
  expect_error(
    rsdt_dp(n_per_cell = 50, n_subjects = 1, dprime = 1.5, criterion = 0,
            Ro = -0.1, Rn = 0, n_ratings = 4, spacing = 0.5),
    "Ro must be between 0 and 1"
  )
  expect_error(
    rsdt_dp(n_per_cell = 50, n_subjects = 1, dprime = 1.5, criterion = 0,
            Ro = 1.1, Rn = 0, n_ratings = 4, spacing = 0.5),
    "Ro must be between 0 and 1"
  )
  expect_error(
    rsdt_dp(n_per_cell = 50, n_subjects = 1, dprime = 1.5, criterion = 0,
            Ro = 0, Rn = -0.1, n_ratings = 4, spacing = 0.5),
    "Rn must be between 0 and 1"
  )
  expect_error(
    rsdt_dp(n_per_cell = 50, n_subjects = 1, dprime = 1.5, criterion = 0,
            Ro = 0, Rn = 1.1, n_ratings = 4, spacing = 0.5),
    "Rn must be between 0 and 1"
  )
})

test_that("rsdt_dp Rn shifts new items to lowest confidence", {
  set.seed(42)
  dat_no_r <- rsdt_dp(n_per_cell = 1000, n_subjects = 1, dprime = 1.0, criterion = 0,
                      Ro = 0, Rn = 0, n_ratings = 4, spacing = 0.5)
  set.seed(42)
  dat_hi_r <- rsdt_dp(n_per_cell = 1000, n_subjects = 1, dprime = 1.0, criterion = 0,
                      Ro = 0, Rn = 0.5, n_ratings = 4, spacing = 0.5)

  new_no_r <- dat_no_r[dat_no_r$stimulus == 0, ]
  new_hi_r <- dat_hi_r[dat_hi_r$stimulus == 0, ]
  expect_true(new_hi_r$r1 > new_no_r$r1)
})

test_that("dsdt_dp computes density", {
  thresholds <- c(-0.5, 0.0, 0.5)

  d1 <- dsdt_dp(counts = c(5, 15, 25, 55), stimulus = 1,
                dprime = 1.5, thresholds = thresholds,
                Ro = 0.3, Rn = 0)
  expect_true(d1 > 0)

  d2 <- dsdt_dp(counts = c(55, 25, 15, 5), stimulus = 0,
                dprime = 1.5, thresholds = thresholds,
                Ro = 0.3, Rn = 0.2)
  expect_true(d2 > 0)
})

test_that("dsdt_dp Ro=0 Rn=0 matches standard rating density", {
  thresholds <- c(-0.5, 0.0, 0.5)
  counts <- c(10, 20, 30, 40)

  d_dp <- dsdt_dp(counts = counts, stimulus = 1,
                  dprime = 1.5, thresholds = thresholds,
                  Ro = 0, Rn = 0, log = TRUE)
  d_std <- dsdt_rating(counts = counts, stimulus = 1,
                       dprime = 1.5, thresholds = thresholds, log = TRUE)
  expect_equal(d_dp, d_std)
})

test_that("dsdt_dp Ro=0 Rn=0 matches for noise trials too", {
  thresholds <- c(-0.5, 0.0, 0.5)
  counts <- c(40, 30, 20, 10)

  d_dp <- dsdt_dp(counts = counts, stimulus = 0,
                  dprime = 1.5, thresholds = thresholds,
                  Ro = 0.3, Rn = 0, log = TRUE)
  d_std <- dsdt_rating(counts = counts, stimulus = 0,
                       dprime = 1.5, thresholds = thresholds, log = TRUE)
  expect_equal(d_dp, d_std)
})


############################################################################# !
# DPSDT CATEGORY PROBABILITY TESTS                                       ####
############################################################################# !

test_that("DPSDT category probs sum to 1", {
  thresholds <- c(-0.5, 0.0, 0.5)

  probs_old <- bmm:::.sdt_dpsdt_category_probs(
    thresholds, shift = 0.75, sdratio = 1, stimulus = 1,
    dist = "normal", Ro = 0.3
  )
  expect_equal(sum(probs_old), 1, tolerance = 1e-10)

  probs_new <- bmm:::.sdt_dpsdt_category_probs(
    thresholds, shift = -0.75, sdratio = 1, stimulus = 0,
    dist = "normal", Ro = 0.3
  )
  expect_equal(sum(probs_new), 1, tolerance = 1e-10)
})

test_that("DPSDT Ro increases p[K] for old items", {
  thresholds <- c(-0.5, 0.0, 0.5)
  shift <- 0.75

  probs_no_r <- bmm:::.sdt_category_probs(
    thresholds, shift, sdratio = 1, stimulus = 1, dist = "normal"
  )
  probs_hi_r <- bmm:::.sdt_dpsdt_category_probs(
    thresholds, shift, sdratio = 1, stimulus = 1, dist = "normal", Ro = 0.5
  )

  K <- length(probs_no_r)
  expect_true(probs_hi_r[K] > probs_no_r[K])

  for (k in 1:(K - 1)) {
    expect_true(probs_hi_r[k] < probs_no_r[k])
  }
})

test_that("DPSDT Ro does not affect new items when Rn=0", {
  thresholds <- c(-0.5, 0.0, 0.5)
  shift <- -0.75

  probs_std <- bmm:::.sdt_category_probs(
    thresholds, shift, sdratio = 1, stimulus = 0, dist = "normal"
  )
  probs_dpsdt <- bmm:::.sdt_dpsdt_category_probs(
    thresholds, shift, sdratio = 1, stimulus = 0, dist = "normal",
    Ro = 0.5, Rn = 0
  )
  expect_equal(probs_dpsdt, probs_std, tolerance = 1e-10)
})

test_that("DPSDT with Ro=1 puts all old items in highest category", {
  thresholds <- c(-0.5, 0.0, 0.5)
  shift <- 0.5

  probs <- bmm:::.sdt_dpsdt_category_probs(
    thresholds, shift, sdratio = 1, stimulus = 1, dist = "normal", Ro = 1.0
  )

  K <- length(probs)
  expect_equal(probs[K], 1.0, tolerance = 1e-10)
  for (k in 1:(K - 1)) {
    expect_true(probs[k] < 1e-8)
  }
})

test_that("DPSDT Rn increases p[1] for new items", {
  thresholds <- c(-0.5, 0.0, 0.5)
  shift <- -0.75

  probs_no_r <- bmm:::.sdt_category_probs(
    thresholds, shift, sdratio = 1, stimulus = 0, dist = "normal"
  )
  probs_hi_r <- bmm:::.sdt_dpsdt_category_probs(
    thresholds, shift, sdratio = 1, stimulus = 0, dist = "normal",
    Ro = 0, Rn = 0.5
  )

  expect_true(probs_hi_r[1] > probs_no_r[1])

  K <- length(probs_no_r)
  for (k in 2:K) {
    expect_true(probs_hi_r[k] < probs_no_r[k])
  }
})

test_that("DPSDT Rn does not affect old items", {
  thresholds <- c(-0.5, 0.0, 0.5)
  shift <- 0.75

  probs_no_rn <- bmm:::.sdt_dpsdt_category_probs(
    thresholds, shift, sdratio = 1, stimulus = 1, dist = "normal",
    Ro = 0.3, Rn = 0
  )
  probs_hi_rn <- bmm:::.sdt_dpsdt_category_probs(
    thresholds, shift, sdratio = 1, stimulus = 1, dist = "normal",
    Ro = 0.3, Rn = 0.5
  )
  expect_equal(probs_hi_rn, probs_no_rn, tolerance = 1e-10)
})

test_that("DPSDT with Rn=1 puts all new items in lowest category", {
  thresholds <- c(-0.5, 0.0, 0.5)
  shift <- -0.5

  probs <- bmm:::.sdt_dpsdt_category_probs(
    thresholds, shift, sdratio = 1, stimulus = 0, dist = "normal",
    Ro = 0, Rn = 1.0
  )

  expect_equal(probs[1], 1.0, tolerance = 1e-10)
  K <- length(probs)
  for (k in 2:K) {
    expect_true(probs[k] < 1e-8)
  }
})

test_that("DPSDT category probs sum to 1 with both Ro and Rn", {
  thresholds <- c(-0.5, 0.0, 0.5)

  probs_old <- bmm:::.sdt_dpsdt_category_probs(
    thresholds, shift = 0.75, sdratio = 1, stimulus = 1,
    dist = "normal", Ro = 0.3, Rn = 0.2
  )
  expect_equal(sum(probs_old), 1, tolerance = 1e-10)

  probs_new <- bmm:::.sdt_dpsdt_category_probs(
    thresholds, shift = -0.75, sdratio = 1, stimulus = 0,
    dist = "normal", Ro = 0.3, Rn = 0.2
  )
  expect_equal(sum(probs_new), 1, tolerance = 1e-10)
})


############################################################################# !
# FORMULA CONSTRUCTION TESTS                                              ####
############################################################################# !

test_that("sdt_dp produces valid stancode with equidistant thresholds", {
  dat <- rsdt_dp(n_per_cell = 50, n_subjects = 3, dprime = 1.0, criterion = 0,
                 Ro = 0.3, Rn = 0.1, n_ratings = 4, spacing = 0.5)
  model <- sdt_dp(c("r1", "r2", "r3", "r4"), "stimulus")
  formula <- bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1, Ro ~ 1, Rn ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
  expect_true(grepl("sdt_dp_lpmf", code, fixed = TRUE))
  expect_true(grepl("inv_logit", code))
  expect_true(grepl("Ro", code))
  expect_true(grepl("Rn", code))
})

test_that("sdt_dp produces valid stancode with log_distance thresholds", {
  dat <- rsdt_dp(n_per_cell = 50, n_subjects = 3, dprime = 1.0, criterion = 0,
                 Ro = 0.3, Rn = 0, n_ratings = 4, spacing = 0.5)
  model <- sdt_dp(c("r1", "r2", "r3", "r4"), "stimulus",
                  threshold_type = "log_distance")
  formula <- bmf(dprime ~ 1, criterion ~ 1, delta1 ~ 1, delta3 ~ 1,
                 Ro ~ 1, Rn ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
  expect_true(grepl("inv_logit", code))
})

test_that("sdt_dp produces valid stancode with K=6", {
  dat <- rsdt_dp(n_per_cell = 50, n_subjects = 2, dprime = 1.0, criterion = 0,
                 Ro = 0.3, Rn = 0, n_ratings = 6, spacing = 0.3)
  model <- sdt_dp(c("r1", "r2", "r3", "r4", "r5", "r6"), "stimulus")
  formula <- bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1, Ro ~ 1, Rn ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
})

test_that("sdt_dp produces valid stancode with predictors", {
  dat <- rsdt_dp(n_per_cell = 50, n_subjects = 3, dprime = 1.0, criterion = 0,
                 Ro = 0.3, Rn = 0, n_ratings = 4, spacing = 0.5)
  dat$condition <- rep(c("A", "B"), length.out = nrow(dat))
  model <- sdt_dp(c("r1", "r2", "r3", "r4"), "stimulus")
  formula <- bmf(dprime ~ condition, criterion ~ 1, spacing ~ 1,
                 Ro ~ condition, Rn ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
})

test_that("sdt_dp + UV-SDT (sdratio ~ 1) produces valid stancode", {
  dat <- rsdt_dp(n_per_cell = 50, n_subjects = 3, dprime = 1.0, criterion = 0,
                 Ro = 0.3, Rn = 0, sdratio = 1.2, n_ratings = 4, spacing = 0.5)
  model <- sdt_dp(c("r1", "r2", "r3", "r4"), "stimulus")
  formula <- bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1,
                 sdratio ~ 1, Ro ~ 1, Rn ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
  expect_true(grepl("sdratio", code))
  expect_true(grepl("inv_logit", code))
})

test_that("sdt_dp with all 4 distributions produces valid stancode", {
  dat <- rsdt_dp(n_per_cell = 50, n_subjects = 2, dprime = 1.0, criterion = 0,
                 Ro = 0.3, Rn = 0, n_ratings = 4, spacing = 0.5)
  for (d in c("normal", "logistic", "gumbel_min", "gumbel_max")) {
    model <- sdt_dp(c("r1", "r2", "r3", "r4"), "stimulus", dist = d)
    formula <- bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1, Ro ~ 1, Rn ~ 1)
    code <- stancode(formula, data = dat, model = model)
    expect_true(nchar(code) > 0)
  }
})
