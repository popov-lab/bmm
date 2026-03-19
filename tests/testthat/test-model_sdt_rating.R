# Test Rating SDT model specification and integration

############################################################################# !
# MODEL CONSTRUCTOR TESTS                                                ####
############################################################################# !

test_that("sdt_rating model can be created with vector response", {
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  expect_s3_class(model, "bmmodel")
  expect_s3_class(model, "sdt")
  expect_s3_class(model, "sdt_rating")
  expect_equal(model$other_vars$n_ratings, 4L)
})

test_that("sdt_rating model infers n_ratings from response length", {
  model <- sdt_rating(c("r1", "r2", "r3", "r4", "r5", "r6"), "stimulus")
  expect_equal(model$other_vars$n_ratings, 6L)
})

test_that("sdt_rating model rejects mismatched n_ratings", {
  expect_error(
    sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus", n_ratings = 6),
    "must match"
  )
})

test_that("sdt_rating model has parsimonious threshold params by default", {
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  expect_true("spacing" %in% names(model$parameters))
  expect_true("dprime" %in% names(model$parameters))
  expect_true("criterion" %in% names(model$parameters))
  expect_false("mu" %in% names(model$parameters))
})

test_that("sdt_rating model has log_distance threshold params", {
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus",
                      threshold_type = "log_distance")
  expect_true("delta1" %in% names(model$parameters))
  expect_true("delta3" %in% names(model$parameters))
  expect_false("delta2" %in% names(model$parameters))
  expect_false("spacing" %in% names(model$parameters))
})

test_that("sdt_rating model has void_mu = TRUE", {
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  expect_true(model$void_mu)
})

test_that("sdt_rating stores threshold_type correctly", {
  model_pa <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  expect_equal(model_pa$other_vars$threshold_type, "parsimonious")

  model_ld <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus",
                         threshold_type = "log_distance")
  expect_equal(model_ld$other_vars$threshold_type, "log_distance")
})

test_that("sdt_rating model has parsimonious threshold params", {
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus",
                      threshold_type = "parsimonious")
  expect_true("spacing" %in% names(model$parameters))
  expect_true("criterion" %in% names(model$parameters))
  expect_true("dprime" %in% names(model$parameters))
  expect_false("delta1" %in% names(model$parameters))
  expect_false("delta3" %in% names(model$parameters))
})

test_that("sdt_rating model stores log_ratio threshold_type", {
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus",
                      threshold_type = "log_ratio")
  expect_equal(model$other_vars$threshold_type, "log_ratio")
})

test_that("sdt_rating model has log_ratio threshold params", {
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus",
                      threshold_type = "log_ratio")
  expect_true("delta1" %in% names(model$parameters))
  expect_true("delta3" %in% names(model$parameters))
  expect_true("criterion" %in% names(model$parameters))
  expect_false("spacing" %in% names(model$parameters))
  expect_false("delta2" %in% names(model$parameters))
})

test_that("sdt_rating always has sdratio as parameter (fixed to 0 by default)", {
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  expect_true("sdratio" %in% names(model$parameters))
  expect_true("sdratio" %in% names(model$links))
  expect_true("sdratio" %in% names(model$default_priors))
  expect_equal(model$fixed_parameters$sdratio, 0)
})

test_that("sdt_rating model has NULL init_ranges", {
  model <- sdt_rating(paste0("r", 1:6), "stimulus")
  expect_null(model$init_ranges)
})


############################################################################# !
# CHECK_DATA TESTS                                                       ####
############################################################################# !

test_that("sdt_rating check_data validates response columns", {
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  formula <- bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1)

  valid_data <- data.frame(
    r1 = c(10, 5), r2 = c(20, 10), r3 = c(15, 30), r4 = c(5, 55),
    stimulus = c(0L, 1L)
  )
  result <- check_data(model, valid_data, formula)
  expect_true("Y" %in% colnames(result))
  expect_true("nTrials" %in% colnames(result))
  expect_equal(result$nTrials, c(50, 100))
})

test_that("sdt_rating check_data rejects missing response columns", {
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  formula <- bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1)

  invalid_data <- data.frame(
    r1 = c(10, 5), r2 = c(20, 10),
    stimulus = c(0L, 1L)
  )
  expect_error(check_data(model, invalid_data, formula), "missing in the data")
})

test_that("sdt_rating check_data rejects negative counts", {
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  formula <- bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1)

  invalid_data <- data.frame(
    r1 = c(-1, 5), r2 = c(20, 10), r3 = c(15, 30), r4 = c(5, 55),
    stimulus = c(0L, 1L)
  )
  expect_error(check_data(model, invalid_data, formula), "non-negative")
})

test_that("sdt_rating check_data validates stimulus coding", {
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  formula <- bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1)

  invalid_data <- data.frame(
    r1 = c(10, 5), r2 = c(20, 10), r3 = c(15, 30), r4 = c(5, 55),
    stimulus = c(0L, 2L)
  )
  expect_error(check_data(model, invalid_data, formula), "must be coded as 0")
})


############################################################################# !
# DISTRIBUTION FUNCTION TESTS                                            ####
############################################################################# !

test_that("rsdt_rating generates correct structure", {
  dat <- rsdt_rating(n_per_cell = 50, n_subjects = 5, dprime = 1.5, criterion = 0,
                     n_ratings = 4, spacing = 0.5)
  expect_true(is.data.frame(dat))
  expect_equal(nrow(dat), 10)
  expect_true(all(c("id", "stimulus", "r1", "r2", "r3", "r4", "nTrials") %in%
                    colnames(dat)))
  expect_equal(unique(dat$stimulus), c(0L, 1L))
  expect_true(all(rowSums(dat[, paste0("r", 1:4)]) == 50))
})

test_that("rsdt_rating works with log_distance thresholds", {
  dat <- rsdt_rating(n_per_cell = 50, n_subjects = 3, dprime = 1.5, criterion = 0,
                     n_ratings = 4, deltas = c(0.5, 0.5),
                     threshold_type = "log_distance")
  expect_equal(nrow(dat), 6)
  expect_true(all(rowSums(dat[, paste0("r", 1:4)]) == 50))
})

test_that("rsdt_rating works with parsimonious thresholds", {
  dat <- rsdt_rating(n_per_cell = 50, n_subjects = 3, dprime = 1.5, criterion = 0,
                     n_ratings = 4, spacing = 0.5,
                     threshold_type = "parsimonious")
  expect_equal(nrow(dat), 6)
  expect_true(all(rowSums(dat[, paste0("r", 1:4)]) == 50))
})

test_that("rsdt_rating requires n_ratings >= 3", {
  expect_error(
    rsdt_rating(n_per_cell = 50, n_subjects = 1, dprime = 1, criterion = 0,
                n_ratings = 2, spacing = 0.5),
    "n_ratings must be >= 3"
  )
})

test_that("rsdt_rating works with different distributions", {
  for (di in c("normal", "logistic", "gumbel_min", "gumbel_max")) {
    dat <- rsdt_rating(n_per_cell = 50, n_subjects = 2, dprime = 1.5, criterion = 0,
                       n_ratings = 4, spacing = 0.5, dist = di)
    expect_equal(nrow(dat), 4, info = paste("dist:", di))
    expect_true(all(rowSums(dat[, paste0("r", 1:4)]) == 50),
                info = paste("dist:", di))
  }
})

test_that("rsdt_rating with K=6 generates 6 response columns", {
  dat <- rsdt_rating(n_per_cell = 100, n_subjects = 2, dprime = 1.5, criterion = 0,
                     n_ratings = 6, spacing = 0.3)
  expect_true(all(paste0("r", 1:6) %in% colnames(dat)))
  expect_true(all(rowSums(dat[, paste0("r", 1:6)]) == 100))
})

test_that("rsdt_rating works with log_ratio thresholds", {
  dat <- rsdt_rating(n_per_cell = 50, n_subjects = 3, dprime = 1.5, criterion = 0,
                     n_ratings = 4, deltas = c(0, 0),
                     threshold_type = "log_ratio")
  expect_equal(nrow(dat), 6)
  expect_true(all(rowSums(dat[, paste0("r", 1:4)]) == 50))
})

test_that("rsdt_rating rejects vector dprime", {
  expect_error(
    rsdt_rating(n_per_cell = 50, n_subjects = 2, dprime = c(1, 2),
                criterion = 0, n_ratings = 4, spacing = 0.5),
    "dprime must be a single value"
  )
})

test_that("rsdt_rating rejects vector criterion", {
  expect_error(
    rsdt_rating(n_per_cell = 50, n_subjects = 2, dprime = 1.5,
                criterion = c(0, 0.5), n_ratings = 4, spacing = 0.5),
    "criterion must be a single value"
  )
})

test_that("rsdt_rating rejects vector sdratio", {
  expect_error(
    rsdt_rating(n_per_cell = 50, n_subjects = 2, dprime = 1.5,
                criterion = 0, sdratio = c(1.0, 1.3),
                n_ratings = 4, spacing = 0.5),
    "sdratio must be a single value"
  )
})

test_that("dsdt_rating returns valid density", {
  d <- dsdt_rating(counts = c(10, 20, 30, 40), stimulus = 1,
                   dprime = 1.5, thresholds = c(-0.5, 0, 0.5))
  expect_true(d > 0)
  expect_true(d <= 1)
})

test_that("dsdt_rating log matches log of density", {
  d <- dsdt_rating(counts = c(10, 20, 30, 40), stimulus = 1,
                   dprime = 1.5, thresholds = c(-0.5, 0, 0.5))
  ld <- dsdt_rating(counts = c(10, 20, 30, 40), stimulus = 1,
                    dprime = 1.5, thresholds = c(-0.5, 0, 0.5), log = TRUE)
  expect_equal(log(d), ld, tolerance = 1e-10)
})

test_that("dsdt_rating validates inputs", {
  expect_error(
    dsdt_rating(counts = c(-1, 20, 30, 40), stimulus = 1,
                dprime = 1.5, thresholds = c(-0.5, 0, 0.5)),
    "non-negative"
  )
  expect_error(
    dsdt_rating(counts = c(10, 20, 30), stimulus = 1,
                dprime = 1.5, thresholds = c(-0.5, 0, 0.5)),
    "K - 1"
  )
})

test_that("dsdt_rating works for all distributions", {
  for (di in c("normal", "logistic", "gumbel_min", "gumbel_max")) {
    d <- dsdt_rating(counts = c(10, 20, 30, 40), stimulus = 1,
                     dprime = 1.5, thresholds = c(-0.5, 0, 0.5), dist = di)
    expect_true(d > 0, info = paste("dist:", di))
  }
})

test_that("dsdt_rating with sdratio != 1 changes density", {
  d_ev <- dsdt_rating(counts = c(10, 20, 30, 40), stimulus = 1,
                      dprime = 1.5, thresholds = c(-0.5, 0, 0.5), sdratio = 1)
  d_uv <- dsdt_rating(counts = c(10, 20, 30, 40), stimulus = 1,
                      dprime = 1.5, thresholds = c(-0.5, 0, 0.5), sdratio = 1.3)
  expect_false(d_ev == d_uv)
})


############################################################################# !
# FORMULA CONSTRUCTION TESTS                                              ####
############################################################################# !

test_that("sdt_rating produces valid stancode with equidistant thresholds", {
  dat <- rsdt_rating(n_per_cell = 50, n_subjects = 3, dprime = 1.5, criterion = 0,
                     n_ratings = 4, spacing = 0.5)
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  formula <- bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
  expect_true(grepl("multinomial", code, ignore.case = TRUE))
})

test_that("sdt_rating produces valid stancode with log_distance thresholds", {
  dat <- rsdt_rating(n_per_cell = 50, n_subjects = 3, dprime = 1.5, criterion = 0,
                     n_ratings = 4, deltas = c(0.5, 0.5),
                     threshold_type = "log_distance")
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus",
                      threshold_type = "log_distance")
  formula <- bmf(dprime ~ 1, criterion ~ 1, delta1 ~ 1, delta3 ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
  expect_true(grepl("multinomial", code, ignore.case = TRUE))
})

test_that("sdt_rating produces valid stancode with parsimonious thresholds", {
  dat <- rsdt_rating(n_per_cell = 50, n_subjects = 3, dprime = 1.5, criterion = 0,
                     n_ratings = 4, spacing = 0.5,
                     threshold_type = "parsimonious")
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus",
                      threshold_type = "parsimonious")
  formula <- bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
  expect_true(grepl("multinomial", code, ignore.case = TRUE))
  expect_true(grepl("1.0986122887", code))
})

test_that("sdt_rating parsimonious K=6 produces valid stancode", {
  dat <- rsdt_rating(n_per_cell = 50, n_subjects = 2, dprime = 1.5, criterion = 0,
                     n_ratings = 6, spacing = 0.5,
                     threshold_type = "parsimonious")
  model <- sdt_rating(c("r1", "r2", "r3", "r4", "r5", "r6"), "stimulus",
                      threshold_type = "parsimonious")
  formula <- bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
})

test_that("sdt_rating produces valid stancode with log_ratio thresholds", {
  dat <- rsdt_rating(n_per_cell = 50, n_subjects = 3, dprime = 1.5, criterion = 0,
                     n_ratings = 4, deltas = c(0, 0),
                     threshold_type = "log_ratio")
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus",
                      threshold_type = "log_ratio")
  formula <- bmf(dprime ~ 1, criterion ~ 1, delta1 ~ 1, delta3 ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
  expect_true(grepl("multinomial", code, ignore.case = TRUE))
  expect_true(grepl("exp\\(nlp_delta", code))
})

test_that("sdt_rating log_ratio K=6 produces valid stancode", {
  dat <- rsdt_rating(n_per_cell = 50, n_subjects = 3, dprime = 1.5, criterion = 0,
                     n_ratings = 6, deltas = c(0, 0, 0, 0),
                     threshold_type = "log_ratio")
  model <- sdt_rating(paste0("r", 1:6), "stimulus", threshold_type = "log_ratio")
  formula <- bmf(dprime ~ 1, criterion ~ 1, delta1 ~ 1, delta2 ~ 1,
                 delta4 ~ 1, delta5 ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
})

test_that("sdt_rating handles K=6 equidistant", {
  dat <- rsdt_rating(n_per_cell = 50, n_subjects = 2, dprime = 1.5, criterion = 0,
                     n_ratings = 6, spacing = 0.3)
  model <- sdt_rating(c("r1", "r2", "r3", "r4", "r5", "r6"), "stimulus")
  formula <- bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
})

test_that("sdt_rating handles predictors on dprime", {
  dat <- rsdt_rating(n_per_cell = 50, n_subjects = 3, dprime = 1.5, criterion = 0,
                     n_ratings = 4, spacing = 0.5)
  dat$condition <- rep(c("A", "B"), length.out = nrow(dat))
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  formula <- bmf(dprime ~ condition, criterion ~ 1, spacing ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
})


############################################################################# !
# UV-SDT TESTS (sdratio overridable fixed parameter)                     ####
############################################################################# !

test_that("sdt_rating with sdratio ~ 1 produces valid stancode", {
  dat <- rsdt_rating(n_per_cell = 50, n_subjects = 3, dprime = 1.5, criterion = 0,
                     sdratio = 1.3, n_ratings = 4, spacing = 0.5)
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  formula <- bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1, sdratio ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
  expect_true(grepl("multinomial", code, ignore.case = TRUE))
  expect_true(grepl("sdratio", code))
})

test_that("sdt_rating UV-SDT with log_distance thresholds produces valid stancode", {
  dat <- rsdt_rating(n_per_cell = 50, n_subjects = 3, dprime = 1.5, criterion = 0,
                     sdratio = 1.3, n_ratings = 4, deltas = c(0.5, 0.5),
                     threshold_type = "log_distance")
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus",
                      threshold_type = "log_distance")
  formula <- bmf(dprime ~ 1, criterion ~ 1, delta1 ~ 1, delta3 ~ 1,
                 sdratio ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
  expect_true(grepl("sdratio", code))
})

test_that("sdt_rating UV-SDT with K=6 produces valid stancode", {
  dat <- rsdt_rating(n_per_cell = 50, n_subjects = 2, dprime = 1.5, criterion = 0,
                     sdratio = 1.3, n_ratings = 6, spacing = 0.3)
  model <- sdt_rating(c("r1", "r2", "r3", "r4", "r5", "r6"), "stimulus")
  formula <- bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1, sdratio ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
})


############################################################################# !
# LOG-SCALE NUMERICS TESTS                                                ####
############################################################################# !

test_that("sdt_rating with log_scale=TRUE produces valid stancode (normal)", {
  dat <- rsdt_rating(n_per_cell = 100, n_subjects = 3, dprime = 1.5,
                     criterion = 0, n_ratings = 6, spacing = 0.5)
  model <- sdt_rating(response = paste0("r", 1:6), stimulus = "stimulus",
                      dist = "normal", log_scale = TRUE)
  formula <- bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
  expect_true(grepl("std_normal_lcdf", code))
  expect_true(grepl("std_normal_lccdf", code))
  expect_true(grepl("log_diff_exp", code))
})

test_that("sdt_rating with log_scale=TRUE works for gumbel_min", {
  dat <- rsdt_rating(n_per_cell = 100, n_subjects = 3, dprime = 1.5,
                     criterion = 0, n_ratings = 6, spacing = 0.5,
                     dist = "gumbel_min")
  model <- sdt_rating(response = paste0("r", 1:6), stimulus = "stimulus",
                      dist = "gumbel_min", log_scale = TRUE)
  formula <- bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
  expect_true(grepl("log_diff_exp", code))
})

test_that("sdt_rating with log_scale=TRUE works for logistic", {
  dat <- rsdt_rating(n_per_cell = 100, n_subjects = 3, dprime = 1.5,
                     criterion = 0, n_ratings = 6, spacing = 0.5,
                     dist = "logistic")
  model <- sdt_rating(response = paste0("r", 1:6), stimulus = "stimulus",
                      dist = "logistic", log_scale = TRUE)
  formula <- bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
  expect_true(grepl("log1p_exp", code))
})

test_that("sdt_rating with log_scale=TRUE works for gumbel_max", {
  dat <- rsdt_rating(n_per_cell = 100, n_subjects = 3, dprime = 1.5,
                     criterion = 0, n_ratings = 6, spacing = 0.5,
                     dist = "gumbel_max")
  model <- sdt_rating(response = paste0("r", 1:6), stimulus = "stimulus",
                      dist = "gumbel_max", log_scale = TRUE)
  formula <- bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
  expect_true(grepl("log_diff_exp", code))
})

test_that("sdt_rating with log_scale=FALSE is unchanged (default)", {
  dat <- rsdt_rating(n_per_cell = 100, n_subjects = 3, dprime = 1.5,
                     criterion = 0, n_ratings = 6, spacing = 0.5)
  model <- sdt_rating(response = paste0("r", 1:6), stimulus = "stimulus")
  formula <- bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
  expect_true(grepl("Phi\\(", code))
  expect_false(grepl("log_Phi", code))
})

test_that("sdt_rating UV-SDT with log_scale=TRUE produces valid stancode", {
  dat <- rsdt_rating(n_per_cell = 100, n_subjects = 3, dprime = 1.5,
                     criterion = 0, sdratio = 1.2, n_ratings = 6,
                     spacing = 0.5)
  model <- sdt_rating(response = paste0("r", 1:6), stimulus = "stimulus",
                      log_scale = TRUE)
  formula <- bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1, sdratio ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
  expect_true(grepl("std_normal_lcdf", code))
  expect_true(grepl("sdratio", code))
})
