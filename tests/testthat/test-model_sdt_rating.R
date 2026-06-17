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

test_that("sdt_rating model rejects n_ratings <= 2", {
  expect_error(sdt_rating(c("r1", "r2"), "stimulus"), "n_ratings > 2")
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

test_that("sdt_rating model does not estimate a base mu", {
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  expect_false("mu" %in% names(model$parameters))
})

test_that("sdt_rating stores threshold_type correctly", {
  model_pa <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  expect_equal(model_pa$other_vars$threshold_type, "parsimonious")

  model_ld <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus",
                         threshold_type = "log_distance")
  expect_equal(model_ld$other_vars$threshold_type, "log_distance")
})

test_that("sdt_rating model stores all distribution options", {
  for (di in c("normal", "logistic", "gumbel_min", "gumbel_max")) {
    model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus", dist = di)
    expect_equal(model$other_vars$dist, di)
    expect_equal(model$other_vars$dist_int, .sdt_dist_id(di))
  }
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

test_that("sdt_rating model has softmax threshold params", {
  model <- sdt_rating(paste0("r", 1:6), "stimulus", threshold_type = "softmax")
  expect_true("spacing" %in% names(model$parameters))
  # K = 6 -> K - 3 = 3 allocation deltas
  expect_true(all(c("delta1", "delta2", "delta3") %in% names(model$parameters)))
})

test_that("sdt_rating always has sdratio as parameter (fixed to 0 by default)", {
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  expect_true("sdratio" %in% names(model$parameters))
  expect_true("sdratio" %in% names(model$links))
  expect_true("sdratio" %in% names(model$default_priors))
  expect_equal(model$fixed_parameters$sdratio, 0)
})

test_that("sdt_rating model accepts custom links", {
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus",
                      links = list(dprime = "log"))
  expect_equal(model$links$dprime, "log")
  expect_equal(model$links$criterion, "identity")
})

test_that("sdt_rating supplies init_ranges for every estimated parameter", {
  # create_initfun() looks up init_ranges per parameter; a missing entry yields
  # NA inits, and the flexible threshold types reject brms' default random init.
  for (tt in c("parsimonious", "equidistant", "log_distance",
               "log_ratio", "softmax")) {
    model <- sdt_rating(paste0("r", 1:6), "stimulus", threshold_type = tt)
    est_pars <- setdiff(names(model$parameters), "sdratio")
    expect_true(all(est_pars %in% names(model$init_ranges)),
                info = paste("threshold_type:", tt))
    expect_true("sdratio" %in% names(model$init_ranges))
    for (r in model$init_ranges) {
      expect_length(r, 2)
      expect_lte(r[1], r[2])
    }
  }
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
  # multinomial encoding: one row per observation with a count matrix Y
  expect_equal(nrow(result), nrow(valid_data))
  expect_true(all(c("Y", "nTrials", "stimulus") %in% colnames(result)))
  expect_equal(ncol(result$Y), 4)
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

test_that("rsdt_rating works with softmax thresholds", {
  dat <- rsdt_rating(n_per_cell = 50, n_subjects = 3, dprime = 1.5, criterion = 0,
                     n_ratings = 6, spacing = 0.3, deltas = c(0, 0, 0),
                     threshold_type = "softmax")
  expect_equal(nrow(dat), 6)
  expect_true(all(rowSums(dat[, paste0("r", 1:6)]) == 50))
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

test_that("rsdt_rating rejects an unknown distribution", {
  expect_error(
    rsdt_rating(n_per_cell = 50, n_subjects = 2, dprime = 1.5, criterion = 0,
                n_ratings = 4, spacing = 0.5, dist = "cauchy"),
    "dist must be one of"
  )
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

test_that("sdt_rating produces valid stancode with parsimonious thresholds", {
  dat <- rsdt_rating(n_per_cell = 50, n_subjects = 3, dprime = 1.5, criterion = 0,
                     n_ratings = 4, spacing = 0.5)
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  formula <- bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
  expect_true(grepl("sdt_rating_logmu", code, fixed = TRUE))
  expect_true(grepl("sdt_thresholds_parsimonious_rating", code, fixed = TRUE))
})

test_that("sdt_rating stancode loads the shared distribution dispatch", {
  # the noise-distribution dispatch lives in sdt_dist_funs.stan; the rating
  # likelihood relies on the log-scale dispatchers from that shared chunk.
  dat <- rsdt_rating(n_per_cell = 50, n_subjects = 3, dprime = 1.5, criterion = 0,
                     n_ratings = 4, spacing = 0.5)
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  formula <- bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(grepl("sdt_log_cumprob", code, fixed = TRUE))
  expect_true(grepl("sdt_log_one_minus_cumprob", code, fixed = TRUE))
  expect_true(grepl("log_diff_exp", code, fixed = TRUE))
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
  expect_true(grepl("sdt_rating_logmu", code, fixed = TRUE))
  expect_true(grepl("sdt_thresholds_log_distance_rating", code, fixed = TRUE))
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
  expect_true(grepl("sdt_thresholds_log_ratio_rating", code, fixed = TRUE))
  expect_true(grepl("sdt_rating_logmu", code, fixed = TRUE))
})

test_that("sdt_rating produces valid stancode with softmax thresholds", {
  dat <- rsdt_rating(n_per_cell = 50, n_subjects = 3, dprime = 1.5, criterion = 0,
                     n_ratings = 6, spacing = 0.3, deltas = c(0, 0, 0),
                     threshold_type = "softmax")
  model <- sdt_rating(paste0("r", 1:6), "stimulus", threshold_type = "softmax")
  formula <- bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1,
                 delta1 ~ 1, delta2 ~ 1, delta3 ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
  expect_true(grepl("sdt_thresholds_softmax_rating", code, fixed = TRUE))
})

test_that("sdt_rating handles K=6 parsimonious", {
  dat <- rsdt_rating(n_per_cell = 50, n_subjects = 2, dprime = 1.5, criterion = 0,
                     n_ratings = 6, spacing = 0.3)
  model <- sdt_rating(paste0("r", 1:6), "stimulus")
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
  expect_true(grepl("sdt_rating_logmu", code, fixed = TRUE))
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


############################################################################# !
# PIPELINE INTEGRATION TESTS (mock backend)                              ####
############################################################################# !

test_that("sdt_rating integrates with the bmm pipeline via mock backend", {
  dat <- rsdt_rating(n_per_cell = 80, n_subjects = 6, dprime = 1.2, criterion = 0,
                     n_ratings = 4, spacing = 0.4)
  for (di in c("normal", "logistic", "gumbel_min", "gumbel_max")) {
    model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus", dist = di)
    expect_silent(
      bmm(bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1),
          dat, model, backend = "mock", mock_fit = 1, rename = FALSE)
    )
  }
})

test_that("sdt_rating UV-SDT integrates with the bmm pipeline via mock backend", {
  dat <- rsdt_rating(n_per_cell = 80, n_subjects = 6, dprime = 1.2, criterion = 0,
                     sdratio = 1.3, n_ratings = 4, spacing = 0.4)
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  expect_silent(
    bmm(bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1, sdratio ~ 1),
        dat, model, backend = "mock", mock_fit = 1, rename = FALSE)
  )
})

test_that("sdt_rating log_distance integrates with the bmm pipeline via mock backend", {
  # log_distance has variable-arity delta parameters, so this exercises the
  # generated multinomial formula and Stan function for that case.
  dat <- rsdt_rating(n_per_cell = 60, n_subjects = 5, dprime = 1.5, criterion = 0,
                     n_ratings = 4, deltas = c(0.5, 0.5),
                     threshold_type = "log_distance")
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus",
                      threshold_type = "log_distance")
  expect_silent(
    bmm(bmf(dprime ~ 1, criterion ~ 1, delta1 ~ 1, delta3 ~ 1),
        dat, model, backend = "mock", mock_fit = 1, rename = FALSE)
  )
})

test_that("sdt_rating posterior_predict draws joint multinomial counts (real fit)", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("cmdstanr")
  skip_if(is.null(tryCatch(cmdstanr::cmdstan_version(), error = function(e) NULL)))

  dat <- rsdt_rating(n_per_cell = 120, n_subjects = 6, dprime = 1.4, criterion = 0,
                     n_ratings = 4, spacing = 0.5)
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  fit <- bmm(bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1), dat, model,
             backend = "cmdstanr", chains = 1, iter = 300, warmup = 150,
             refresh = 0, silent = 2)

  pp <- brms::posterior_predict(fit)
  expect_length(dim(pp), 3)
  row_sums <- apply(pp, c(1, 2), sum)
  expect_true(all(row_sums == 120))
})

test_that("sdt_rating log_distance fits with finite likelihood at K=6 (real fit)", {
  # Regression: at K=6 the log_distance threshold builder has mid=3, so its
  # lower thresholds were filled by a descending Stan for-loop. Stan for-loops
  # only count up, so the loop was skipped and thresholds[1:2] stayed NaN,
  # failing every chain. Earlier tests only used K=4 (mid=2), where the loop ran.
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("cmdstanr")
  skip_if(is.null(tryCatch(cmdstanr::cmdstan_version(), error = function(e) NULL)))

  dat <- rsdt_rating(n_per_cell = 150, n_subjects = 6, dprime = 1.4, criterion = 0,
                     n_ratings = 6, deltas = c(-0.4, -0.4, -0.4, -0.4),
                     threshold_type = "log_distance")
  model <- sdt_rating(paste0("r", 1:6), "stimulus", threshold_type = "log_distance")
  fit <- bmm(bmf(dprime ~ 1, criterion ~ 1, delta1 ~ 1, delta2 ~ 1,
                 delta4 ~ 1, delta5 ~ 1),
             dat, model, backend = "cmdstanr", chains = 1, iter = 300, warmup = 150,
             refresh = 0, silent = 2)

  expect_true(all(is.finite(brms::log_lik(fit))))
})
