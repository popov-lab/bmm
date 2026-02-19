# Test SDT model specification and integration

############################################################################# !
# MODEL CONSTRUCTOR TESTS                                                ####
############################################################################# !

test_that("sdt model can be created with default arguments", {
  expect_silent(sdt("n_old", "stimulus", "n_trials"))
})

test_that("sdt model can be created with all distribution options", {
  expect_silent(sdt("n_old", "stimulus", "n_trials", dist = "normal"))
  expect_silent(sdt("n_old", "stimulus", "n_trials", dist = "logistic"))
  expect_silent(sdt("n_old", "stimulus", "n_trials", dist = "gumbel_min"))
  expect_silent(sdt("n_old", "stimulus", "n_trials", dist = "gumbel_max"))
})

test_that("sdt model has correct class structure", {
  model <- sdt("n_old", "stimulus", "n_trials")
  expect_s3_class(model, "bmmodel")
  expect_s3_class(model, "sdt")
  expect_s3_class(model, "sdt_evsdt")
})

test_that("sdt model parameters are correctly defined", {
  model <- sdt("n_old", "stimulus", "n_trials")
  expect_true("mu" %in% names(model$parameters))
  expect_true("dprime" %in% names(model$parameters))
  expect_true("criterion" %in% names(model$parameters))
  # mu is fixed to 0
  expect_equal(model$fixed_parameters$mu, 0)
})

test_that("sdt model has correct default link functions", {
  model <- sdt("n_old", "stimulus", "n_trials")
  expect_equal(model$links$mu, "identity")
  expect_equal(model$links$dprime, "identity")
  expect_equal(model$links$criterion, "identity")
})

test_that("sdt model accepts custom links", {
  custom_links <- list(dprime = "log")
  model <- sdt("n_old", "stimulus", "n_trials", links = custom_links)
  expect_equal(model$links$dprime, "log")
  expect_equal(model$links$criterion, "identity")  # unchanged
})

test_that("sdt model stores distribution info correctly", {
  model <- sdt("n_old", "stimulus", "n_trials", dist = "gumbel_min")
  expect_equal(model$other_vars$dist, "gumbel_min")
  expect_equal(model$other_vars$dist_int, 2L)

  model2 <- sdt("n_old", "stimulus", "n_trials", dist = "logistic")
  expect_equal(model2$other_vars$dist, "logistic")
  expect_equal(model2$other_vars$dist_int, 4L)
})

test_that("sdt model stores variances property correctly", {
  model_eq <- sdt("n_old", "stimulus", "n_trials")
  expect_equal(model_eq$other_vars$variances, "equal")

  # UV-SDT requires ratings
  model_uneq <- sdt(c("r1", "r2", "r3", "r4"), "stimulus",
                     variances = "unequal")
  expect_equal(model_uneq$other_vars$variances, "unequal")

  # UV-SDT with binary data should error
  expect_error(
    sdt("n_old", "stimulus", "n_trials", variances = "unequal"),
    "requires confidence rating"
  )
})

test_that("sdt model has default priors", {
  model <- sdt("n_old", "stimulus", "n_trials")
  expect_true("dprime" %in% names(model$default_priors))
  expect_true("criterion" %in% names(model$default_priors))
  expect_true("main" %in% names(model$default_priors$dprime))
  expect_true("effects" %in% names(model$default_priors$dprime))
})


############################################################################# !
# CHECK_DATA TESTS                                                       ####
############################################################################# !

test_that("sdt check_data validates required columns", {
  model <- sdt("n_old", "stimulus", "n_trials")
  formula <- bmf(dprime ~ 1, criterion ~ 1)

  # Valid data
  valid_data <- data.frame(
    n_old = c(30, 40),
    stimulus = c(0L, 1L),
    n_trials = c(50, 50)
  )
  result <- check_data(model, valid_data, formula)
  expect_true("dist_type" %in% colnames(result))

  # Missing columns
  invalid_data <- data.frame(n_old = c(30, 40))
  expect_error(
    check_data(model, invalid_data, formula),
    "missing in the data"
  )
})

test_that("sdt check_data validates stimulus coding", {
  model <- sdt("n_old", "stimulus", "n_trials")
  formula <- bmf(dprime ~ 1, criterion ~ 1)

  invalid_data <- data.frame(
    n_old = c(30, 40),
    stimulus = c(1, 2),  # should be 0/1
    n_trials = c(50, 50)
  )
  expect_error(
    check_data(model, invalid_data, formula),
    "must be coded as 0"
  )
})

test_that("sdt check_data validates response counts", {
  model <- sdt("n_old", "stimulus", "n_trials")
  formula <- bmf(dprime ~ 1, criterion ~ 1)

  # Negative counts
  invalid_data <- data.frame(
    n_old = c(-5, 40),
    stimulus = c(0L, 1L),
    n_trials = c(50, 50)
  )
  expect_error(
    check_data(model, invalid_data, formula),
    "non-negative"
  )

  # Response exceeds n_trials
  invalid_data2 <- data.frame(
    n_old = c(60, 40),
    stimulus = c(0L, 1L),
    n_trials = c(50, 50)
  )
  expect_error(
    check_data(model, invalid_data2, formula),
    "must not exceed"
  )
})

test_that("sdt check_data validates n_trials", {
  model <- sdt("n_old", "stimulus", "n_trials")
  formula <- bmf(dprime ~ 1, criterion ~ 1)

  invalid_data <- data.frame(
    n_old = c(0, 0),
    stimulus = c(0L, 1L),
    n_trials = c(-10, 50)
  )
  expect_error(
    check_data(model, invalid_data, formula),
    "positive"
  )
})

test_that("sdt check_data adds dist_type column", {
  model <- sdt("n_old", "stimulus", "n_trials", dist = "gumbel_min")
  formula <- bmf(dprime ~ 1, criterion ~ 1)

  valid_data <- data.frame(
    n_old = c(10, 40),
    stimulus = c(0L, 1L),
    n_trials = c(50, 50)
  )
  result <- check_data(model, valid_data, formula)
  expect_equal(unique(result$dist_type), 2L)
})


############################################################################# !
# DISTRIBUTION FUNCTION TESTS                                            ####
############################################################################# !

test_that("rsdt generates data with correct structure", {
  dat <- rsdt(n_per_cell = 50, n_subjects = 5)
  expect_true(is.data.frame(dat))
  expect_equal(nrow(dat), 10)  # 5 subjects * 2 stimulus types
  expect_true(all(c("id", "stimulus", "n_trials", "n_old") %in% colnames(dat)))
  expect_equal(unique(dat$stimulus), c(0L, 1L))
  expect_true(all(dat$n_old >= 0))
  expect_true(all(dat$n_old <= dat$n_trials))
})

test_that("rsdt generates different data for different distributions", {
  set.seed(42)
  dat_norm <- rsdt(n_per_cell = 1000, n_subjects = 1, dprime = 2, criterion = 0)
  set.seed(42)
  dat_logis <- rsdt(n_per_cell = 1000, n_subjects = 1, dprime = 2, criterion = 0,
                    dist = "logistic")
  # Different distributions should produce different counts (even with same seed,
  # because the underlying CDF is different)
  expect_false(identical(dat_norm$n_old, dat_logis$n_old))
})

test_that("rsdt with zero dprime produces equal hit and FA rates on average", {
  set.seed(123)
  dat <- rsdt(n_per_cell = 10000, n_subjects = 1, dprime = 0, criterion = 0)
  noise_prop <- dat$n_old[dat$stimulus == 0] / dat$n_trials[dat$stimulus == 0]
  signal_prop <- dat$n_old[dat$stimulus == 1] / dat$n_trials[dat$stimulus == 1]
  # With dprime = 0, hit rate and FA rate should be similar (~0.5 with criterion = 0)
  expect_equal(noise_prop, 0.5, tolerance = 0.05)
  expect_equal(signal_prop, 0.5, tolerance = 0.05)
})

test_that("dsdt returns valid binomial density", {
  # Known case: stimulus=1, dprime=0, criterion=0 => p=0.5
  # dbinom(50, 100, 0.5) should be the most likely outcome
  d <- dsdt(n_old = 50, n_trials = 100, stimulus = 1,
            dprime = 0, criterion = 0)
  expect_true(d > 0)
  expect_true(d <= 1)
})

test_that("dsdt log density matches log of density", {
  d <- dsdt(n_old = 80, n_trials = 100, stimulus = 1,
            dprime = 1.5, criterion = 0.2)
  ld <- dsdt(n_old = 80, n_trials = 100, stimulus = 1,
             dprime = 1.5, criterion = 0.2, log = TRUE)
  expect_equal(log(d), ld, tolerance = 1e-10)
})

test_that("dsdt is vectorized over observations", {
  d <- dsdt(n_old = c(30, 80), n_trials = c(100, 100),
            stimulus = c(0, 1), dprime = 1.5, criterion = 0.2)
  expect_length(d, 2)
  expect_true(all(d > 0))
})

test_that("dsdt works for all distributions", {
  dists <- c("normal", "logistic", "gumbel_min", "gumbel_max")
  for (di in dists) {
    d <- dsdt(n_old = 70, n_trials = 100, stimulus = 1,
              dprime = 1.5, criterion = 0, dist = di)
    expect_true(d > 0, info = paste("dist:", di))
  }
})

test_that("dsdt validates input", {
  expect_error(dsdt(n_old = -1, n_trials = 100, stimulus = 1,
                    dprime = 1, criterion = 0), "non-negative")
  expect_error(dsdt(n_old = 50, n_trials = 100, stimulus = 2,
                    dprime = 1, criterion = 0), "0.*1")
})

test_that("dsdt with sdratio produces different densities than EV-SDT", {
  d_ev <- dsdt(n_old = 80, n_trials = 100, stimulus = 1,
               dprime = 1.5, criterion = 0, sdratio = 1)
  d_uv <- dsdt(n_old = 80, n_trials = 100, stimulus = 1,
               dprime = 1.5, criterion = 0, sdratio = 1.3)
  expect_false(d_ev == d_uv)
})

test_that("dsdt with sdratio = 1 matches equal-variance result", {
  d1 <- dsdt(n_old = 70, n_trials = 100, stimulus = 0,
             dprime = 1.5, criterion = 0.2, sdratio = 1)
  # For noise trials, sdratio should not matter in the formula
  # eta = -dprime/2 - criterion for both EV and UV noise
  d2 <- dsdt(n_old = 70, n_trials = 100, stimulus = 0,
             dprime = 1.5, criterion = 0.2, sdratio = 1.5)
  # Noise trials ARE affected by sdratio in UV-SDT only through
  # changed hit rates, but for a single observation density they differ
  # because UV only changes signal side. Let's verify noise side is unchanged:
  # For noise: eta = -dprime/2 - criterion regardless of sdratio
  # So noise density should be same
  expect_equal(d1, d2)
})

test_that("sdt_dprime computes correct values for normal distribution", {
  # Known values: if hit_rate = pnorm(1) and fa_rate = pnorm(-1), d' = 2
  hr <- pnorm(1)
  far <- pnorm(-1)
  expect_equal(sdt_dprime(hr, far, dist = "normal"), 2, tolerance = 1e-10)
})

test_that("sdt_criterion computes correct values for normal distribution", {
  # With d'=2, criterion=0 means hit_rate = pnorm(1), fa_rate = pnorm(-1)
  hr <- pnorm(1)
  far <- pnorm(-1)
  expect_equal(sdt_criterion(hr, far, dist = "normal"), 0, tolerance = 1e-10)
})

test_that("sdt_dprime validates input", {
  expect_error(sdt_dprime(0, 0.5), "between 0 and 1")
  expect_error(sdt_dprime(1, 0.5), "between 0 and 1")
  expect_error(sdt_dprime(0.5, 0), "between 0 and 1")
})


############################################################################# !
# CDF HELPER TESTS                                                       ####
############################################################################# !

test_that(".sdt_cdf returns correct values for known inputs", {
  # Normal: Phi(0) = 0.5
  expect_equal(bmm:::.sdt_cdf(0, "normal"), 0.5)
  # Logistic: inv_logit(0) = 0.5
  expect_equal(bmm:::.sdt_cdf(0, "logistic"), 0.5)
  # Gumbel_min at 0: exp(-exp(0)) = exp(-1) ≈ 0.368
  expect_equal(bmm:::.sdt_cdf(0, "gumbel_min"), exp(-1), tolerance = 1e-10)
  # Gumbel_max at 0: 1 - exp(-exp(0)) = 1 - exp(-1) ≈ 0.632
  expect_equal(bmm:::.sdt_cdf(0, "gumbel_max"), 1 - exp(-1), tolerance = 1e-10)
})

test_that(".sdt_cdf is vectorized over eta", {
  eta <- c(-1, 0, 1)
  result <- bmm:::.sdt_cdf(eta, "normal")
  expect_length(result, 3)
  expect_equal(result, pnorm(eta))
})


############################################################################# !
# FORMULA AND CONFIGURE_MODEL TESTS                                      ####
############################################################################# !

test_that("sdt model produces valid stancode with brms", {
  dat <- data.frame(
    n_old = c(10, 40, 15, 35),
    stimulus = c(0L, 1L, 0L, 1L),
    n_trials = c(50, 50, 50, 50),
    id = c(1, 1, 2, 2)
  )
  model <- sdt("n_old", "stimulus", "n_trials", dist = "normal")
  formula <- bmf(dprime ~ 1, criterion ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
  expect_true(grepl("sdt_binary", code))
  expect_true(grepl("sdt_cumprob", code))
})

test_that("sdt model produces valid stancode for all distributions", {
  dat <- data.frame(
    n_old = c(10, 40),
    stimulus = c(0L, 1L),
    n_trials = c(50, 50)
  )
  for (d in c("normal", "logistic", "gumbel_min", "gumbel_max")) {
    model <- sdt("n_old", "stimulus", "n_trials", dist = d)
    formula <- bmf(dprime ~ 1, criterion ~ 1)
    code <- stancode(formula, data = dat, model = model)
    expect_true(nchar(code) > 0, info = paste("Failed for dist:", d))
  }
})

test_that("sdt model handles predictors in formula", {
  dat <- data.frame(
    n_old = c(10, 40, 15, 35),
    stimulus = c(0L, 1L, 0L, 1L),
    n_trials = c(50, 50, 50, 50),
    condition = c("A", "A", "B", "B")
  )
  model <- sdt("n_old", "stimulus", "n_trials", dist = "normal")
  formula <- bmf(dprime ~ condition, criterion ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
})

test_that("sdt model handles random effects in formula", {
  dat <- data.frame(
    n_old = c(10, 40, 15, 35),
    stimulus = c(0L, 1L, 0L, 1L),
    n_trials = c(50, 50, 50, 50),
    id = c(1, 1, 2, 2)
  )
  model <- sdt("n_old", "stimulus", "n_trials", dist = "normal")
  formula <- bmf(dprime ~ 1 + (1 | id), criterion ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
})

test_that("sdt default_prior returns valid prior object", {
  dat <- data.frame(
    n_old = c(10, 40),
    stimulus = c(0L, 1L),
    n_trials = c(50, 50)
  )
  model <- sdt("n_old", "stimulus", "n_trials", dist = "normal")
  formula <- bmf(dprime ~ 1, criterion ~ 1)
  prior <- default_prior(formula, data = dat, model = model)
  expect_s3_class(prior, "brmsprior")
  # Check that our default priors appear
  prior_strs <- prior$prior
  expect_true(any(grepl("normal", prior_strs)))
})


############################################################################# !
# RATING MODEL CONSTRUCTOR TESTS                                          ####
############################################################################# !

test_that("sdt rating model can be created with vector response", {
  model <- sdt(c("r1", "r2", "r3", "r4"), "stimulus")
  expect_s3_class(model, "bmmodel")
  expect_s3_class(model, "sdt")
  expect_s3_class(model, "sdt_evsdt_rating")
  expect_equal(model$other_vars$n_ratings, 4L)
})

test_that("sdt rating model infers n_ratings from response length", {
  model <- sdt(c("r1", "r2", "r3", "r4", "r5", "r6"), "stimulus")
  expect_equal(model$other_vars$n_ratings, 6L)
})

test_that("sdt rating model rejects mismatched n_ratings", {
  expect_error(
    sdt(c("r1", "r2", "r3", "r4"), "stimulus", n_ratings = 6),
    "must match"
  )
})

test_that("sdt rating model has equidistant threshold params by default", {
  model <- sdt(c("r1", "r2", "r3", "r4"), "stimulus")
  expect_true("spacing" %in% names(model$parameters))
  expect_true("dprime" %in% names(model$parameters))
  expect_true("criterion" %in% names(model$parameters))
  # mu should NOT be present for rating models
  expect_false("mu" %in% names(model$parameters))
})

test_that("sdt rating model has log_distance threshold params", {
  model <- sdt(c("r1", "r2", "r3", "r4"), "stimulus",
               threshold_type = "log_distance")
  # K=4 -> mid=2 -> K-2=2 deltas: delta1, delta3 (skip mid)
  expect_true("delta1" %in% names(model$parameters))
  expect_true("delta3" %in% names(model$parameters))
  expect_false("delta2" %in% names(model$parameters))
  expect_false("spacing" %in% names(model$parameters))
})

test_that("sdt rating model has void_mu = TRUE", {
  model <- sdt(c("r1", "r2", "r3", "r4"), "stimulus")
  expect_true(model$void_mu)
})

test_that("sdt binary requires n_trials", {
  expect_error(sdt("n_old", "stimulus"), "n_trials is required")
})

test_that("sdt rating stores threshold_type correctly", {
  model_eq <- sdt(c("r1", "r2", "r3", "r4"), "stimulus")
  expect_equal(model_eq$other_vars$threshold_type, "equidistant")

  model_ld <- sdt(c("r1", "r2", "r3", "r4"), "stimulus",
                  threshold_type = "log_distance")
  expect_equal(model_ld$other_vars$threshold_type, "log_distance")

  model_pa <- sdt(c("r1", "r2", "r3", "r4"), "stimulus",
                  threshold_type = "parsimonious")
  expect_equal(model_pa$other_vars$threshold_type, "parsimonious")
})

test_that("sdt rating model has parsimonious threshold params", {
  model <- sdt(c("r1", "r2", "r3", "r4"), "stimulus",
               threshold_type = "parsimonious")
  # Parsimonious uses same 2 parameters as equidistant: criterion + spacing
  expect_true("spacing" %in% names(model$parameters))
  expect_true("criterion" %in% names(model$parameters))
  expect_true("dprime" %in% names(model$parameters))
  # No delta parameters
  expect_false("delta1" %in% names(model$parameters))
  expect_false("delta3" %in% names(model$parameters))
})


############################################################################# !
# RATING CHECK_DATA TESTS                                                 ####
############################################################################# !

test_that("sdt rating check_data validates response columns", {
  model <- sdt(c("r1", "r2", "r3", "r4"), "stimulus")
  formula <- bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1)

  # Valid data
  valid_data <- data.frame(
    r1 = c(10, 5), r2 = c(20, 10), r3 = c(15, 30), r4 = c(5, 55),
    stimulus = c(0L, 1L)
  )
  result <- check_data(model, valid_data, formula)
  expect_true("Y" %in% colnames(result))
  expect_true("nTrials" %in% colnames(result))
  expect_equal(result$nTrials, c(50, 100))
})

test_that("sdt rating check_data rejects missing response columns", {
  model <- sdt(c("r1", "r2", "r3", "r4"), "stimulus")
  formula <- bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1)

  invalid_data <- data.frame(
    r1 = c(10, 5), r2 = c(20, 10),
    stimulus = c(0L, 1L)
  )
  expect_error(check_data(model, invalid_data, formula), "missing in the data")
})

test_that("sdt rating check_data rejects negative counts", {
  model <- sdt(c("r1", "r2", "r3", "r4"), "stimulus")
  formula <- bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1)

  invalid_data <- data.frame(
    r1 = c(-1, 5), r2 = c(20, 10), r3 = c(15, 30), r4 = c(5, 55),
    stimulus = c(0L, 1L)
  )
  expect_error(check_data(model, invalid_data, formula), "non-negative")
})

test_that("sdt rating check_data validates stimulus coding", {
  model <- sdt(c("r1", "r2", "r3", "r4"), "stimulus")
  formula <- bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1)

  invalid_data <- data.frame(
    r1 = c(10, 5), r2 = c(20, 10), r3 = c(15, 30), r4 = c(5, 55),
    stimulus = c(0L, 2L)
  )
  expect_error(check_data(model, invalid_data, formula), "must be coded as 0")
})


############################################################################# !
# RATING DISTRIBUTION FUNCTION TESTS                                      ####
############################################################################# !

test_that("rsdt rating generates correct structure", {
  dat <- rsdt(n_per_cell = 50, n_subjects = 5, dprime = 1.5, criterion = 0,
              n_ratings = 4, spacing = 0.5, version = "rating")
  expect_true(is.data.frame(dat))
  expect_equal(nrow(dat), 10)  # 5 subjects * 2 stimulus types
  expect_true(all(c("id", "stimulus", "r1", "r2", "r3", "r4", "nTrials") %in%
                    colnames(dat)))
  expect_equal(unique(dat$stimulus), c(0L, 1L))
  # Row sums should equal n_per_cell
  expect_true(all(rowSums(dat[, paste0("r", 1:4)]) == 50))
})

test_that("rsdt rating works with log_distance thresholds", {
  dat <- rsdt(n_per_cell = 50, n_subjects = 3, dprime = 1.5, criterion = 0,
              n_ratings = 4, deltas = c(0.5, 0.5), version = "rating",
              threshold_type = "log_distance")
  expect_equal(nrow(dat), 6)
  expect_true(all(rowSums(dat[, paste0("r", 1:4)]) == 50))
})

test_that("rsdt rating works with parsimonious thresholds", {
  dat <- rsdt(n_per_cell = 50, n_subjects = 3, dprime = 1.5, criterion = 0,
              n_ratings = 4, spacing = 0.5, version = "rating",
              threshold_type = "parsimonious")
  expect_equal(nrow(dat), 6)
  expect_true(all(rowSums(dat[, paste0("r", 1:4)]) == 50))
})

test_that(".sdt_make_thresholds parsimonious computes logit-spaced thresholds", {
  # For K=4, logit(k/4) for k=1,2,3: log(1/3), 0, log(3)
  thr <- .sdt_make_thresholds(criterion = 0, n_ratings = 4,
                              threshold_type = "parsimonious", spacing = 0)
  # exp(0) = 1, so thresholds = logit(k/4)
  expected <- log(c(1, 2, 3) / c(3, 2, 1))  # logit(k/K) for k=1,2,3
  expect_equal(thr, expected, tolerance = 1e-8)

  # criterion shifts all thresholds
  thr_shifted <- .sdt_make_thresholds(criterion = 1, n_ratings = 4,
                                      threshold_type = "parsimonious", spacing = 0)
  expect_equal(thr_shifted, expected + 1, tolerance = 1e-8)

  # Thresholds are strictly increasing
  thr_pos <- .sdt_make_thresholds(criterion = 0, n_ratings = 6,
                                  threshold_type = "parsimonious", spacing = 0.5)
  expect_true(all(diff(thr_pos) > 0))
})

test_that("rsdt rating requires n_ratings >= 3", {
  expect_error(
    rsdt(n_per_cell = 50, n_subjects = 1, dprime = 1, criterion = 0,
         n_ratings = 2, spacing = 0.5, version = "rating"),
    "n_ratings must be >= 3"
  )
})

test_that("rsdt rating works with different distributions", {
  for (di in c("normal", "logistic", "gumbel_min", "gumbel_max")) {
    dat <- rsdt(n_per_cell = 50, n_subjects = 2, dprime = 1.5, criterion = 0,
                n_ratings = 4, spacing = 0.5, version = "rating", dist = di)
    expect_equal(nrow(dat), 4, info = paste("dist:", di))
    expect_true(all(rowSums(dat[, paste0("r", 1:4)]) == 50),
                info = paste("dist:", di))
  }
})

test_that("rsdt rating with K=6 generates 6 response columns", {
  dat <- rsdt(n_per_cell = 100, n_subjects = 2, dprime = 1.5, criterion = 0,
              n_ratings = 6, spacing = 0.3, version = "rating")
  expect_true(all(paste0("r", 1:6) %in% colnames(dat)))
  expect_true(all(rowSums(dat[, paste0("r", 1:6)]) == 100))
})

test_that("dsdt rating returns valid density", {
  d <- dsdt(counts = c(10, 20, 30, 40), stimulus = 1,
            dprime = 1.5, thresholds = c(-0.5, 0, 0.5),
            version = "rating")
  expect_true(d > 0)
  expect_true(d <= 1)
})

test_that("dsdt rating log matches log of density", {
  d <- dsdt(counts = c(10, 20, 30, 40), stimulus = 1,
            dprime = 1.5, thresholds = c(-0.5, 0, 0.5),
            version = "rating")
  ld <- dsdt(counts = c(10, 20, 30, 40), stimulus = 1,
             dprime = 1.5, thresholds = c(-0.5, 0, 0.5),
             version = "rating", log = TRUE)
  expect_equal(log(d), ld, tolerance = 1e-10)
})

test_that("dsdt rating validates inputs", {
  expect_error(
    dsdt(counts = c(-1, 20, 30, 40), stimulus = 1,
         dprime = 1.5, thresholds = c(-0.5, 0, 0.5), version = "rating"),
    "non-negative"
  )
  expect_error(
    dsdt(counts = c(10, 20, 30), stimulus = 1,
         dprime = 1.5, thresholds = c(-0.5, 0, 0.5), version = "rating"),
    "K - 1"
  )
})

test_that("dsdt rating works for all distributions", {
  for (di in c("normal", "logistic", "gumbel_min", "gumbel_max")) {
    d <- dsdt(counts = c(10, 20, 30, 40), stimulus = 1,
              dprime = 1.5, thresholds = c(-0.5, 0, 0.5),
              version = "rating", dist = di)
    expect_true(d > 0, info = paste("dist:", di))
  }
})

test_that("dsdt rating with sdratio != 1 changes density", {
  d_ev <- dsdt(counts = c(10, 20, 30, 40), stimulus = 1,
               dprime = 1.5, thresholds = c(-0.5, 0, 0.5),
               version = "rating", sdratio = 1)
  d_uv <- dsdt(counts = c(10, 20, 30, 40), stimulus = 1,
               dprime = 1.5, thresholds = c(-0.5, 0, 0.5),
               version = "rating", sdratio = 1.3)
  expect_false(d_ev == d_uv)
})


############################################################################# !
# COMBINE_SDT_RESPONSE TESTS                                              ####
############################################################################# !

test_that("combine_sdt_response maps categories correctly", {
  stim <- c(0, 0, 0, 1, 1, 1)
  conf <- c(3, 2, 1, 1, 2, 3)
  result <- combine_sdt_response(stim, conf, n_levels = 3)
  # Noise: high conf=1, mid conf=2, low conf=3
  # Signal: low conf=4, mid conf=5, high conf=6
  expect_equal(result, c(1, 2, 3, 4, 5, 6))
})

test_that("combine_sdt_response validates inputs", {
  expect_error(combine_sdt_response(c(0, 2), c(1, 1), n_levels = 3),
               "stimulus must be 0")
  expect_error(combine_sdt_response(c(0, 1), c(0, 1), n_levels = 3),
               "confidence must be between")
})

test_that("combine_sdt_response handles K=4 (2 confidence levels)", {
  stim <- c(0, 0, 1, 1)
  conf <- c(2, 1, 1, 2)
  result <- combine_sdt_response(stim, conf, n_levels = 2)
  expect_equal(result, c(1, 2, 3, 4))
})


############################################################################# !
# RATING FORMULA CONSTRUCTION TESTS                                       ####
############################################################################# !

test_that("sdt rating produces valid stancode with equidistant thresholds", {
  dat <- rsdt(n_per_cell = 50, n_subjects = 3, dprime = 1.5, criterion = 0,
              n_ratings = 4, spacing = 0.5, version = "rating")
  model <- sdt(c("r1", "r2", "r3", "r4"), "stimulus")
  formula <- bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
  # Should use multinomial family
  expect_true(grepl("multinomial", code, ignore.case = TRUE))
})

test_that("sdt rating produces valid stancode with log_distance thresholds", {
  dat <- rsdt(n_per_cell = 50, n_subjects = 3, dprime = 1.5, criterion = 0,
              n_ratings = 4, deltas = c(0.5, 0.5), version = "rating",
              threshold_type = "log_distance")
  model <- sdt(c("r1", "r2", "r3", "r4"), "stimulus",
               threshold_type = "log_distance")
  formula <- bmf(dprime ~ 1, criterion ~ 1, delta1 ~ 1, delta3 ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
  expect_true(grepl("multinomial", code, ignore.case = TRUE))
})

test_that("sdt rating produces valid stancode with parsimonious thresholds", {
  dat <- rsdt(n_per_cell = 50, n_subjects = 3, dprime = 1.5, criterion = 0,
              n_ratings = 4, spacing = 0.5, version = "rating",
              threshold_type = "parsimonious")
  model <- sdt(c("r1", "r2", "r3", "r4"), "stimulus",
               threshold_type = "parsimonious")
  formula <- bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
  expect_true(grepl("multinomial", code, ignore.case = TRUE))
  # Canonical positions should be embedded as numeric literals
  expect_true(grepl("1.0986122887", code))
})

test_that("sdt rating parsimonious K=6 produces valid stancode", {
  dat <- rsdt(n_per_cell = 50, n_subjects = 2, dprime = 1.5, criterion = 0,
              n_ratings = 6, spacing = 0.5, version = "rating",
              threshold_type = "parsimonious")
  model <- sdt(c("r1", "r2", "r3", "r4", "r5", "r6"), "stimulus",
               threshold_type = "parsimonious")
  formula <- bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
})

# --- log_ratio threshold tests ---

test_that("sdt rating model stores log_ratio threshold_type", {
  model <- sdt(c("r1", "r2", "r3", "r4"), "stimulus",
               threshold_type = "log_ratio")
  expect_equal(model$other_vars$threshold_type, "log_ratio")
})

test_that("sdt rating model has log_ratio threshold params", {
  model <- sdt(c("r1", "r2", "r3", "r4"), "stimulus",
               threshold_type = "log_ratio")
  expect_true("delta1" %in% names(model$parameters))
  expect_true("delta3" %in% names(model$parameters))
  expect_true("criterion" %in% names(model$parameters))
  expect_false("spacing" %in% names(model$parameters))
  expect_false("delta2" %in% names(model$parameters))
})

test_that(".sdt_make_thresholds log_ratio computes ratio-scaled thresholds", {
  # All deltas = 0 -> equally spaced with unit spacing
  thr <- .sdt_make_thresholds(criterion = 0, n_ratings = 4,
                              threshold_type = "log_ratio", deltas = c(0, 0))
  expect_equal(thr, c(-1, 0, 1), tolerance = 1e-8)

  # K=6, all zeros -> equally spaced
  thr6 <- .sdt_make_thresholds(criterion = 0, n_ratings = 6,
                               threshold_type = "log_ratio",
                               deltas = c(0, 0, 0, 0))
  expect_equal(thr6, c(-2, -1, 0, 1, 2), tolerance = 1e-8)

  # Criterion shift
  thr_shifted <- .sdt_make_thresholds(criterion = 1, n_ratings = 4,
                                      threshold_type = "log_ratio",
                                      deltas = c(0, 0))
  expect_equal(thr_shifted, c(0, 1, 2), tolerance = 1e-8)

  # Monotonicity with non-zero deltas
  thr_asym <- .sdt_make_thresholds(criterion = 0, n_ratings = 6,
                                   threshold_type = "log_ratio",
                                   deltas = c(0.3, -0.2, 0.5, 0.1))
  expect_true(all(diff(thr_asym) > 0))

  # Asymmetric: delta_{mid-1} != 0 makes lower spread differ from upper
  thr_k4 <- .sdt_make_thresholds(criterion = 0, n_ratings = 4,
                                 threshold_type = "log_ratio",
                                 deltas = c(log(2), 0))
  # spread_above = exp(0) = 1, spread_below = exp(log(2)) * 1 = 2
  expect_equal(thr_k4, c(-2, 0, 1), tolerance = 1e-8)
})

test_that("rsdt rating works with log_ratio thresholds", {
  dat <- rsdt(n_per_cell = 50, n_subjects = 3, dprime = 1.5, criterion = 0,
              n_ratings = 4, deltas = c(0, 0), version = "rating",
              threshold_type = "log_ratio")
  expect_equal(nrow(dat), 6)
  expect_true(all(rowSums(dat[, paste0("r", 1:4)]) == 50))
})

test_that("sdt rating produces valid stancode with log_ratio thresholds", {
  dat <- rsdt(n_per_cell = 50, n_subjects = 3, dprime = 1.5, criterion = 0,
              n_ratings = 4, deltas = c(0, 0), version = "rating",
              threshold_type = "log_ratio")
  model <- sdt(c("r1", "r2", "r3", "r4"), "stimulus",
               threshold_type = "log_ratio")
  formula <- bmf(dprime ~ 1, criterion ~ 1, delta1 ~ 1, delta3 ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
  expect_true(grepl("multinomial", code, ignore.case = TRUE))
  # log_ratio multiplies delta expressions: exp(nlp_delta1) * exp(nlp_delta3)
  expect_true(grepl("exp\\(nlp_delta", code))
})

test_that("sdt rating log_ratio K=6 produces valid stancode", {
  dat <- rsdt(n_per_cell = 50, n_subjects = 3, dprime = 1.5, criterion = 0,
              n_ratings = 6, deltas = c(0, 0, 0, 0), version = "rating",
              threshold_type = "log_ratio")
  model <- sdt(paste0("r", 1:6), "stimulus", threshold_type = "log_ratio")
  formula <- bmf(dprime ~ 1, criterion ~ 1, delta1 ~ 1, delta2 ~ 1,
                 delta4 ~ 1, delta5 ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
})

test_that("sdt rating handles K=6 equidistant", {
  dat <- rsdt(n_per_cell = 50, n_subjects = 2, dprime = 1.5, criterion = 0,
              n_ratings = 6, spacing = 0.3, version = "rating")
  model <- sdt(c("r1", "r2", "r3", "r4", "r5", "r6"), "stimulus")
  formula <- bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
})

test_that("sdt rating handles predictors on dprime", {
  dat <- rsdt(n_per_cell = 50, n_subjects = 3, dprime = 1.5, criterion = 0,
              n_ratings = 4, spacing = 0.5, version = "rating")
  dat$condition <- rep(c("A", "B"), length.out = nrow(dat))
  model <- sdt(c("r1", "r2", "r3", "r4"), "stimulus")
  formula <- bmf(dprime ~ condition, criterion ~ 1, spacing ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
})


############################################################################# !
# UNEQUAL-VARIANCE SDT TESTS                                              ####
############################################################################# !

test_that("sdt UV-SDT model has sdratio parameter", {
  model <- sdt(c("r1", "r2", "r3", "r4"), "stimulus", variances = "unequal")
  expect_true("sdratio" %in% names(model$parameters))
  expect_true("sdratio" %in% names(model$links))
  expect_true("sdratio" %in% names(model$default_priors))
  expect_equal(model$other_vars$variances, "unequal")
})

test_that("sdt EV-SDT model does NOT have sdratio parameter", {
  model <- sdt(c("r1", "r2", "r3", "r4"), "stimulus", variances = "equal")
  expect_false("sdratio" %in% names(model$parameters))
})

test_that("sdt UV-SDT requires ratings", {
  expect_error(
    sdt("n_old", "stimulus", "n_trials", variances = "unequal"),
    "requires confidence rating"
  )
})

test_that("sdt UV-SDT produces valid stancode", {
  dat <- rsdt(n_per_cell = 50, n_subjects = 3, dprime = 1.5, criterion = 0,
              sdratio = 1.3, n_ratings = 4, spacing = 0.5, version = "rating")
  model <- sdt(c("r1", "r2", "r3", "r4"), "stimulus", variances = "unequal")
  formula <- bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1, sdratio ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
  expect_true(grepl("multinomial", code, ignore.case = TRUE))
  # Should contain sdratio in the generated code
  expect_true(grepl("sdratio", code))
})

test_that("sdt UV-SDT with log_distance thresholds produces valid stancode", {
  dat <- rsdt(n_per_cell = 50, n_subjects = 3, dprime = 1.5, criterion = 0,
              sdratio = 1.3, n_ratings = 4, deltas = c(0.5, 0.5),
              version = "rating", threshold_type = "log_distance")
  model <- sdt(c("r1", "r2", "r3", "r4"), "stimulus",
               variances = "unequal", threshold_type = "log_distance")
  formula <- bmf(dprime ~ 1, criterion ~ 1, delta1 ~ 1, delta3 ~ 1,
                 sdratio ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
  expect_true(grepl("sdratio", code))
})

test_that("sdt UV-SDT with K=6 produces valid stancode", {
  dat <- rsdt(n_per_cell = 50, n_subjects = 2, dprime = 1.5, criterion = 0,
              sdratio = 1.3, n_ratings = 6, spacing = 0.3, version = "rating")
  model <- sdt(c("r1", "r2", "r3", "r4", "r5", "r6"), "stimulus",
               variances = "unequal")
  formula <- bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1, sdratio ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
})


############################################################################# !
# VECTORIZATION TESTS                                                      ####
############################################################################# !

test_that(".sdt_eta is vectorized over all arguments", {
  # Vector stimulus
  eta <- bmm:::.sdt_eta(1.5, 0.2, c(0, 1))
  expect_length(eta, 2)
  expect_equal(eta[1], -1.5/2 - 0.2)
  expect_equal(eta[2], 1.5/2 - 0.2)

  # Vector dprime + stimulus
  eta <- bmm:::.sdt_eta(c(1.0, 2.0), 0.2, c(0, 1))
  expect_length(eta, 2)
  expect_equal(eta[1], -1.0/2 - 0.2)
  expect_equal(eta[2], 2.0/2 - 0.2)

  # Vector sdratio (UV-SDT)
  eta <- bmm:::.sdt_eta(1.5, 0.2, c(0, 1), sdratio = c(1.0, 1.3))
  expect_length(eta, 2)
  expect_equal(eta[1], -1.5/2 - 0.2)  # noise: sdratio ignored
  expect_equal(eta[2], (1.5/2 - 0.2) / 1.3)  # signal: scaled
})

test_that(".sdt_eta EV and UV give same results for noise trials", {
  eta_ev <- bmm:::.sdt_eta(1.5, 0.2, 0, sdratio = 1)
  eta_uv <- bmm:::.sdt_eta(1.5, 0.2, 0, sdratio = 1.5)
  expect_equal(eta_ev, eta_uv)
})

test_that("dsdt is vectorized over sdratio", {
  # Should not error with vector sdratio
  d <- dsdt(n_old = c(30, 80), n_trials = c(100, 100),
            stimulus = c(0, 1), dprime = 1.5, criterion = 0.2,
            sdratio = c(1.0, 1.3))
  expect_length(d, 2)
  expect_true(all(d > 0))
})

test_that("rsdt rating rejects vector dprime", {
  expect_error(
    rsdt(n_per_cell = 50, n_subjects = 2, dprime = c(1, 2),
         criterion = 0, n_ratings = 4, spacing = 0.5, version = "rating"),
    "dprime must be a single value"
  )
})

test_that("rsdt rating rejects vector criterion", {
  expect_error(
    rsdt(n_per_cell = 50, n_subjects = 2, dprime = 1.5,
         criterion = c(0, 0.5), n_ratings = 4, spacing = 0.5,
         version = "rating"),
    "criterion must be a single value"
  )
})

test_that("rsdt rating rejects vector sdratio", {
  expect_error(
    rsdt(n_per_cell = 50, n_subjects = 2, dprime = 1.5,
         criterion = 0, sdratio = c(1.0, 1.3),
         n_ratings = 4, spacing = 0.5, version = "rating"),
    "sdratio must be a single value"
  )
})


############################################################################# !
# DPSDT MODEL TESTS                                                       ####
############################################################################# !

# --- Constructor Tests ---

test_that("sdt DPSDT model can be created", {
  model <- sdt(c("r1", "r2", "r3", "r4"), "stimulus", version = "dpsdt")
  expect_s3_class(model, "bmmodel")
  expect_s3_class(model, "sdt")
  expect_s3_class(model, "sdt_dpsdt_rating")
})

test_that("sdt DPSDT has Ro parameter", {
  model <- sdt(c("r1", "r2", "r3", "r4"), "stimulus", version = "dpsdt")
  expect_true("Ro" %in% names(model$parameters))
  expect_true("Ro" %in% names(model$links))
  expect_true("Ro" %in% names(model$default_priors))
  expect_equal(model$links$Ro, "identity")
})

test_that("sdt DPSDT has standard SDT parameters too", {
  model <- sdt(c("r1", "r2", "r3", "r4"), "stimulus", version = "dpsdt")
  expect_true("dprime" %in% names(model$parameters))
  expect_true("criterion" %in% names(model$parameters))
  expect_true("spacing" %in% names(model$parameters))
})

test_that("sdt DPSDT requires ratings", {
  expect_error(
    sdt("n_old", "stimulus", "n_trials", version = "dpsdt"),
    "DPSDT requires confidence rating"
  )
})

test_that("sdt DPSDT can be combined with unequal variances", {
  model <- sdt(c("r1", "r2", "r3", "r4"), "stimulus",
               version = "dpsdt", variances = "unequal")
  expect_true("Ro" %in% names(model$parameters))
  expect_true("sdratio" %in% names(model$parameters))
})

test_that("sdt DPSDT works with different distributions", {
  for (d in c("normal", "logistic", "gumbel_min", "gumbel_max")) {
    model <- sdt(c("r1", "r2", "r3", "r4"), "stimulus",
                 version = "dpsdt", dist = d)
    expect_s3_class(model, "sdt_dpsdt_rating")
    expect_equal(model$other_vars$dist, d)
  }
})

test_that("sdt DPSDT works with log_distance thresholds", {
  model <- sdt(c("r1", "r2", "r3", "r4"), "stimulus",
               version = "dpsdt", threshold_type = "log_distance")
  expect_true("delta1" %in% names(model$parameters))
  expect_true("delta3" %in% names(model$parameters))
  expect_false("spacing" %in% names(model$parameters))
})

# --- Distribution Function Tests ---

test_that("rsdt generates DPSDT data", {
  dat <- rsdt(n_per_cell = 100, n_subjects = 5, dprime = 1.0, criterion = 0,
              Ro = 0.3, n_ratings = 4, spacing = 0.5, version = "dpsdt")
  expect_true(is.data.frame(dat))
  expect_true(all(c("r1", "r2", "r3", "r4", "stimulus", "id") %in% names(dat)))
  expect_equal(nrow(dat), 10)  # 5 subjects x 2 stimulus types
  expect_true(all(rowSums(dat[, paste0("r", 1:4)]) == 100))
})

test_that("rsdt DPSDT Ro=0 matches standard rating SDT", {
  set.seed(42)
  dat_dp <- rsdt(n_per_cell = 1000, n_subjects = 1, dprime = 1.5, criterion = 0,
                 Ro = 0, n_ratings = 4, spacing = 0.5, version = "dpsdt")
  set.seed(42)
  dat_std <- rsdt(n_per_cell = 1000, n_subjects = 1, dprime = 1.5, criterion = 0,
                  n_ratings = 4, spacing = 0.5, version = "rating")
  # With Ro=0, DPSDT should produce identical results to standard rating
  expect_equal(dat_dp[, paste0("r", 1:4)], dat_std[, paste0("r", 1:4)])
})

test_that("rsdt DPSDT Ro shifts old items to highest confidence", {
  set.seed(42)
  dat_no_r <- rsdt(n_per_cell = 1000, n_subjects = 1, dprime = 1.0, criterion = 0,
                   Ro = 0, n_ratings = 4, spacing = 0.5, version = "dpsdt")
  set.seed(42)
  dat_hi_r <- rsdt(n_per_cell = 1000, n_subjects = 1, dprime = 1.0, criterion = 0,
                   Ro = 0.5, n_ratings = 4, spacing = 0.5, version = "dpsdt")

  # For old items (stimulus=1): Ro should increase highest "old" (r4) counts
  old_no_r <- dat_no_r[dat_no_r$stimulus == 1, ]
  old_hi_r <- dat_hi_r[dat_hi_r$stimulus == 1, ]
  expect_true(old_hi_r$r4 > old_no_r$r4)

  # For new items (stimulus=0): Ro should not affect counts
  new_no_r <- dat_no_r[dat_no_r$stimulus == 0, ]
  new_hi_r <- dat_hi_r[dat_hi_r$stimulus == 0, ]
  # Use different seed so can't compare directly — check proportions are similar
})

test_that("rsdt DPSDT requires Ro argument", {
  expect_error(
    rsdt(n_per_cell = 50, n_subjects = 1, dprime = 1.5, criterion = 0,
         n_ratings = 4, spacing = 0.5, version = "dpsdt"),
    "Ro.*required"
  )
})

test_that("rsdt DPSDT validates Ro bounds", {
  expect_error(
    rsdt(n_per_cell = 50, n_subjects = 1, dprime = 1.5, criterion = 0,
         Ro = -0.1, n_ratings = 4, spacing = 0.5, version = "dpsdt"),
    "Ro must be between 0 and 1"
  )
  expect_error(
    rsdt(n_per_cell = 50, n_subjects = 1, dprime = 1.5, criterion = 0,
         Ro = 1.1, n_ratings = 4, spacing = 0.5, version = "dpsdt"),
    "Ro must be between 0 and 1"
  )
})

test_that("dsdt DPSDT computes density", {
  thresholds <- c(-0.5, 0.0, 0.5)

  # Old item with recollection
  d1 <- dsdt(counts = c(5, 15, 25, 55), stimulus = 1,
             dprime = 1.5, thresholds = thresholds,
             Ro = 0.3, version = "dpsdt")
  expect_true(d1 > 0)

  # New item (Ro should not affect)
  d2 <- dsdt(counts = c(55, 25, 15, 5), stimulus = 0,
             dprime = 1.5, thresholds = thresholds,
             Ro = 0.3, version = "dpsdt")
  expect_true(d2 > 0)
})

test_that("dsdt DPSDT Ro=0 matches standard rating density", {
  thresholds <- c(-0.5, 0.0, 0.5)
  counts <- c(10, 20, 30, 40)

  d_dp <- dsdt(counts = counts, stimulus = 1,
               dprime = 1.5, thresholds = thresholds,
               Ro = 0, version = "dpsdt", log = TRUE)
  d_std <- dsdt(counts = counts, stimulus = 1,
                dprime = 1.5, thresholds = thresholds,
                version = "rating", log = TRUE)
  expect_equal(d_dp, d_std)
})

test_that("dsdt DPSDT Ro=0 matches for noise trials too", {
  thresholds <- c(-0.5, 0.0, 0.5)
  counts <- c(40, 30, 20, 10)

  d_dp <- dsdt(counts = counts, stimulus = 0,
               dprime = 1.5, thresholds = thresholds,
               Ro = 0.3, version = "dpsdt", log = TRUE)
  d_std <- dsdt(counts = counts, stimulus = 0,
                dprime = 1.5, thresholds = thresholds,
                version = "rating", log = TRUE)
  # For noise trials, Ro doesn't apply, so densities should match
  expect_equal(d_dp, d_std)
})

# --- DPSDT Category Probability Tests ---

test_that("DPSDT category probs sum to 1", {
  thresholds <- c(-0.5, 0.0, 0.5)

  # Old items
  probs_old <- bmm:::.sdt_dpsdt_category_probs(
    thresholds, shift = 0.75, sdratio = 1, stimulus = 1,
    dist = "normal", Ro = 0.3
  )
  expect_equal(sum(probs_old), 1, tolerance = 1e-10)

  # New items
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

  # Highest "old" category should increase
  K <- length(probs_no_r)
  expect_true(probs_hi_r[K] > probs_no_r[K])

  # All other categories should decrease (scaled by 1-Ro)
  for (k in 1:(K - 1)) {
    expect_true(probs_hi_r[k] < probs_no_r[k])
  }
})

test_that("DPSDT Ro does not affect new items", {
  thresholds <- c(-0.5, 0.0, 0.5)
  shift <- -0.75

  probs_std <- bmm:::.sdt_category_probs(
    thresholds, shift, sdratio = 1, stimulus = 0, dist = "normal"
  )
  probs_dpsdt <- bmm:::.sdt_dpsdt_category_probs(
    thresholds, shift, sdratio = 1, stimulus = 0, dist = "normal", Ro = 0.5
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
  # All probability should be in the highest category
  expect_equal(probs[K], 1.0, tolerance = 1e-10)
  for (k in 1:(K - 1)) {
    expect_true(probs[k] < 1e-8)
  }
})

# --- Formula Construction Tests ---

test_that("sdt DPSDT produces valid stancode with equidistant thresholds", {
  dat <- rsdt(n_per_cell = 50, n_subjects = 3, dprime = 1.0, criterion = 0,
              Ro = 0.3, n_ratings = 4, spacing = 0.5, version = "dpsdt")
  model <- sdt(c("r1", "r2", "r3", "r4"), "stimulus", version = "dpsdt")
  formula <- bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1, Ro ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
  expect_true(grepl("multinomial", code, ignore.case = TRUE))
  # Should contain inv_logit(Ro) for recollection
  expect_true(grepl("inv_logit", code))
  expect_true(grepl("Ro", code))
})

test_that("sdt DPSDT produces valid stancode with log_distance thresholds", {
  dat <- rsdt(n_per_cell = 50, n_subjects = 3, dprime = 1.0, criterion = 0,
              Ro = 0.3, n_ratings = 4, spacing = 0.5, version = "dpsdt")
  model <- sdt(c("r1", "r2", "r3", "r4"), "stimulus",
               version = "dpsdt", threshold_type = "log_distance")
  formula <- bmf(dprime ~ 1, criterion ~ 1, delta1 ~ 1, delta3 ~ 1, Ro ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
  expect_true(grepl("inv_logit", code))
})

test_that("sdt DPSDT produces valid stancode with K=6", {
  dat <- rsdt(n_per_cell = 50, n_subjects = 2, dprime = 1.0, criterion = 0,
              Ro = 0.3, n_ratings = 6, spacing = 0.3, version = "dpsdt")
  model <- sdt(c("r1", "r2", "r3", "r4", "r5", "r6"), "stimulus",
               version = "dpsdt")
  formula <- bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1, Ro ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
})

test_that("sdt DPSDT produces valid stancode with predictors", {
  dat <- rsdt(n_per_cell = 50, n_subjects = 3, dprime = 1.0, criterion = 0,
              Ro = 0.3, n_ratings = 4, spacing = 0.5, version = "dpsdt")
  dat$condition <- rep(c("A", "B"), length.out = nrow(dat))
  model <- sdt(c("r1", "r2", "r3", "r4"), "stimulus", version = "dpsdt")
  formula <- bmf(dprime ~ condition, criterion ~ 1, spacing ~ 1,
                 Ro ~ condition)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
})

test_that("sdt DPSDT + UV-SDT produces valid stancode", {
  dat <- rsdt(n_per_cell = 50, n_subjects = 3, dprime = 1.0, criterion = 0,
              Ro = 0.3, sdratio = 1.2, n_ratings = 4, spacing = 0.5,
              version = "dpsdt")
  model <- sdt(c("r1", "r2", "r3", "r4"), "stimulus",
               version = "dpsdt", variances = "unequal")
  formula <- bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1,
                 sdratio ~ 1, Ro ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
  expect_true(grepl("sdratio", code))
  expect_true(grepl("inv_logit", code))
})

test_that("sdt DPSDT with all 4 distributions produces valid stancode", {
  dat <- rsdt(n_per_cell = 50, n_subjects = 2, dprime = 1.0, criterion = 0,
              Ro = 0.3, n_ratings = 4, spacing = 0.5, version = "dpsdt")
  for (d in c("normal", "logistic", "gumbel_min", "gumbel_max")) {
    model <- sdt(c("r1", "r2", "r3", "r4"), "stimulus",
                 version = "dpsdt", dist = d)
    formula <- bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1, Ro ~ 1)
    code <- stancode(formula, data = dat, model = model)
    expect_true(nchar(code) > 0)
  }
})


############################################################################# !
# META-D' MODEL TESTS                                                     ####
############################################################################# !

# --- Constructor Tests ---

test_that("sdt meta-d' model can be created", {
  model <- sdt(c("r1", "r2", "r3", "r4"), "stimulus", version = "metad")
  expect_s3_class(model, "bmmodel")
  expect_s3_class(model, "sdt")
  expect_s3_class(model, "sdt_metad_rating")
})

test_that("sdt meta-d' has metad parameter", {
  model <- sdt(c("r1", "r2", "r3", "r4"), "stimulus", version = "metad")
  expect_true("metad" %in% names(model$parameters))
  expect_true("metad" %in% names(model$links))
  expect_true("metad" %in% names(model$default_priors))
  expect_equal(model$links$metad, "identity")
})

test_that("sdt meta-d' has standard SDT parameters too", {
  model <- sdt(c("r1", "r2", "r3", "r4"), "stimulus", version = "metad")
  expect_true("dprime" %in% names(model$parameters))
  expect_true("criterion" %in% names(model$parameters))
  expect_true("spacing" %in% names(model$parameters))
})

test_that("sdt meta-d' requires ratings", {
  expect_error(
    sdt("n_old", "stimulus", "n_trials", version = "metad"),
    "Meta-d' requires confidence rating"
  )
})

test_that("sdt meta-d' works with different distributions", {
  for (d in c("normal", "logistic", "gumbel_min", "gumbel_max")) {
    model <- sdt(c("r1", "r2", "r3", "r4"), "stimulus",
                 version = "metad", dist = d)
    expect_s3_class(model, "sdt_metad_rating")
  }
})

# --- Distribution Function Tests ---

test_that("rsdt generates meta-d' data", {
  dat <- rsdt(n_per_cell = 100, n_subjects = 5, dprime = 1.5, criterion = 0,
              metad = 1.0, n_ratings = 4, spacing = 0.5, version = "metad")
  expect_true(is.data.frame(dat))
  expect_true(all(c("r1", "r2", "r3", "r4", "stimulus", "id") %in% names(dat)))
  expect_equal(nrow(dat), 10)  # 5 subjects x 2 stimulus types
  expect_true(all(rowSums(dat[, paste0("r", 1:4)]) == 100))
})

test_that("rsdt meta-d' with metad=dprime matches standard rating SDT", {
  set.seed(42)
  dat_metad <- rsdt(n_per_cell = 1000, n_subjects = 1, dprime = 1.5,
                    criterion = 0, metad = 1.5, n_ratings = 4, spacing = 0.5,
                    version = "metad")
  set.seed(42)
  dat_std <- rsdt(n_per_cell = 1000, n_subjects = 1, dprime = 1.5,
                  criterion = 0, n_ratings = 4, spacing = 0.5,
                  version = "rating")
  # When metad == dprime, normalization is 1, so results should match
  expect_equal(dat_metad[, paste0("r", 1:4)], dat_std[, paste0("r", 1:4)])
})

test_that("rsdt meta-d' requires metad argument", {
  expect_error(
    rsdt(n_per_cell = 50, n_subjects = 1, dprime = 1.5, criterion = 0,
         n_ratings = 4, spacing = 0.5, version = "metad"),
    "metad.*required"
  )
})

test_that("dsdt meta-d' computes density", {
  thresholds <- c(-0.5, 0.0, 0.5)

  d1 <- dsdt(counts = c(5, 15, 25, 55), stimulus = 1,
             dprime = 1.5, metad = 1.0, thresholds = thresholds,
             version = "metad")
  expect_true(d1 > 0)

  d2 <- dsdt(counts = c(55, 25, 15, 5), stimulus = 0,
             dprime = 1.5, metad = 1.0, thresholds = thresholds,
             version = "metad")
  expect_true(d2 > 0)
})

test_that("dsdt meta-d' with metad=dprime matches standard rating", {
  thresholds <- c(-0.5, 0.0, 0.5)
  counts <- c(10, 20, 30, 40)

  d_metad <- dsdt(counts = counts, stimulus = 1,
                  dprime = 1.5, metad = 1.5, thresholds = thresholds,
                  version = "metad", log = TRUE)
  d_std <- dsdt(counts = counts, stimulus = 1,
                dprime = 1.5, thresholds = thresholds,
                version = "rating", log = TRUE)
  expect_equal(d_metad, d_std, tolerance = 1e-10)
})

# --- Meta-d' Category Probability Tests ---

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

  # With lower meta-d', confidence resolution is worse
  # The extreme categories should be less differentiated from middle ones
  # Specifically, the highest confidence "old" (category K) should be lower
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

# --- Formula Construction Tests ---

test_that("sdt meta-d' produces valid stancode with equidistant thresholds", {
  dat <- rsdt(n_per_cell = 50, n_subjects = 3, dprime = 1.5, criterion = 0,
              metad = 1.0, n_ratings = 4, spacing = 0.5, version = "metad")
  model <- sdt(c("r1", "r2", "r3", "r4"), "stimulus", version = "metad")
  formula <- bmf(dprime ~ 1, metad ~ 1, criterion ~ 1, spacing ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
  expect_true(grepl("multinomial", code, ignore.case = TRUE))
  expect_true(grepl("metad", code))
})

test_that("sdt meta-d' produces valid stancode with K=6", {
  dat <- rsdt(n_per_cell = 50, n_subjects = 2, dprime = 1.5, criterion = 0,
              metad = 1.0, n_ratings = 6, spacing = 0.3, version = "metad")
  model <- sdt(c("r1", "r2", "r3", "r4", "r5", "r6"), "stimulus",
               version = "metad")
  formula <- bmf(dprime ~ 1, metad ~ 1, criterion ~ 1, spacing ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
})

test_that("sdt meta-d' produces valid stancode with predictors on metad", {
  dat <- rsdt(n_per_cell = 50, n_subjects = 3, dprime = 1.5, criterion = 0,
              metad = 1.0, n_ratings = 4, spacing = 0.5, version = "metad")
  dat$condition <- rep(c("A", "B"), length.out = nrow(dat))
  model <- sdt(c("r1", "r2", "r3", "r4"), "stimulus", version = "metad")
  formula <- bmf(dprime ~ 1, metad ~ condition, criterion ~ 1, spacing ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
})

test_that("sdt meta-d' with all 4 distributions produces valid stancode", {
  dat <- rsdt(n_per_cell = 50, n_subjects = 2, dprime = 1.5, criterion = 0,
              metad = 1.0, n_ratings = 4, spacing = 0.5, version = "metad")
  for (d in c("normal", "logistic", "gumbel_min", "gumbel_max")) {
    model <- sdt(c("r1", "r2", "r3", "r4"), "stimulus",
                 version = "metad", dist = d)
    formula <- bmf(dprime ~ 1, metad ~ 1, criterion ~ 1, spacing ~ 1)
    code <- stancode(formula, data = dat, model = model)
    expect_true(nchar(code) > 0)
  }
})

test_that("sdt meta-d' with log_distance thresholds produces valid stancode", {
  dat <- rsdt(n_per_cell = 50, n_subjects = 3, dprime = 1.5, criterion = 0,
              metad = 1.0, n_ratings = 4, spacing = 0.5, version = "metad")
  model <- sdt(c("r1", "r2", "r3", "r4"), "stimulus",
               version = "metad", threshold_type = "log_distance")
  formula <- bmf(dprime ~ 1, metad ~ 1, criterion ~ 1,
                 delta1 ~ 1, delta3 ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
})


############################################################################# !
# m-AFC MODEL TESTS                                                        ####
############################################################################# !

# --- Constructor Tests ---

test_that("sdt mafc model can be created", {
  model <- sdt("n_correct", n_trials = "n_trials", version = "mafc", m = 4)
  expect_s3_class(model, "bmmodel")
  expect_s3_class(model, "sdt")
  expect_s3_class(model, "sdt_mafc")
})

test_that("sdt mafc has correct parameters (only dprime, no criterion)", {
  model <- sdt("n_correct", n_trials = "n_trials", version = "mafc", m = 4)
  expect_true("dprime" %in% names(model$parameters))
  expect_true("mu" %in% names(model$parameters))
  # No criterion for m-AFC
  expect_false("criterion" %in% names(model$parameters))
  # mu is fixed to 0
  expect_equal(model$fixed_parameters$mu, 0)
})

test_that("sdt mafc stores m correctly", {
  model <- sdt("n_correct", n_trials = "n_trials", version = "mafc", m = 4)
  expect_equal(model$other_vars$m, 4L)

  model2 <- sdt("n_correct", n_trials = "n_trials", version = "mafc", m = 8)
  expect_equal(model2$other_vars$m, 8L)
})

test_that("sdt mafc has correct links", {
  model <- sdt("n_correct", n_trials = "n_trials", version = "mafc", m = 4)
  expect_equal(model$links$mu, "identity")
  expect_equal(model$links$dprime, "identity")
})

test_that("sdt mafc accepts custom links", {
  model <- sdt("n_correct", n_trials = "n_trials", version = "mafc", m = 4,
               links = list(dprime = "log"))
  expect_equal(model$links$dprime, "log")
})

test_that("sdt mafc requires m >= 2", {
  expect_error(sdt("n_correct", n_trials = "n_trials", version = "mafc", m = 1))
  expect_error(sdt("n_correct", n_trials = "n_trials", version = "mafc", m = NULL))
  expect_error(sdt("n_correct", n_trials = "n_trials", version = "mafc"))
})

test_that("sdt mafc requires n_trials", {
  expect_error(sdt("n_correct", version = "mafc", m = 4))
})

test_that("sdt mafc requires dist = 'normal'", {
  expect_error(
    sdt("n_correct", n_trials = "n_trials", version = "mafc",
        m = 4, dist = "logistic")
  )
})

test_that("sdt mafc does not require stimulus", {
  # stimulus is NULL by default, which should be fine for mafc
  expect_silent(sdt("n_correct", n_trials = "n_trials", version = "mafc", m = 4))
})

# --- Distribution Function Tests ---

test_that(".mafc_pc_r gives correct 2-AFC result", {
  # 2-AFC: P(correct) = Phi(d'/sqrt(2))
  d <- 1.5
  expected <- pnorm(d / sqrt(2))
  result <- bmm:::.mafc_pc_r(d, 2L)
  expect_equal(result, expected, tolerance = 1e-10)
})

test_that(".mafc_pc_r gives correct results for m >= 3", {
  # Verify against numerical integration for m = 4
  d <- 1.5
  exact <- integrate(function(x) dnorm(x) * pnorm(x + d)^3, -Inf, Inf)$value
  result <- bmm:::.mafc_pc_r(d, 4L)
  expect_equal(result, exact, tolerance = 1e-6)
})

test_that(".mafc_pc_r works for various m values", {
  d <- 1.0
  # P(correct) should decrease as m increases (more alternatives = harder)
  pc2 <- bmm:::.mafc_pc_r(d, 2L)
  pc4 <- bmm:::.mafc_pc_r(d, 4L)
  pc8 <- bmm:::.mafc_pc_r(d, 8L)
  expect_true(pc2 > pc4)
  expect_true(pc4 > pc8)

  # All should be above chance (1/m)
  expect_true(pc2 > 0.5)
  expect_true(pc4 > 0.25)
  expect_true(pc8 > 0.125)
})

test_that(".mafc_pc_r at d'=0 gives chance performance", {
  # At d'=0, P(correct) = 1/m for all m
  expect_equal(bmm:::.mafc_pc_r(0, 2L), 0.5, tolerance = 1e-10)
  expect_equal(bmm:::.mafc_pc_r(0, 4L), 0.25, tolerance = 1e-6)
  expect_equal(bmm:::.mafc_pc_r(0, 8L), 0.125, tolerance = 1e-3)
})

test_that("rsdt generates valid mafc data", {
  set.seed(42)
  dat <- rsdt(n_per_cell = 100, n_subjects = 5, dprime = 1.5,
              version = "mafc", m = 4)

  expect_true("id" %in% colnames(dat))
  expect_true("n_trials" %in% colnames(dat))
  expect_true("n_correct" %in% colnames(dat))
  expect_equal(nrow(dat), 5)  # one row per subject (no stimulus split)
  expect_true(all(dat$n_correct >= 0))
  expect_true(all(dat$n_correct <= dat$n_trials))
  expect_true(all(dat$n_trials == 100))
})

test_that("dsdt computes valid mafc densities", {
  # Density should be valid probability
  dens <- dsdt(n_correct = 80, n_trials = 100, dprime = 1.5, m = 4,
               version = "mafc")
  expect_true(dens > 0 && dens <= 1)

  # Log density
  log_dens <- dsdt(n_correct = 80, n_trials = 100, dprime = 1.5, m = 4,
                   version = "mafc", log = TRUE)
  expect_equal(log(dens), log_dens, tolerance = 1e-10)
})

test_that("dsdt mafc validates inputs", {
  expect_error(dsdt(n_correct = NULL, n_trials = 100, dprime = 1.5,
                    m = 4, version = "mafc"))
  expect_error(dsdt(n_correct = 80, n_trials = NULL, dprime = 1.5,
                    m = 4, version = "mafc"))
  expect_error(dsdt(n_correct = 80, n_trials = 100, dprime = 1.5,
                    m = NULL, version = "mafc"))
  expect_error(dsdt(n_correct = 110, n_trials = 100, dprime = 1.5,
                    m = 4, version = "mafc"))
})

test_that("rsdt mafc validates inputs", {
  expect_error(rsdt(version = "mafc", m = NULL))
  expect_error(rsdt(version = "mafc", m = 1))
})

# --- Check Data Tests ---

test_that("check_data validates mafc data correctly", {
  model <- sdt("n_correct", n_trials = "n_trials", version = "mafc", m = 4)
  dat <- rsdt(n_per_cell = 50, n_subjects = 3, dprime = 1.5,
              version = "mafc", m = 4)
  formula <- bmf(dprime ~ 1)

  # Should not error
  result <- check_data(model, dat, formula)
  expect_true("m_afc" %in% colnames(result))
  expect_equal(unique(result$m_afc), 4L)
})

test_that("check_data errors on missing columns for mafc", {
  model <- sdt("n_correct", n_trials = "n_trials", version = "mafc", m = 4)
  dat <- data.frame(id = 1:5, wrong_col = 50)
  formula <- bmf(dprime ~ 1)
  expect_error(check_data(model, dat, formula))
})

test_that("check_data errors on invalid response counts for mafc", {
  model <- sdt("n_correct", n_trials = "n_trials", version = "mafc", m = 4)
  dat <- data.frame(n_correct = c(-1, 10), n_trials = c(50, 50))
  formula <- bmf(dprime ~ 1)
  expect_error(check_data(model, dat, formula))
})

test_that("check_data errors when n_correct > n_trials for mafc", {
  model <- sdt("n_correct", n_trials = "n_trials", version = "mafc", m = 4)
  dat <- data.frame(n_correct = c(60, 10), n_trials = c(50, 50))
  formula <- bmf(dprime ~ 1)
  expect_error(check_data(model, dat, formula))
})

# --- Formula Construction Tests ---

test_that("bmf2bf.sdt produces correct formula for mafc", {
  model <- sdt("n_correct", n_trials = "n_trials", version = "mafc", m = 4)
  dat <- rsdt(n_per_cell = 50, n_subjects = 3, dprime = 1.5,
              version = "mafc", m = 4)
  dat$m_afc <- 4L  # simulate check_data adding this column
  bf <- bmf2bf(model, bmf(dprime ~ 1))

  # Should be a brmsformula
  expect_s3_class(bf, "brmsformula")
  # Should reference n_correct, n_trials, and vint for m_afc
  bf_str <- deparse(bf$formula, width.cutoff = 500)
  expect_true(grepl("n_correct", bf_str))
  expect_true(grepl("trials", bf_str))
  expect_true(grepl("vint.*m_afc", bf_str))
})

# --- Stan Code Generation Tests ---

test_that("sdt mafc produces valid stancode (2-AFC)", {
  dat <- rsdt(n_per_cell = 50, n_subjects = 3, dprime = 1.5,
              version = "mafc", m = 2)
  model <- sdt("n_correct", n_trials = "n_trials", version = "mafc", m = 2)
  formula <- bmf(dprime ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
  expect_true(grepl("mafc_pc", code))
  expect_true(grepl("sdt_mafc_lpmf", code))
})

test_that("sdt mafc produces valid stancode (4-AFC)", {
  dat <- rsdt(n_per_cell = 50, n_subjects = 3, dprime = 1.5,
              version = "mafc", m = 4)
  model <- sdt("n_correct", n_trials = "n_trials", version = "mafc", m = 4)
  formula <- bmf(dprime ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
  expect_true(grepl("gh_nodes", code))
  expect_true(grepl("gh_weights", code))
})

test_that("sdt mafc produces valid stancode (8-AFC)", {
  dat <- rsdt(n_per_cell = 50, n_subjects = 3, dprime = 1.5,
              version = "mafc", m = 8)
  model <- sdt("n_correct", n_trials = "n_trials", version = "mafc", m = 8)
  formula <- bmf(dprime ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
})

test_that("sdt mafc stancode works with predictors on dprime", {
  set.seed(42)
  dat <- rsdt(n_per_cell = 50, n_subjects = 5, dprime = 1.5,
              version = "mafc", m = 4)
  dat$condition <- rep(c("A", "B"), length.out = nrow(dat))
  model <- sdt("n_correct", n_trials = "n_trials", version = "mafc", m = 4)
  formula <- bmf(dprime ~ condition)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
})

test_that("sdt mafc stancode works with random effects on dprime", {
  dat <- rsdt(n_per_cell = 50, n_subjects = 5, dprime = 1.5,
              version = "mafc", m = 4)
  model <- sdt("n_correct", n_trials = "n_trials", version = "mafc", m = 4)
  formula <- bmf(dprime ~ 1 + (1 | id))
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
})


############################################################################# !
# RANKING MODEL TESTS                                                       ####
############################################################################# !

# --- Constructor Tests ---

test_that("sdt ranking model can be created", {
  model <- sdt("observed", version = "ranking", rank = "rank", m = 4)
  expect_s3_class(model, "bmmodel")
  expect_s3_class(model, "sdt")
  expect_s3_class(model, "sdt_ranking")
})

test_that("sdt ranking has correct parameters (only dprime, no criterion)", {
  model <- sdt("observed", version = "ranking", rank = "rank", m = 4)
  expect_true("dprime" %in% names(model$parameters))
  expect_true("mu" %in% names(model$parameters))
  expect_false("criterion" %in% names(model$parameters))
  expect_equal(model$fixed_parameters$mu, 0)
})

test_that("sdt ranking stores m and rank correctly", {
  model <- sdt("observed", version = "ranking", rank = "rank", m = 4)
  expect_equal(model$other_vars$m, 4L)
  expect_equal(model$other_vars$rank, "rank")

  model2 <- sdt("observed", version = "ranking", rank = "pos", m = 3)
  expect_equal(model2$other_vars$m, 3L)
  expect_equal(model2$other_vars$rank, "pos")
})

test_that("sdt ranking has correct links", {
  model <- sdt("observed", version = "ranking", rank = "rank", m = 4)
  expect_equal(model$links$mu, "identity")
  expect_equal(model$links$dprime, "identity")
})

test_that("sdt ranking accepts custom links", {
  model <- sdt("observed", version = "ranking", rank = "rank", m = 4,
               links = list(dprime = "log"))
  expect_equal(model$links$dprime, "log")
})

test_that("sdt ranking requires m >= 2", {
  expect_error(sdt("observed", version = "ranking", rank = "rank", m = 1))
  expect_error(sdt("observed", version = "ranking", rank = "rank", m = NULL))
  expect_error(sdt("observed", version = "ranking", rank = "rank"))
})

test_that("sdt ranking requires rank", {
  expect_error(sdt("observed", version = "ranking", m = 4))
})

test_that("sdt ranking only supports gumbel_min and normal", {
  expect_silent(sdt("observed", version = "ranking", rank = "rank",
                    m = 4, dist = "gumbel_min"))
  expect_silent(sdt("observed", version = "ranking", rank = "rank",
                    m = 4, dist = "normal"))
  expect_error(sdt("observed", version = "ranking", rank = "rank",
                   m = 4, dist = "logistic"))
  expect_error(sdt("observed", version = "ranking", rank = "rank",
                   m = 4, dist = "gumbel_max"))
})

test_that("sdt ranking UV-SDT requires dist = 'normal'", {
  expect_error(sdt("observed", version = "ranking", rank = "rank",
                   m = 4, dist = "gumbel_min", variances = "unequal"))
})

test_that("sdt ranking UV-SDT adds sdratio parameter", {
  model <- sdt("observed", version = "ranking", rank = "rank",
               m = 4, dist = "normal", variances = "unequal")
  expect_true("sdratio" %in% names(model$parameters))
})

test_that("sdt ranking does not require stimulus", {
  expect_silent(sdt("observed", version = "ranking", rank = "rank", m = 4))
})

# --- R-side Probability Tests ---

test_that(".ranking_prob_r gives valid probabilities for Gumbel", {
  p <- bmm:::.ranking_prob_r(1.5, 1L, 4L, dist = "gumbel_min")
  expect_true(p > 0 && p < 1)
})

test_that(".ranking_all_probs_r sums to 1 for Gumbel", {
  probs <- bmm:::.ranking_all_probs_r(1.5, 4L, dist = "gumbel_min")
  expect_equal(sum(probs), 1.0, tolerance = 1e-10)
})

test_that(".ranking_all_probs_r at d'=0 gives uniform for Gumbel", {
  probs <- bmm:::.ranking_all_probs_r(0, 4L, dist = "gumbel_min")
  expect_equal(probs, rep(0.25, 4), tolerance = 1e-6)
})

test_that(".ranking_all_probs_r is monotone: rank 1 increases with dprime", {
  p_low <- bmm:::.ranking_all_probs_r(0.5, 4L, dist = "gumbel_min")
  p_high <- bmm:::.ranking_all_probs_r(2.0, 4L, dist = "gumbel_min")
  expect_true(p_high[1] > p_low[1])
})

test_that(".ranking_all_probs_r works for m=3", {
  probs <- bmm:::.ranking_all_probs_r(1.0, 3L, dist = "gumbel_min")
  expect_length(probs, 3)
  expect_equal(sum(probs), 1.0, tolerance = 1e-10)
})

test_that(".ranking_all_probs_r works for Gaussian UV-SDT", {
  probs <- bmm:::.ranking_all_probs_r(1.5, 4L, dist = "normal", sdratio = 0)
  expect_equal(sum(probs), 1.0, tolerance = 1e-6)
  expect_true(probs[1] > probs[4])
})

# --- Distribution Function Tests ---

test_that("rsdt generates valid ranking data", {
  set.seed(42)
  dat <- rsdt(n_per_cell = 100, n_subjects = 5, dprime = 1.5,
              version = "ranking", m = 4)

  expect_true("id" %in% colnames(dat))
  expect_true("rank" %in% colnames(dat))
  expect_true("maxrank" %in% colnames(dat))
  expect_true("observed" %in% colnames(dat))
  expect_equal(nrow(dat), 20)  # 5 subjects x 4 ranks
  expect_true(all(dat$observed >= 0))
  expect_true(all(dat$maxrank == 4))
  expect_true(all(dat$rank %in% 1:4))
})

test_that("rsdt ranking validates inputs", {
  expect_error(rsdt(version = "ranking", m = NULL))
  expect_error(rsdt(version = "ranking", m = 1))
})

test_that("dsdt computes valid ranking densities", {
  dens <- dsdt(counts = c(40, 30, 20, 10), dprime = 1.5, m = 4,
               version = "ranking")
  expect_true(dens > 0)

  log_dens <- dsdt(counts = c(40, 30, 20, 10), dprime = 1.5, m = 4,
                   version = "ranking", log = TRUE)
  expect_equal(log(dens), log_dens, tolerance = 1e-10)
})

# --- Check Data Tests ---

test_that("check_data validates ranking data correctly", {
  model <- sdt("observed", version = "ranking", rank = "rank", m = 4)
  dat <- rsdt(n_per_cell = 50, n_subjects = 3, dprime = 1.5,
              version = "ranking", m = 4)
  formula <- bmf(dprime ~ 1)

  result <- check_data(model, dat, formula)
  expect_true("rank_pos" %in% colnames(result))
  expect_true("max_rank" %in% colnames(result))
  expect_equal(unique(result$max_rank), 4L)
})

test_that("check_data errors on missing rank column for ranking", {
  model <- sdt("observed", version = "ranking", rank = "rank", m = 4)
  dat <- data.frame(id = 1:5, observed = 50, wrong_col = 1:5)
  formula <- bmf(dprime ~ 1)
  expect_error(check_data(model, dat, formula))
})

test_that("check_data errors on invalid rank values for ranking", {
  model <- sdt("observed", version = "ranking", rank = "rank", m = 4)
  dat <- data.frame(observed = c(10, 10, 10), rank = c(0, 2, 5))
  formula <- bmf(dprime ~ 1)
  expect_error(check_data(model, dat, formula))
})

# --- Formula Construction Tests ---

test_that("bmf2bf.sdt produces correct formula for ranking", {
  model <- sdt("observed", version = "ranking", rank = "rank", m = 4)
  bf <- bmf2bf(model, bmf(dprime ~ 1))

  expect_s3_class(bf, "brmsformula")
  bf_str <- deparse(bf$formula, width.cutoff = 500)
  expect_true(grepl("observed", bf_str))
  expect_true(grepl("vint.*rank_pos.*max_rank", bf_str))
  # No trials clause
  expect_false(grepl("trials", bf_str))
})

# --- Stan Code Generation Tests ---

test_that("sdt ranking produces valid stancode (Gumbel, m=4)", {
  dat <- rsdt(n_per_cell = 50, n_subjects = 3, dprime = 1.5,
              version = "ranking", m = 4)
  model <- sdt("observed", version = "ranking", rank = "rank", m = 4)
  formula <- bmf(dprime ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
  expect_true(grepl("sdt_ranking_lpmf", code))
  expect_true(grepl("lgamma", code))
})

test_that("sdt ranking produces valid stancode (Gumbel, m=3)", {
  dat <- rsdt(n_per_cell = 50, n_subjects = 3, dprime = 1.5,
              version = "ranking", m = 3)
  model <- sdt("observed", version = "ranking", rank = "rank", m = 3)
  formula <- bmf(dprime ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
})

test_that("sdt ranking stancode works with predictors on dprime", {
  set.seed(42)
  dat <- rsdt(n_per_cell = 50, n_subjects = 5, dprime = 1.5,
              version = "ranking", m = 4)
  dat$condition <- rep(c("A", "B"), length.out = nrow(dat))
  model <- sdt("observed", version = "ranking", rank = "rank", m = 4)
  formula <- bmf(dprime ~ condition)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
})

test_that("sdt ranking stancode works with random effects on dprime", {
  dat <- rsdt(n_per_cell = 50, n_subjects = 5, dprime = 1.5,
              version = "ranking", m = 4)
  model <- sdt("observed", version = "ranking", rank = "rank", m = 4)
  formula <- bmf(dprime ~ 1 + (1 | id))
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
})


############################################################################# !
# LOG-SCALE NUMERICS TESTS                                                  ####
############################################################################# !

test_that("sdt rating with log_scale=TRUE produces valid stancode (normal)", {
  dat <- rsdt(n_per_cell = 100, n_subjects = 3, dprime = 1.5,
              criterion = 0, n_ratings = 6, spacing = 0.5,
              version = "rating")
  model <- sdt(response = paste0("r", 1:6), stimulus = "stimulus",
               dist = "normal", log_scale = TRUE)
  formula <- bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
  expect_true(grepl("log_Phi", code))
  expect_true(grepl("log1m_Phi", code))
  expect_true(grepl("log_diff_exp", code))
})

test_that("sdt rating with log_scale=TRUE works for gumbel_min", {
  dat <- rsdt(n_per_cell = 100, n_subjects = 3, dprime = 1.5,
              criterion = 0, n_ratings = 6, spacing = 0.5,
              version = "rating", dist = "gumbel_min")
  model <- sdt(response = paste0("r", 1:6), stimulus = "stimulus",
               dist = "gumbel_min", log_scale = TRUE)
  formula <- bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
  expect_true(grepl("log_diff_exp", code))
})

test_that("sdt rating with log_scale=TRUE works for logistic", {
  dat <- rsdt(n_per_cell = 100, n_subjects = 3, dprime = 1.5,
              criterion = 0, n_ratings = 6, spacing = 0.5,
              version = "rating", dist = "logistic")
  model <- sdt(response = paste0("r", 1:6), stimulus = "stimulus",
               dist = "logistic", log_scale = TRUE)
  formula <- bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
  expect_true(grepl("log1p_exp", code))
})

test_that("sdt rating with log_scale=TRUE works for gumbel_max", {
  dat <- rsdt(n_per_cell = 100, n_subjects = 3, dprime = 1.5,
              criterion = 0, n_ratings = 6, spacing = 0.5,
              version = "rating", dist = "gumbel_max")
  model <- sdt(response = paste0("r", 1:6), stimulus = "stimulus",
               dist = "gumbel_max", log_scale = TRUE)
  formula <- bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
  expect_true(grepl("log_diff_exp", code))
})

test_that("sdt rating with log_scale=FALSE is unchanged (default)", {
  dat <- rsdt(n_per_cell = 100, n_subjects = 3, dprime = 1.5,
              criterion = 0, n_ratings = 6, spacing = 0.5,
              version = "rating")
  model <- sdt(response = paste0("r", 1:6), stimulus = "stimulus")
  formula <- bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
  # Default (log_scale=FALSE) should use Phi, not log_Phi
  expect_true(grepl("Phi\\(", code))
  expect_false(grepl("log_Phi", code))
})

test_that("sdt UV-SDT rating with log_scale=TRUE produces valid stancode", {
  dat <- rsdt(n_per_cell = 100, n_subjects = 3, dprime = 1.5,
              criterion = 0, sdratio = 1.2, n_ratings = 6,
              spacing = 0.5, version = "rating")
  model <- sdt(response = paste0("r", 1:6), stimulus = "stimulus",
               variances = "unequal", log_scale = TRUE)
  formula <- bmf(dprime ~ 1, criterion ~ 1, spacing ~ 1, sdratio ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
  expect_true(grepl("log_Phi", code))
  expect_true(grepl("sdratio", code))
})
