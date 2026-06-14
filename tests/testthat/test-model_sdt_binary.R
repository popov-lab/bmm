# Test Binary SDT model specification and integration

############################################################################# !
# MODEL CONSTRUCTOR TESTS                                                ####
############################################################################# !

test_that("sdt_binary model can be created with default arguments", {
  expect_silent(sdt_binary("n_old", "stimulus", "n_trials"))
})

test_that("sdt_binary model can be created with all distribution options", {
  expect_silent(sdt_binary("n_old", "stimulus", "n_trials", dist = "normal"))
  expect_silent(sdt_binary("n_old", "stimulus", "n_trials", dist = "logistic"))
  expect_silent(sdt_binary("n_old", "stimulus", "n_trials", dist = "gumbel_min"))
  expect_silent(sdt_binary("n_old", "stimulus", "n_trials", dist = "gumbel_max"))
})

test_that("sdt_binary model has correct class structure", {
  model <- sdt_binary("n_old", "stimulus", "n_trials")
  expect_s3_class(model, "bmmodel")
  expect_s3_class(model, "sdt")
  expect_s3_class(model, "sdt_binary")
})

test_that("sdt_binary model parameters are correctly defined", {
  model <- sdt_binary("n_old", "stimulus", "n_trials")
  expect_true("dprime" %in% names(model$parameters))
  expect_true("criterion" %in% names(model$parameters))
  expect_true("sdratio" %in% names(model$parameters))
  expect_equal(model$fixed_parameters$sdratio, 0)
})

test_that("sdt_binary model has correct default link functions", {
  model <- sdt_binary("n_old", "stimulus", "n_trials")
  expect_equal(model$links$dprime, "identity")
  expect_equal(model$links$criterion, "identity")
  expect_equal(model$links$sdratio, "identity")
})

test_that("sdt_binary model accepts custom links", {
  custom_links <- list(dprime = "log")
  model <- sdt_binary("n_old", "stimulus", "n_trials", links = custom_links)
  expect_equal(model$links$dprime, "log")
  expect_equal(model$links$criterion, "identity")
})

test_that("sdt_binary model stores distribution info correctly", {
  model <- sdt_binary("n_old", "stimulus", "n_trials", dist = "gumbel_min")
  expect_equal(model$other_vars$dist, "gumbel_min")
  expect_equal(model$other_vars$dist_int, 2L)

  model2 <- sdt_binary("n_old", "stimulus", "n_trials", dist = "logistic")
  expect_equal(model2$other_vars$dist, "logistic")
  expect_equal(model2$other_vars$dist_int, 4L)
})

test_that("sdt_binary model has default priors", {
  model <- sdt_binary("n_old", "stimulus", "n_trials")
  expect_true("dprime" %in% names(model$default_priors))
  expect_true("criterion" %in% names(model$default_priors))
  expect_true("main" %in% names(model$default_priors$dprime))
  expect_true("effects" %in% names(model$default_priors$dprime))
})

test_that("sdt_binary requires n_trials argument", {
  expect_error(sdt_binary("n_old", "stimulus"))
})

test_that("sdt_binary model has init_ranges for estimated SDT parameters", {
  model <- sdt_binary("n_old", "stimulus", "n_trials")
  expect_false(is.null(model$init_ranges))
  expect_true(all(c("dprime", "criterion", "sdratio") %in% names(model$init_ranges)))
  expect_equal(model$init_ranges$dprime, c(0.5, 1.5))
  expect_equal(model$init_ranges$criterion, c(-0.5, 0.5))
  expect_equal(model$init_ranges$sdratio, c(-0.3, 0.3))
})


############################################################################# !
# CHECK_DATA TESTS                                                       ####
############################################################################# !

test_that("sdt_binary check_data validates required columns", {
  model <- sdt_binary("n_old", "stimulus", "n_trials")
  formula <- bmf(dprime ~ 1, criterion ~ 1)

  valid_data <- data.frame(
    n_old = c(30, 40),
    stimulus = c(0L, 1L),
    n_trials = c(50, 50)
  )
  result <- check_data(model, valid_data, formula)
  expect_true("dist_type" %in% colnames(result))

  invalid_data <- data.frame(n_old = c(30, 40))
  expect_error(
    check_data(model, invalid_data, formula),
    "missing in the data"
  )
})

test_that("sdt_binary check_data validates stimulus coding", {
  model <- sdt_binary("n_old", "stimulus", "n_trials")
  formula <- bmf(dprime ~ 1, criterion ~ 1)

  invalid_data <- data.frame(
    n_old = c(30, 40),
    stimulus = c(1, 2),
    n_trials = c(50, 50)
  )
  expect_error(
    check_data(model, invalid_data, formula),
    "must be coded as 0"
  )
})

test_that("sdt_binary check_data validates response counts", {
  model <- sdt_binary("n_old", "stimulus", "n_trials")
  formula <- bmf(dprime ~ 1, criterion ~ 1)

  invalid_data <- data.frame(
    n_old = c(-5, 40),
    stimulus = c(0L, 1L),
    n_trials = c(50, 50)
  )
  expect_error(
    check_data(model, invalid_data, formula),
    "non-negative"
  )

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

test_that("sdt_binary check_data validates n_trials", {
  model <- sdt_binary("n_old", "stimulus", "n_trials")
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

test_that("sdt_binary check_data adds dist_type column", {
  model <- sdt_binary("n_old", "stimulus", "n_trials", dist = "gumbel_min")
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

test_that("rsdt_binary generates data with correct structure", {
  dat <- rsdt_binary(n_per_cell = 50, n_subjects = 5, dprime = 1.5, criterion = 0)
  expect_true(is.data.frame(dat))
  expect_equal(nrow(dat), 10)
  expect_true(all(c("id", "stimulus", "n_trials", "n_old") %in% colnames(dat)))
  expect_equal(unique(dat$stimulus), c(0L, 1L))
  expect_true(all(dat$n_old >= 0))
  expect_true(all(dat$n_old <= dat$n_trials))
})

test_that("rsdt_binary generates different data for different distributions", {
  set.seed(42)
  dat_norm <- rsdt_binary(n_per_cell = 1000, n_subjects = 1, dprime = 2,
                          criterion = 0)
  set.seed(42)
  dat_logis <- rsdt_binary(n_per_cell = 1000, n_subjects = 1, dprime = 2,
                           criterion = 0, dist = "logistic")
  expect_false(identical(dat_norm$n_old, dat_logis$n_old))
})

test_that("rsdt_binary with zero dprime produces equal hit and FA rates on average", {
  set.seed(123)
  dat <- rsdt_binary(n_per_cell = 10000, n_subjects = 1, dprime = 0, criterion = 0)
  noise_prop <- dat$n_old[dat$stimulus == 0] / dat$n_trials[dat$stimulus == 0]
  signal_prop <- dat$n_old[dat$stimulus == 1] / dat$n_trials[dat$stimulus == 1]
  expect_equal(noise_prop, 0.5, tolerance = 0.05)
  expect_equal(signal_prop, 0.5, tolerance = 0.05)
})

test_that("dsdt_binary returns valid binomial density", {
  d <- dsdt_binary(n_old = 50, n_trials = 100, stimulus = 1,
                   dprime = 0, criterion = 0)
  expect_true(d > 0)
  expect_true(d <= 1)
})

test_that("dsdt_binary log density matches log of density", {
  d <- dsdt_binary(n_old = 80, n_trials = 100, stimulus = 1,
                   dprime = 1.5, criterion = 0.2)
  ld <- dsdt_binary(n_old = 80, n_trials = 100, stimulus = 1,
                    dprime = 1.5, criterion = 0.2, log = TRUE)
  expect_equal(log(d), ld, tolerance = 1e-10)
})

test_that("dsdt_binary is vectorized over observations", {
  d <- dsdt_binary(n_old = c(30, 80), n_trials = c(100, 100),
                   stimulus = c(0, 1), dprime = 1.5, criterion = 0.2)
  expect_length(d, 2)
  expect_true(all(d > 0))
})

test_that("dsdt_binary works for all distributions", {
  dists <- c("normal", "logistic", "gumbel_min", "gumbel_max")
  for (di in dists) {
    d <- dsdt_binary(n_old = 70, n_trials = 100, stimulus = 1,
                     dprime = 1.5, criterion = 0, dist = di)
    expect_true(d > 0, info = paste("dist:", di))
  }
})

test_that("dsdt_binary validates input", {
  expect_error(dsdt_binary(n_old = -1, n_trials = 100, stimulus = 1,
                           dprime = 1, criterion = 0), "non-negative")
  expect_error(dsdt_binary(n_old = 50, n_trials = 100, stimulus = 2,
                           dprime = 1, criterion = 0), "0.*1")
})

test_that("dsdt_binary with sdratio produces different densities than EV-SDT", {
  d_ev <- dsdt_binary(n_old = 80, n_trials = 100, stimulus = 1,
                      dprime = 1.5, criterion = 0, sdratio = 1)
  d_uv <- dsdt_binary(n_old = 80, n_trials = 100, stimulus = 1,
                      dprime = 1.5, criterion = 0, sdratio = 1.3)
  expect_false(d_ev == d_uv)
})

test_that("dsdt_binary with sdratio = 1 matches equal-variance for noise trials", {
  d1 <- dsdt_binary(n_old = 70, n_trials = 100, stimulus = 0,
                    dprime = 1.5, criterion = 0.2, sdratio = 1)
  d2 <- dsdt_binary(n_old = 70, n_trials = 100, stimulus = 0,
                    dprime = 1.5, criterion = 0.2, sdratio = 1.5)
  expect_equal(d1, d2)
})

test_that("dsdt_binary is vectorized over sdratio", {
  d <- dsdt_binary(n_old = c(30, 80), n_trials = c(100, 100),
                   stimulus = c(0, 1), dprime = 1.5, criterion = 0.2,
                   sdratio = c(1.0, 1.3))
  expect_length(d, 2)
  expect_true(all(d > 0))
})


############################################################################# !
# FORMULA AND CONFIGURE_MODEL TESTS                                      ####
############################################################################# !

test_that("sdt_binary model produces valid stancode with brms", {
  dat <- data.frame(
    n_old = c(10, 40, 15, 35),
    stimulus = c(0L, 1L, 0L, 1L),
    n_trials = c(50, 50, 50, 50),
    id = c(1, 1, 2, 2)
  )
  model <- sdt_binary("n_old", "stimulus", "n_trials", dist = "normal")
  formula <- bmf(dprime ~ 1, criterion ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
  expect_true(grepl("sdt_binary", code))
  expect_true(grepl("sdt_cumprob", code))
})

test_that("sdt_binary model produces valid stancode for all distributions", {
  dat <- data.frame(
    n_old = c(10, 40),
    stimulus = c(0L, 1L),
    n_trials = c(50, 50)
  )
  for (d in c("normal", "logistic", "gumbel_min", "gumbel_max")) {
    model <- sdt_binary("n_old", "stimulus", "n_trials", dist = d)
    formula <- bmf(dprime ~ 1, criterion ~ 1)
    code <- stancode(formula, data = dat, model = model)
    expect_true(nchar(code) > 0, info = paste("Failed for dist:", d))
  }
})

test_that("sdt_binary model handles predictors in formula", {
  dat <- data.frame(
    n_old = c(10, 40, 15, 35),
    stimulus = c(0L, 1L, 0L, 1L),
    n_trials = c(50, 50, 50, 50),
    condition = c("A", "A", "B", "B")
  )
  model <- sdt_binary("n_old", "stimulus", "n_trials", dist = "normal")
  formula <- bmf(dprime ~ condition, criterion ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
})

test_that("sdt_binary model handles random effects in formula", {
  dat <- data.frame(
    n_old = c(10, 40, 15, 35),
    stimulus = c(0L, 1L, 0L, 1L),
    n_trials = c(50, 50, 50, 50),
    id = c(1, 1, 2, 2)
  )
  model <- sdt_binary("n_old", "stimulus", "n_trials", dist = "normal")
  formula <- bmf(dprime ~ 1 + (1 | id), criterion ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(nchar(code) > 0)
})

test_that("sdt_binary default_prior returns valid prior object", {
  dat <- data.frame(
    n_old = c(10, 40),
    stimulus = c(0L, 1L),
    n_trials = c(50, 50)
  )
  model <- sdt_binary("n_old", "stimulus", "n_trials", dist = "normal")
  formula <- bmf(dprime ~ 1, criterion ~ 1)
  prior <- default_prior(formula, data = dat, model = model)
  expect_s3_class(prior, "brmsprior")
  prior_strs <- prior$prior
  expect_true(any(grepl("normal", prior_strs)))
})
