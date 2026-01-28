# Test DDM model specification and integration

test_that("ddm model can be created with all versions", {
  # All DDM versions require cmdstanr (use Wiener likelihood not in rstan)
  skip_if_not(
    requireNamespace("cmdstanr", quietly = TRUE),
    "cmdstanr is required for DDM models"
  )
  
  expect_silent(ddm("rt", "response", version = "3par"))
  expect_silent(ddm("rt", "response", version = "4par"))
})

test_that("ddm model has correct class structure", {
  skip_if_not(
    requireNamespace("cmdstanr", quietly = TRUE),
    "cmdstanr is required for DDM models"
  )
  
  model <- ddm("rt", "response", version = "3par")
  expect_s3_class(model, "bmmodel")
  expect_s3_class(model, "ddm")
  expect_s3_class(model, "ddm_3par")
})

test_that("ddm model parameters are correctly defined for 3par version", {
  skip_if_not(
    requireNamespace("cmdstanr", quietly = TRUE),
    "cmdstanr is required for DDM models"
  )
  
  model <- ddm("rt", "response", version = "3par")
  expect_true("drift" %in% names(model$parameters))
  expect_true("bound" %in% names(model$parameters))
  expect_true("ndt" %in% names(model$parameters))
  expect_equal(model$fixed_parameters$zr, 0.5)
})

test_that("ddm model parameters are correctly defined for 4par version", {
  skip_if_not(
    requireNamespace("cmdstanr", quietly = TRUE),
    "cmdstanr is required for DDM models"
  )
  
  model <- ddm("rt", "response", version = "4par")
  expect_true("drift" %in% names(model$parameters))
  expect_true("bound" %in% names(model$parameters))
  expect_true("ndt" %in% names(model$parameters))
  expect_true("zr" %in% names(model$parameters))
  expect_length(model$fixed_parameters, 0)
})

test_that("ddm model has correct link functions", {
  skip_if_not(
    requireNamespace("cmdstanr", quietly = TRUE),
    "cmdstanr is required for DDM models"
  )
  
  model <- ddm("rt", "response", version = "4par")
  expect_equal(model$links$drift, "identity")
  expect_equal(model$links$bound, "log")
  expect_equal(model$links$ndt, "log")
  expect_equal(model$links$zr, "logit")
})

test_that("ddm model accepts custom links", {
  skip_if_not(
    requireNamespace("cmdstanr", quietly = TRUE),
    "cmdstanr is required for DDM models"
  )
  
  custom_links <- list(drift = "log")
  model <- ddm("rt", "response", version = "3par", links = custom_links)
  expect_equal(model$links$drift, "log")
})

test_that("ddm legacy version names work with deprecation warning", {
  skip_if_not(
    requireNamespace("cmdstanr", quietly = TRUE),
    "cmdstanr is required for DDM models"
  )
  
  expect_warning(
    ddm("rt", "response", version = "three_par"),
    "outdated version label"
  )
  expect_warning(
    ddm("rt", "response", version = "four_par"),
    "outdated version label"
  )
})

test_that("ddm check_data validates rt variable", {
  skip_if_not(
    requireNamespace("cmdstanr", quietly = TRUE),
    "cmdstanr is required for DDM models"
  )
  
  model <- ddm("rt", "response", version = "3par")
  
  # Valid data
  valid_data <- data.frame(rt = c(0.5, 0.6, 0.7), response = c(1, 1, 0))
  expect_silent(check_data(model, valid_data, bmf(drift ~ 1, bound ~ 1, ndt ~ 1)))
  
  # Negative RTs should error
  invalid_data <- data.frame(rt = c(-0.5, 0.6), response = c(1, 0))
  expect_error(
    check_data(model, invalid_data, bmf(drift ~ 1, bound ~ 1, ndt ~ 1)),
    "lower than zero"
  )
})

test_that("ddm check_data validates response variable", {
  skip_if_not(
    requireNamespace("cmdstanr", quietly = TRUE),
    "cmdstanr is required for DDM models"
  )
  
  model <- ddm("rt", "response", version = "3par")
  
  # Invalid response codes
  invalid_data <- data.frame(rt = c(0.5, 0.6), response = c(2, 3))
  expect_error(
    check_data(model, invalid_data, bmf(drift ~ 1, bound ~ 1, ndt ~ 1)),
    "only contain values of zero and one"
  )
})

test_that("ddm check_data handles missing values", {
  skip_if_not(
    requireNamespace("cmdstanr", quietly = TRUE),
    "cmdstanr is required for DDM models"
  )
  
  model <- ddm("rt", "response", version = "3par")
  
  # Missing values produce a warning and are removed, not an error
  data_with_na <- data.frame(rt = c(0.5, NA, 0.7), response = c(1, 1, 0))
  expect_warning(
    check_data(model, data_with_na, bmf(drift ~ 1, bound ~ 1, ndt ~ 1)),
    "were NA"
  )
})

test_that("ddm works with mock backend - 3par version", {
  skip_on_cran()
  skip_if_not(
    requireNamespace("cmdstanr", quietly = TRUE),
    "cmdstanr is required for DDM models"
  )
  
  sim_data <- rddm(50, drift = 2, bound = 1.5, ndt = 0.3)
  model <- ddm("rt", "response", version = "3par")
  formula <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1)
  
  expect_silent(bmm(formula, sim_data, model, backend = "mock", mock = 1, rename = FALSE))
})

test_that("ddm works with mock backend - 4par version", {
  skip_on_cran()
  skip_if_not(
    requireNamespace("cmdstanr", quietly = TRUE),
    "cmdstanr is required for DDM models"
  )
  
  sim_data <- rddm(50, drift = 2, bound = 1.5, ndt = 0.3, zr = 0.5)
  model <- ddm("rt", "response", version = "4par")
  formula <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1, zr ~ 1)
  
  expect_silent(bmm(formula, sim_data, model, backend = "mock", mock = 1, rename = FALSE))
})

test_that("ddm formula conversion works correctly", {
  skip_on_cran()
  skip_if_not(
    requireNamespace("cmdstanr", quietly = TRUE),
    "cmdstanr is required for DDM models"
  )
  
  sim_data <- rddm(50, drift = 2, bound = 1.5, ndt = 0.3)
  model <- ddm("rt", "response", version = "3par")
  formula <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1)
  
  fit <- bmm(formula, sim_data, model, backend = "mock", mock = 1, rename = FALSE)
  
  # Check that formula was converted properly
  expect_s3_class(fit$formula, "brmsformula")
})

test_that("ddm with condition effects works", {
  skip_on_cran()
  skip_if_not(
    requireNamespace("cmdstanr", quietly = TRUE),
    "cmdstanr is required for DDM models"
  )
  
  # Simulate data with condition effects
  n_per_cond <- 25
  data_a <- rddm(n_per_cond, drift = 2.5, bound = 1.5, ndt = 0.3)
  data_a$condition <- "A"
  data_b <- rddm(n_per_cond, drift = 1.5, bound = 1.5, ndt = 0.3)
  data_b$condition <- "B"
  sim_data <- rbind(data_a, data_b)
  sim_data$condition <- factor(sim_data$condition)
  
  model <- ddm("rt", "response", version = "3par")
  formula <- bmf(drift ~ 0 + condition, bound ~ 1, ndt ~ 1)
  
  expect_silent(bmm(formula, sim_data, model, backend = "mock", mock = 1, rename = FALSE))
})

test_that("ddm with hierarchical structure works", {
  skip_on_cran()
  skip_if_not(
    requireNamespace("cmdstanr", quietly = TRUE),
    "cmdstanr is required for DDM models"
  )
  
  # Simulate hierarchical data
  n_subjects <- 3
  n_per_subject <- 20
  
  data_list <- lapply(1:n_subjects, function(i) {
    d <- rddm(n_per_subject, drift = rnorm(1, 2, 0.3), bound = 1.5, ndt = 0.3)
    d$id <- paste0("S", i)
    d
  })
  sim_data <- do.call(rbind, data_list)
  sim_data$id <- factor(sim_data$id)
  
  model <- ddm("rt", "response", version = "3par")
  formula <- bmf(drift ~ 1 + (1 | id), bound ~ 1, ndt ~ 1)
  
  expect_silent(bmm(formula, sim_data, model, backend = "mock", mock = 1, rename = FALSE))
})

test_that("ddm allows missing parameters with message", {
  skip_on_cran()
  skip_if_not(
    requireNamespace("cmdstanr", quietly = TRUE),
    "cmdstanr is required for DDM models"
  )
  
  sim_data <- rddm(50, drift = 2, bound = 1.5, ndt = 0.3)
  model <- ddm("rt", "response", version = "3par")
  
  # Missing ndt parameter should work with message (not error)
  formula_incomplete <- bmf(drift ~ 1, bound ~ 1)
  expect_message(
    bmm(formula_incomplete, sim_data, model, backend = "mock", mock = 1, rename = FALSE),
    "No formula for parameter ndt"
  )
})

test_that("ddm default priors are correctly set", {
  skip_on_cran()
  skip_if_not(
    requireNamespace("cmdstanr", quietly = TRUE),
    "cmdstanr is required for DDM models"
  )
  
  sim_data <- rddm(50, drift = 2, bound = 1.5, ndt = 0.3)
  model <- ddm("rt", "response", version = "3par")
  formula <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1)
  
  fit <- bmm(formula, sim_data, model, backend = "mock", mock = 1, rename = FALSE)
  prior_summary <- brms::prior_summary(fit)
  
  # Check that priors are set for main parameters (in dpar column, not nlpar)
  expect_true(any(grepl("drift", prior_summary$dpar)))
  expect_true(any(grepl("bound", prior_summary$dpar)))
  expect_true(any(grepl("ndt", prior_summary$dpar)))
})

test_that("ddm stanvars are correctly added", {
  skip_on_cran()
  skip_if_not(
    requireNamespace("cmdstanr", quietly = TRUE),
    "cmdstanr is required for DDM models"
  )
  
  sim_data <- rddm(50, drift = 2, bound = 1.5, ndt = 0.3)
  model <- ddm("rt", "response", version = "3par")
  formula <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1)
  
  fit <- bmm(formula, sim_data, model, backend = "mock", mock = 1, rename = FALSE)
  
  # Check that custom Stan functions were added
  expect_true(!is.null(fit$stanvars))
})
