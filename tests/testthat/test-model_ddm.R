# Test DDM model specification and integration

test_that("ddm model can be created", {
  # All DDM versions require cmdstanr (use Wiener likelihood not in rstan)
  skip_if_not(
    requireNamespace("cmdstanr", quietly = TRUE),
    "cmdstanr is required for DDM models"
  )

  expect_silent(ddm("rt", "response"))
})

test_that("ddm model has correct class structure", {
  skip_if_not(
    requireNamespace("cmdstanr", quietly = TRUE),
    "cmdstanr is required for DDM models"
  )

  model <- ddm("rt", "response")
  expect_s3_class(model, "bmmodel")
  expect_s3_class(model, "ddm")
})

test_that("ddm model parameters are correctly defined", {
  skip_if_not(
    requireNamespace("cmdstanr", quietly = TRUE),
    "cmdstanr is required for DDM models"
  )

  model <- ddm("rt", "response")
  expect_true("drift" %in% names(model$parameters))
  expect_true("bound" %in% names(model$parameters))
  expect_true("ndt" %in% names(model$parameters))
  expect_true("zr" %in% names(model$parameters))
  expect_equal(model$fixed_parameters$zr, 0)
})

test_that("ddm model has correct link functions", {
  skip_if_not(
    requireNamespace("cmdstanr", quietly = TRUE),
    "cmdstanr is required for DDM models"
  )
  
  model <- ddm("rt", "response")
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
  model <- ddm("rt", "response", links = custom_links)
  expect_equal(model$links$drift, "log")
})

test_that("ddm fixed parameters are updated by formula", {
  skip_if_not(
    requireNamespace("cmdstanr", quietly = TRUE),
    "cmdstanr is required for DDM models"
  )

  model <- ddm("rt", "response")
  expect_equal(model$fixed_parameters$zr, 0)

  formula_free_zr <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1, zr ~ 1)
  updated <- update_model_fixed_parameters(model, formula_free_zr)
  expect_false("zr" %in% names(updated$fixed_parameters))
})

test_that("ddm check_data validates rt variable", {
  skip_if_not(
    requireNamespace("cmdstanr", quietly = TRUE),
    "cmdstanr is required for DDM models"
  )
  
  model <- ddm("rt", "response")
  
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
  
  model <- ddm("rt", "response")
  
  # Invalid response codes
  invalid_data <- data.frame(rt = c(0.5, 0.6), response = c(2, 3))
  expect_error(
    check_data(model, invalid_data, bmf(drift ~ 1, bound ~ 1, ndt ~ 1)),
    "only contain values of 0 and 1"
  )
})

test_that("ddm check_data handles missing values", {
  skip_if_not(
    requireNamespace("cmdstanr", quietly = TRUE),
    "cmdstanr is required for DDM models"
  )
  
  model <- ddm("rt", "response")
  
  # Missing values produce a warning and are removed, not an error
  data_with_na <- data.frame(rt = c(0.5, NA, 0.7), response = c(1, 1, 0))
  expect_warning(
    check_data(model, data_with_na, bmf(drift ~ 1, bound ~ 1, ndt ~ 1)),
    "were NA"
  )
})

test_that("ddm works with mock backend - fixed zr", {
  skip_on_cran()
  skip_if_not(
    requireNamespace("cmdstanr", quietly = TRUE),
    "cmdstanr is required for DDM models"
  )

  sim_data <- rddm(50, drift = 2, bound = 1.5, ndt = 0.3)
  model <- ddm("rt", "response")
  formula <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1)

  expect_silent(bmm(formula, sim_data, model, backend = "mock", mock = 1, rename = FALSE))
})

test_that("ddm works with mock backend - free zr", {
  skip_on_cran()
  skip_if_not(
    requireNamespace("cmdstanr", quietly = TRUE),
    "cmdstanr is required for DDM models"
  )

  sim_data <- rddm(50, drift = 2, bound = 1.5, ndt = 0.3, zr = 0.5)
  model <- ddm("rt", "response")
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
  model <- ddm("rt", "response")
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
  
  model <- ddm("rt", "response")
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
  
  model <- ddm("rt", "response")
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
  model <- ddm("rt", "response")
  
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
  model <- ddm("rt", "response")
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
  model <- ddm("rt", "response")
  formula <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1)
  
  fit <- bmm(formula, sim_data, model, backend = "mock", mock = 1, rename = FALSE)
  
  # Check that custom Stan functions were added
  expect_true(!is.null(fit$stanvars))
})

# posterior_epred (#475) ------------------------------------------------------

ddm_epred_sets <- data.frame(
  drift = c(1.5, -0.8, 0, 3),
  bound = c(1.2, 1.8, 1, 2),
  ndt = c(0.3, 0.2, 0.25, 0.4),
  zr = c(0.5, 0.65, 0.3, 0.4)
)

test_that("posterior_epred_ddm() is the mean RT posterior_predict_ddm() simulates", {
  skip_on_cran()
  withr::local_seed(475)
  res <- epred_vs_predict(ddm_epred_sets, posterior_epred_ddm, posterior_predict_ddm)
  # the Monte-Carlo SE of each mean is below 0.5% of it at 20000 draws
  expect_lt(max(abs(res[, "rel_error"])), 0.02)
})

test_that("posterior_epred_ddm() returns one column per observation", {
  dpars <- lapply(ddm_epred_sets, matrix, nrow = 2)
  expect_epred_by_cell(posterior_epred_ddm, dpars)
})

test_that("a ddm fit stores posterior_epred_ddm() in its family", {
  skip_on_cran()
  sim_data <- rddm(50, drift = 2, bound = 1.5, ndt = 0.3)
  fit <- bmm(bmf(drift ~ 1, bound ~ 1, ndt ~ 1), sim_data, ddm("rt", "response"),
             backend = "mock", mock = 1, rename = FALSE)
  expect_identical(fit$formula$family$posterior_epred, posterior_epred_ddm)
})

test_that("posterior_epred() works on a ddm fit saved without the function", {
  skip_on_cran()
  fit <- load_fixture_fit("bmmfit_ddm_ppcheck.rds")
  epred <- brms::posterior_epred(fit, ndraws = 20)
  expect_equal(dim(epred), c(20L, nrow(fit$data)))
  expect_true(all(is.finite(epred)))
})
