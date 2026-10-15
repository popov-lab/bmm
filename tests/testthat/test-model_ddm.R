# Test DDM model specification and integration

test_that("ddm model can be created", {
  expect_silent(ddm("rt", "response"))
})

test_that("ddm model has correct class structure", {
  model <- ddm("rt", "response")
  expect_s3_class(model, "bmmodel")
  expect_s3_class(model, "ddm")
})

test_that("ddm model parameters are correctly defined", {
  model <- ddm("rt", "response")
  expect_true("drift" %in% names(model$parameters))
  expect_true("bound" %in% names(model$parameters))
  expect_true("ndt" %in% names(model$parameters))
  expect_true("zr" %in% names(model$parameters))
  expect_equal(model$fixed_parameters$zr, 0)
})

test_that("ddm model has correct link functions", {
  model <- ddm("rt", "response")
  expect_equal(model$links$drift, "identity")
  expect_equal(model$links$bound, "log")
  expect_equal(model$links$ndt, "log")
  expect_equal(model$links$zr, "logit")
})

test_that("ddm model accepts custom links", {
  custom_links <- list(drift = "log")
  model <- ddm("rt", "response", links = custom_links)
  expect_equal(model$links$drift, "log")
})

test_that("ddm fixed parameters are updated by formula", {
  model <- ddm("rt", "response")
  expect_equal(model$fixed_parameters$zr, 0)

  formula_free_zr <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1, zr ~ 1)
  updated <- update_model_fixed_parameters(model, formula_free_zr)
  expect_false("zr" %in% names(updated$fixed_parameters))
})

test_that("ddm check_data validates rt variable", {
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
  model <- ddm("rt", "response")
  
  # Invalid response codes
  invalid_data <- data.frame(rt = c(0.5, 0.6), response = c(2, 3))
  expect_error(
    check_data(model, invalid_data, bmf(drift ~ 1, bound ~ 1, ndt ~ 1)),
    "only contain values of 0 and 1"
  )
})

test_that("ddm check_data handles missing values", {
  model <- ddm("rt", "response")
  
  # Missing values produce a warning and are removed, not an error
  data_with_na <- data.frame(rt = c(0.5, NA, 0.7), response = c(1, 1, 0))
  expect_warning(
    check_data(model, data_with_na, bmf(drift ~ 1, bound ~ 1, ndt ~ 1)),
    "were NA"
  )
})

test_that("ddm check_data drops rows with an NA response and says so", {
  model <- ddm("rt", "response")
  dat <- data.frame(rt = c(0.5, 0.6, 0.7), response = c(1, NA, 0))
  expect_warning(out <- check_data(model, dat, bmf(drift ~ 1, bound ~ 1, ndt ~ 1)),
                 "Some values in response were NA")
  expect_equal(nrow(out), 2)
})

test_that("ddm check_data recodes logical and upper/lower responses to 1 and 0", {
  model <- ddm("rt", "response")
  formula <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1)
  logical <- data.frame(rt = c(0.5, 0.6), response = c(TRUE, FALSE))
  expect_warning(out <- check_data(model, logical, formula), "response variable is boolean")
  expect_equal(out$response, c(1, 0))
  # case does not matter
  character <- data.frame(rt = c(0.5, 0.6, 0.7), response = c("Upper", "lower", "UPPER"))
  expect_warning(out <- check_data(model, character, formula), "is a character variable")
  expect_equal(out$response, c(1, 0, 1))
})

test_that("ddm check_data refuses RTs and responses it cannot read", {
  model <- ddm("rt", "response")
  formula <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1)
  expect_error(
    check_data(model, data.frame(rt = c("0.5", "0.6"), response = c(1, 0)), formula),
    "The RT variable 'rt' needs to be of type double or integer"
  )
  expect_error(
    check_data(model, data.frame(rt = c(0.5, 0.6), response = c("upper", "up")), formula),
    "contains invalid character values"
  )
  expect_error(
    check_data(model, data.frame(rt = c(0.5, 0.6), response = as.complex(c(1, 0))), formula),
    "is of type 'complex'"
  )
})

test_that("ddm works with mock backend - fixed zr", {
  skip_on_cran()
  sim_data <- rddm(50, drift = 2, bound = 1.5, ndt = 0.3)
  model <- ddm("rt", "response")
  formula <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1)

  expect_silent(bmm(formula, sim_data, model, backend = "mock", mock = 1, rename = FALSE))
})

test_that("ddm works with mock backend - free zr", {
  skip_on_cran()
  sim_data <- rddm(50, drift = 2, bound = 1.5, ndt = 0.3, zr = 0.5)
  model <- ddm("rt", "response")
  formula <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1, zr ~ 1)

  expect_silent(bmm(formula, sim_data, model, backend = "mock", mock = 1, rename = FALSE))
})

test_that("ddm formula conversion works correctly", {
  skip_on_cran()
  
  sim_data <- rddm(50, drift = 2, bound = 1.5, ndt = 0.3)
  model <- ddm("rt", "response")
  formula <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1)
  
  fit <- bmm(formula, sim_data, model, backend = "mock", mock = 1, rename = FALSE)
  
  # Check that formula was converted properly
  expect_s3_class(fit$formula, "brmsformula")
})

test_that("ddm with condition effects works", {
  skip_on_cran()
  
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
  # under load_all() brms would also find the function by name on the search
  # path, so the stored function is what shows that restructure() added it
  expect_true(is.function(restructure(fit)$formula$family$posterior_epred))
  epred <- brms::posterior_epred(fit, ndraws = 20)
  expect_equal(dim(epred), c(20L, nrow(fit$data)))
  expect_true(all(is.finite(epred)))
})

test_that("log_lik() of a ddm fit is the density of each observation under each draw", {
  # two draws (rows) of two observations (columns), every value different
  dpars <- list(
    drift = matrix(c(1.5, 2, -0.5, 1), 2, 2), bound = matrix(c(1, 1.2, 0.9, 1.1), 2, 2),
    ndt = matrix(c(0.2, 0.25, 0.22, 0.18), 2, 2), zr = matrix(c(0.4, 0.6, 0.5, 0.55), 2, 2)
  )
  prep <- epred_prep(dpars, list(Y = c(0.6, 0.9), dec = c(1, 0)))
  for (i in 1:2) {
    expect_equal(log_lik_ddm(i, prep),
                 dddm(prep$data$Y[i], prep$data$dec[i], dpars$drift[, i], dpars$bound[, i],
                      dpars$ndt[, i], zr = dpars$zr[, i], log = TRUE))
  }
})
