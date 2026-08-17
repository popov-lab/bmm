save_pars <- brms::save_pars

test_that("update.bmmfit works", {
  skip_on_cran()
  path <- test_path("assets/bmmfit_example1.rds")
  skip_if_not(file.exists(path), "SDM fixture not available (excluded by .Rbuildignore)")
  fit1 <- restructure(readRDS(path))
  data <- fit1$data

  # formula is replaced
  up <- suppressMessages(update(fit1, formula. = bmf(c ~ 1, kappa ~ 1), testmode = TRUE))
  expect_true(is(up, "bmmfit"))
  expect_equal(up$bmm$user_formula$c, c ~ 1, ignore_attr = TRUE)

  # data is replaced, old formula is kept
  new_data <- data
  new_data$dev_rad <- rnorm(nrow(new_data), 0, 0.5)
  up <- suppressMessages(
    update(fit1, newdata = new_data, save_pars = save_pars(group = FALSE), testmode = TRUE)
  )
  expect_true(is(up, "bmmfit"))
  expect_equal(attr(up$data, "data_name"), "new_data")
  expect_equal(up$bmm$user_formula$c, c ~ 0 + set_size, ignore_formula_env = T, ignore_attr = TRUE)

  # prior is replaced
  up <- suppressMessages(
    update(
      fit1,
      formula. = bmf(c ~ 1, kappa ~ 1), testmode = TRUE,
      prior = brms::set_prior("normal(0,0.1)", class = "Intercept", dpar = "kappa")
    )
  )
  expect_true(is(up, "bmmfit"))

  # refuse to change model
  expect_error(
    update(fit1, model = mixture2p(resp_error = "dev_rad")),
    "You cannot update with a different model"
  )

  up <- suppressMessages(update(fit1, save_pars = save_pars(group = FALSE), testmode = TRUE))
  expect_true(is(up, "bmmfit"))
  up <- suppressMessages(update(fit1, save_pars = save_pars(latent = FALSE), testmode = TRUE))
  expect_true(is(up, "bmmfit"))
  expect_error(update(fit1, data = data), "use argument 'newdata'")
})

sdm_fixture <- function() {
  path <- test_path("assets/bmmfit_example1.rds")
  skip_if_not(file.exists(path), "SDM fixture not available (excluded by .Rbuildignore)")
  restructure(readRDS(path))
}

update_mock <- function(object, ...) {
  suppressMessages(update(object, ..., backend = "mock", mock_fit = 1, rename = FALSE))
}

constant_priors <- function(fit) {
  sum(grepl("constant", as.data.frame(fit$prior)$prior))
}

test_that("update.bmmfit frees parameters that the new formula predicts", {
  skip_on_cran()
  fit1 <- sdm_fixture()

  # mu is fixed to 0 by default in sdm; predicting it must drop the constant()
  # prior the original fit stored for it
  up <- update_mock(fit1, formula. = bmf(c ~ 1, kappa ~ 1, mu ~ 1))
  expect_length(up$bmm$model$fixed_parameters, 0)
  expect_equal(constant_priors(up), 0)

  # a parameter left alone stays fixed
  expect_equal(constant_priors(update_mock(fit1)), 1)
})

sdm_likelihood_is_sliced <- function(fit) {
  blocks <- vapply(fit$stanvars, function(x) x$block %||% "", character(1))
  scode <- vapply(fit$stanvars[blocks == "likelihood"],
    function(x) paste(x$scode, collapse = "\n"), character(1)
  )
  any(grepl("ldenom_slice", scode, fixed = TRUE))
}

test_that("update.bmmfit configures the likelihood for the effective threading spec", {
  skip_on_cran()
  fit1 <- sdm_fixture()

  # brms::update.brmsfit falls back to the fit's own threading spec when threads
  # is not passed, so the threaded chunk must be emitted without it too
  threaded_fit <- fit1
  threaded_fit$threads <- brms::threading(2)
  expect_true(sdm_likelihood_is_sliced(update_mock(threaded_fit)))

  expect_false(sdm_likelihood_is_sliced(update_mock(fit1)))
  expect_true(sdm_likelihood_is_sliced(update_mock(fit1, threads = brms::threading(2))))
})
