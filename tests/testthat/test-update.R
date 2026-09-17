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
  sum(grepl("^constant\\(", as.data.frame(fit$prior)$prior))
}

test_that("update.bmmfit frees parameters that the new formula predicts", {
  skip_on_cran()
  fit1 <- sdm_fixture()

  # mu is fixed to 0 by default in sdm; predicting it must drop the constant()
  # prior the original fit stored for it
  up <- update_mock(fit1, formula. = bmf(c ~ 1, kappa ~ 1, mu ~ 1))
  expect_length(up$bmm$model$fixed_parameters, 0)
  expect_equal(constant_priors(up), 0)
  expect_match(brms::stancode(up), "real Intercept;", fixed = TRUE)

  # freeing a parameter and replacing the data at once
  new_data <- data.frame(
    dev_rad = rsdm(80, c = 4, kappa = 3),
    set_size = rep(1:2, each = 40)
  )
  up <- update_mock(fit1, formula. = bmf(c ~ 1, kappa ~ 1, mu ~ 1), newdata = new_data)
  expect_length(up$bmm$model$fixed_parameters, 0)
  expect_equal(constant_priors(up), 0)
  expect_equal(nrow(up$data), 80)

  # a parameter left alone stays fixed
  expect_equal(constant_priors(update_mock(fit1)), 1)
})

test_that("update.bmmfit fixes parameters that the new formula sets to a constant", {
  skip_on_cran()
  fit1 <- sdm_fixture()

  # the reverse direction: the old fit's free prior on kappa must not survive and
  # override the constant() that the new formula asks for
  up <- update_mock(fit1, formula. = bmf(c ~ 0 + set_size, kappa = 5))
  expect_equal(up$bmm$model$fixed_parameters$kappa, 5)
  expect_true(any(grepl("^constant\\(5\\)", as.data.frame(up$prior)$prior)))
  expect_match(brms::stancode(up), "Intercept_kappa = 5;", fixed = TRUE)
  expect_false(grepl("student_t_lpdf(Intercept_kappa", brms::stancode(up), fixed = TRUE))

  # changing the value of an already fixed parameter is the same defect: the old
  # constant() row would otherwise survive and pin mu at the original 0
  up <- update_mock(fit1, formula. = bmf(c ~ 0 + set_size, mu = 0.5))
  expect_equal(up$bmm$model$fixed_parameters$mu, 0.5)
  expect_match(brms::stancode(up), "Intercept = 0.5;", fixed = TRUE)
})

# read the emitted program rather than bmm's own stanvars, and assert on the
# call: the sliced function is declared in the functions block either way
sdm_likelihood_is_sliced <- function(fit) {
  grepl("target += sdm_simple_run_ldenom_slice", brms::stancode(fit), fixed = TRUE)
}

test_that("update.bmmfit configures the likelihood for the effective threading spec", {
  skip_on_cran()
  fit1 <- sdm_fixture()

  # brms::update.brmsfit falls back to the fit's own threading spec when threads
  # is not passed, so the threaded chunk must be emitted without it too
  threaded_fit <- fit1
  threaded_fit$threads <- brms::threading(2)
  expect_true(sdm_likelihood_is_sliced(update_mock(threaded_fit)))

  # brms reads an explicit threads = NULL as "turn threading off", which it
  # distinguishes from an absent argument by name presence, not by NULL-ness
  expect_false(sdm_likelihood_is_sliced(update_mock(threaded_fit, threads = NULL)))

  expect_false(sdm_likelihood_is_sliced(update_mock(fit1)))
  expect_true(sdm_likelihood_is_sliced(update_mock(fit1, threads = brms::threading(2))))

  # the spec of the fit also has to win over a global brms.threads, which brms
  # itself ignores once it has fallen back to the fit's own spec
  withr::with_options(list(brms.threads = brms::threading(2)), {
    expect_false(sdm_likelihood_is_sliced(update_mock(fit1)))

    # a fit saved before brmsfit carried a $threads field
    no_threads <- fit1
    no_threads$threads <- NULL
    expect_false(sdm_likelihood_is_sliced(update_mock(no_threads)))
  })
})
