save_pars <- brms::save_pars

test_that("update.bmmfit works", {
  skip_if_not(interactive())
  fit1 <- restructure(readRDS(test_path("assets/bmmfit_example1.rds")))
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
