test_that("add_links works", {
  x <- .model_mixture2p()
  links <- x$links
  x$links <- NULL
  x <- add_links(x)
  expect_equal(x$links, links)
})

# bmm 0.2.1 saved its fits as plain brmsfit objects; the model and the user's
# formula live only in the environment of the custom family
legacy_fit <- function() {
  fit <- load_fixture_fit("bmmfit_legacy_sdm_v0.2.1.rds")
  class(fit) <- c("bmmfit", class(fit))
  fit
}

test_that("a fit saved by bmm 0.2.1 is rebuilt as the current model", {
  fit <- restructure(legacy_fit())
  expect_identical(class(fit$bmm$model), c("bmmodel", "circular", "sdm", "sdm_simple"))
  expect_equal(fit$bmm$model$links, .model_sdm()$links)
  expect_equal(fit$bmm$user_formula, bmf(c ~ 0 + SetSize, kappa ~ 1), ignore_formula_env = TRUE)
  expect_identical(fit$version$bmm_restructure, utils::packageVersion("bmm"))
  expect_s3_class(summary(fit), "bmmsummary")
})

test_that("an old fit whose family kept no model is refused with the advice to refit", {
  fit <- legacy_fit()
  fit$family$env <- NULL
  expect_error(restructure(fit), "Please refit")
})

test_that("restructure() adds posterior_epred to a custom family saved without it", {
  fit <- load_fixture_fit("bmmfit_example1.rds")
  # as a fit saved before the family carried it
  fit$formula$family$posterior_epred <- NULL
  fit$family$posterior_epred <- NULL
  fit <- restructure(fit)
  expect_identical(fit$formula$family$posterior_epred, posterior_epred_sdm_simple)
  expect_identical(fit$family$posterior_epred, posterior_epred_sdm_simple)
})

test_that("restructure() leaves native brms families without posterior_epred", {
  fit <- load_fixture_fit("bmmfit_m3_ppcheck.rds")
  fit <- restructure(fit)
  expect_null(fit$formula$family$posterior_epred)
})
