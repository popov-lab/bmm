test_that("default_prior() works with brmsformula", {
  ff <- brms::bf(count ~ zAge + zBase * Trt + (1 | patient))
  prior <- default_prior(ff, data = brms::epilepsy, family = poisson())
  expect_equal(class(prior)[1], "brmsprior")
})

test_that("default_prior() works with formula", {
  ff <- count ~ zAge + zBase * Trt + (1 | patient)
  prior <- default_prior(ff, data = brms::epilepsy, family = poisson())
  expect_equal(class(prior)[1], "brmsprior")
})

test_that("default_prior() works with bmmformula", {
  ff <- bmmformula(kappa ~ 1, thetat ~ 1, thetant ~ 1)
  prior <- default_prior(ff, oberauer_lin_2017, mixture3p(
    resp_error = "dev_rad",
    nt_features = "col_nt",
    set_size = "set_size",
    regex = T
  ))
  expect_equal(class(prior)[1], "brmsprior")
})

test_that("combine prior returns a brmsprior object", {
  prior1 <- brms::prior(normal(0, 1), class = "sd", dpar = "c")
  prior2 <- brms::prior(normal(0, 1), class = "sd", dpar = "kappa")

  # combine the prior
  prior <- combine_prior(prior1, prior2)
  expect_equal(class(prior)[1], "brmsprior")
})

test_that("in combine prior, prior2 overwrites only shared components with prior1", {
  prior1 <- brms::prior(normal(0, 1), class = "sd", dpar = "c1") +
    brms::prior(normal(0, 1), class = "sd", dpar = "c2") +
    brms::prior(normal(0, 1), class = "sd", dpar = "c3")
  prior2 <- brms::prior(normal(0, 1), class = "sd", dpar = "kappa") +
    brms::prior(normal(0, 2), class = "sd", dpar = "c2")

  # combine the prior
  prior <- combine_prior(prior1, prior2)
  expect_equal(nrow(prior), 4)
  expect_equal(dplyr::filter(prior, dpar == "c1"), dplyr::filter(prior1, dpar == "c1"))
  expect_equal(dplyr::filter(prior, dpar == "c2"), dplyr::filter(prior2, dpar == "c2"))
  expect_equal(dplyr::filter(prior, dpar == "c3"), dplyr::filter(prior1, dpar == "c3"))
  expect_equal(dplyr::filter(prior, dpar == "kappa"), dplyr::filter(prior2, dpar == "kappa"))
})

test_that("combine_prior() handles empty priors", {
  empty <- brms::empty_prior()
  expect_equal(nrow(combine_prior(empty, empty)), 0)

  prior <- brms::prior(normal(0, 1), class = "sd", dpar = "c")
  expect_equal(nrow(combine_prior(empty, prior)), 1)
  expect_equal(nrow(combine_prior(prior, empty)), 1)
})

test_that("validate_default_priors() accepts a formula without parameter formulas", {
  model <- list(default_priors = list())
  expect_equal(validate_default_priors(model, brms::bf(y ~ 1)), list())
})

test_that("fixed_pars_priors errors clearly when a fixed parameter is absent from the formula", {
  model <- list(fixed_parameters = list(sdratio = 0))
  formula <- brms::bf(y ~ 1)
  expect_error(
    fixed_pars_priors(model, formula),
    "sdratio.*not part of the model formula"
  )
})

test_that("the fixed-parameter guard fires through the full default_prior pipeline", {
  # a fixed parameter that no configure_model wires into the formula survives
  # check_model/check_formula reconciliation and must be caught before brms
  model <- sdm(resp_error = "y")
  model$fixed_parameters$sdratio <- 0
  dat <- data.frame(y = rsdm(30, mu = 0, c = 3, kappa = 4))
  expect_error(
    default_prior(bmf(c ~ 1, kappa ~ 1), dat, model),
    "sdratio.*not part of the model formula"
  )
})

test_that("fixed_pars_priors builds constant priors for dpar and nlpar fixed parameters", {
  dpar_prior <- fixed_pars_priors(
    list(fixed_parameters = list(sigma = 1)),
    brms::bf(y ~ 1, sigma ~ 1)
  )
  expect_equal(dpar_prior$prior, "constant(1)")
  expect_equal(dpar_prior$dpar, "sigma")

  nlpar_prior <- fixed_pars_priors(
    list(fixed_parameters = list(b = 0)),
    brms::bf(y ~ a + b, a ~ 1, b ~ 1, nl = TRUE)
  )
  expect_equal(nlpar_prior$prior, "constant(0)")
  expect_equal(nlpar_prior$nlpar, "b")
})

test_that("no check for sort_data with default_priors function", {
  withr::local_options("bmm.sort_data" = "check")
  res <- capture_messages(default_prior(
    bmf(kappa ~ set_size, c ~ set_size),
    oberauer_lin_2017,
    sdm("dev_rad")
  ))
  expect_false(any(grepl("sort", res)))
})
