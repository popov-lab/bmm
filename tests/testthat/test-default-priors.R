# mixture2p's parameters are distributional parameters of a custom family, so
# their priors are addressed by dpar rather than by nlpar, and an intercept
# prior carries the Intercept class instead of being stored as class "b" with
# coef "Intercept". Look priors up by parameter rather than by row order.
prior_for <- function(pr, dpar, class, coef = "") {
  pr$prior[pr$dpar == dpar & pr$class == class & pr$coef == coef]
}

test_that("default priors are set correctly with fixed effects only", {
  data <- oberauer_lin_2017
  model <- mixture2p("dev_rad")

  # Intercept only
  formula <- bmf(kappa ~ 1, thetat ~ 1)
  pr <- default_prior(formula, data, model)
  expect_equal(prior_for(pr, "kappa", "Intercept"), "normal(2, 1)")
  expect_equal(prior_for(pr, "thetat", "Intercept"), "logistic(0, 1)")
  # an intercept-only formula has no population-level coefficients to put an
  # effects prior on, so there is no b row at all
  expect_equal(prior_for(pr, "kappa", "b"), character(0))
  expect_equal(prior_for(pr, "thetat", "b"), character(0))
  # the location parameter is fixed, and brms writes the response parameter's
  # own priors with an empty dpar
  expect_equal(prior_for(pr, "", "Intercept"), "constant(0)")

  # 1 fixed effect + intercept
  formula <- bmf(kappa ~ set_size, thetat ~ set_size)
  pr <- default_prior(formula, data, model)
  expect_equal(prior_for(pr, "kappa", "b"), "normal(0, 1)")
  expect_equal(prior_for(pr, "thetat", "b"), "")
  expect_equal(prior_for(pr, "kappa", "Intercept"), "normal(2, 1)")
  expect_equal(prior_for(pr, "thetat", "Intercept"), "logistic(0, 1)")
  expect_equal(prior_for(pr, "", "Intercept"), "constant(0)")

  # 1 fixed effect intercept suppressed
  formula <- bmf(kappa ~ 0 + set_size, thetat ~ 0 + set_size)
  pr <- default_prior(formula, data, model)
  expect_equal(prior_for(pr, "kappa", "b"), "normal(2, 1)")
  expect_equal(prior_for(pr, "thetat", "b"), "logistic(0, 1)")
  expect_equal(prior_for(pr, "kappa", "Intercept"), character(0))
  expect_equal(prior_for(pr, "", "Intercept"), "constant(0)")

  # 2 fixed effects + intercept
  formula <- bmf(kappa ~ set_size + session, thetat ~ set_size + session)
  pr <- default_prior(formula, data, model)
  expect_equal(prior_for(pr, "kappa", "Intercept"), "normal(2, 1)")
  expect_equal(prior_for(pr, "thetat", "Intercept"), "logistic(0, 1)")
  expect_equal(prior_for(pr, "kappa", "b"), "normal(0, 1)")
  expect_equal(prior_for(pr, "thetat", "b"), "")

  # 2 fixed effects + intercept suppressed
  formula <- bmf(kappa ~ 0 + set_size + session, thetat ~ 0 + set_size + session)
  pr <- default_prior(formula, data, model)
  expect_equal(prior_for(pr, "kappa", "Intercept"), character(0))
  expect_equal(prior_for(pr, "kappa", "b"), "normal(0, 1)")
  expect_equal(prior_for(pr, "thetat", "b"), "")
  expect_equal(prior_for(pr, "kappa", "b", "set_size1"), "normal(2, 1)")
  expect_equal(prior_for(pr, "thetat", "b", "set_size1"), "logistic(0, 1)")
  expect_equal(prior_for(pr, "kappa", "b", "session2"), "")

  # 2 fixed effects + interaction + intercept
  formula <- bmf(kappa ~ set_size * session, thetat ~ set_size * session)
  pr <- default_prior(formula, data, model)
  expect_equal(prior_for(pr, "kappa", "Intercept"), "normal(2, 1)")
  expect_equal(prior_for(pr, "thetat", "Intercept"), "logistic(0, 1)")
  expect_equal(prior_for(pr, "kappa", "b"), "normal(0, 1)")

  # 2 fixed effects + interaction + intercept suppressed
  formula <- bmf(kappa ~ 0 + set_size * session, thetat ~ 0 + set_size * session)
  pr <- default_prior(formula, data, model)
  expect_equal(prior_for(pr, "kappa", "Intercept"), character(0))
  expect_equal(prior_for(pr, "kappa", "b", "set_size1"), "normal(2, 1)")
  expect_equal(prior_for(pr, "thetat", "b", "set_size1"), "logistic(0, 1)")
  expect_equal(prior_for(pr, "kappa", "b", "session2"), "")
  expect_equal(prior_for(pr, "kappa", "b", "set_size2:session2"), "")

  # interaction only between 2 fixed effects
  formula <- bmf(kappa ~ 0 + set_size:session, thetat ~ 0 + set_size:session)
  pr <- default_prior(formula, data, model)
  expect_equal(prior_for(pr, "kappa", "Intercept"), character(0))
  expect_equal(prior_for(pr, "kappa", "b"), "normal(2, 1)")
  expect_equal(prior_for(pr, "thetat", "b"), "logistic(0, 1)")
})


test_that("default priors are set correctly with random effects", {
  data <- oberauer_lin_2017
  model <- mixture2p("dev_rad")

  # Intercept only
  formula <- bmf(kappa ~ 1 + (1 | ID), thetat ~ 1 + (1 | ID))
  pr <- default_prior(formula, data, model)
  expect_equal(prior_for(pr, "kappa", "Intercept"), "normal(2, 1)")
  expect_equal(prior_for(pr, "thetat", "Intercept"), "logistic(0, 1)")
  expect_equal(prior_for(pr, "kappa", "b"), character(0))
  expect_equal(unique(pr[pr$class == "sd", ]$prior), c("student_t(3, 0, 2.5)", ""))

  # 1 fixed effect + intercept
  formula <- bmf(kappa ~ set_size + (1 | ID), thetat ~ set_size + (1 | ID))
  pr <- default_prior(formula, data, model)
  expect_equal(prior_for(pr, "kappa", "b"), "normal(0, 1)")
  expect_equal(prior_for(pr, "thetat", "b"), "")
  expect_equal(prior_for(pr, "kappa", "Intercept"), "normal(2, 1)")
  expect_equal(prior_for(pr, "thetat", "Intercept"), "logistic(0, 1)")
  expect_equal(prior_for(pr, "", "Intercept"), "constant(0)")
  expect_equal(unique(pr[pr$class == "sd", ]$prior), c("student_t(3, 0, 2.5)", ""))

  # 1 fixed effect intercept suppressed
  formula <- bmf(kappa ~ 0 + set_size + (1 | ID), thetat ~ 0 + set_size + (1 | ID))
  pr <- default_prior(formula, data, model)
  expect_equal(prior_for(pr, "kappa", "b"), "normal(2, 1)")
  expect_equal(prior_for(pr, "thetat", "b"), "logistic(0, 1)")
  expect_equal(prior_for(pr, "kappa", "Intercept"), character(0))
  expect_equal(prior_for(pr, "", "Intercept"), "constant(0)")
  expect_equal(unique(pr[pr$class == "sd", ]$prior), c("student_t(3, 0, 2.5)", ""))

  # 2 fixed effects + intercept
  formula <- bmf(kappa ~ set_size + session + (1 | ID), thetat ~ set_size + session + (1 | ID))
  pr <- default_prior(formula, data, model)
  expect_equal(prior_for(pr, "kappa", "Intercept"), "normal(2, 1)")
  expect_equal(prior_for(pr, "thetat", "Intercept"), "logistic(0, 1)")
  expect_equal(prior_for(pr, "kappa", "b"), "normal(0, 1)")
  expect_equal(unique(pr[pr$class == "sd", ]$prior), c("student_t(3, 0, 2.5)", ""))

  # 2 fixed effects + intercept suppressed
  formula <- bmf(kappa ~ 0 + set_size + session + (1 | ID), thetat ~ 0 + set_size + session + (1 | ID))
  pr <- default_prior(formula, data, model)
  expect_equal(prior_for(pr, "kappa", "Intercept"), character(0))
  expect_equal(prior_for(pr, "kappa", "b"), "normal(0, 1)")
  expect_equal(prior_for(pr, "kappa", "b", "set_size1"), "normal(2, 1)")
  expect_equal(prior_for(pr, "thetat", "b", "set_size1"), "logistic(0, 1)")
  expect_equal(prior_for(pr, "kappa", "b", "session2"), "")
  expect_equal(unique(pr[pr$class == "sd", ]$prior), c("student_t(3, 0, 2.5)", ""))

  # 2 fixed effects + interaction + intercept
  formula <- bmf(kappa ~ set_size * session + (1 | ID), thetat ~ set_size * session + (1 | ID))
  pr <- default_prior(formula, data, model)
  expect_equal(prior_for(pr, "kappa", "Intercept"), "normal(2, 1)")
  expect_equal(prior_for(pr, "thetat", "Intercept"), "logistic(0, 1)")
  expect_equal(prior_for(pr, "kappa", "b"), "normal(0, 1)")
  expect_equal(unique(pr[pr$class == "sd", ]$prior), c("student_t(3, 0, 2.5)", ""))

  # 2 fixed effects + interaction + intercept suppressed
  formula <- bmf(kappa ~ 0 + set_size * session + (1 | ID), thetat ~ 0 + set_size * session + (1 | ID))
  pr <- default_prior(formula, data, model)
  expect_equal(prior_for(pr, "kappa", "Intercept"), character(0))
  expect_equal(prior_for(pr, "kappa", "b", "set_size1"), "normal(2, 1)")
  expect_equal(prior_for(pr, "thetat", "b", "set_size1"), "logistic(0, 1)")
  expect_equal(prior_for(pr, "kappa", "b", "session2"), "")
  expect_equal(prior_for(pr, "kappa", "b", "set_size2:session2"), "")
  expect_equal(unique(pr[pr$class == "sd", ]$prior), c("student_t(3, 0, 2.5)", ""))

  # interaction only between 2 fixed effects
  formula <- bmf(kappa ~ 0 + set_size:session + (1 | ID), thetat ~ 0 + set_size:session + (1 | ID))
  pr <- default_prior(formula, data, model)
  expect_equal(prior_for(pr, "kappa", "Intercept"), character(0))
  expect_equal(prior_for(pr, "kappa", "b"), "normal(2, 1)")
  expect_equal(prior_for(pr, "thetat", "b"), "logistic(0, 1)")
  expect_equal(unique(pr[pr$class == "sd", ]$prior), c("student_t(3, 0, 2.5)", ""))
})


test_that("default priors are set correctly with fixed effects only and sdm model", {
  data <- oberauer_lin_2017
  model <- sdm("dev_rad")

  # Intercept only
  formula <- bmf(kappa ~ 1, c ~ 1)
  pr <- default_prior(formula, data, model)
  expect_equal(pr[pr$class == "Intercept", ]$prior, c("student_t(5, 2, 0.75)", "student_t(5, 1.75, 0.75)", "constant(0)"))
  expect_equal(pr[pr$coef == "" & pr$class == "b", ]$prior, character(0))

  # 1 fixed effect + intercept
  formula <- bmf(kappa ~ set_size, c ~ set_size)
  pr <- default_prior(formula, data, model)
  expect_equal(pr[pr$coef == "" & pr$class == "b", ]$prior, c("normal(0, 1)", "normal(0, 1)"))
  expect_equal(pr[pr$class == "Intercept", ]$prior, c("student_t(5, 2, 0.75)", "student_t(5, 1.75, 0.75)", "constant(0)"))

  # 1 fixed effect intercept suppressed
  formula <- bmf(kappa ~ 0 + set_size, c ~ 0 + set_size)
  pr <- default_prior(formula, data, model)
  expect_equal(pr[pr$coef == "" & pr$class == "b", ]$prior, c("student_t(5, 2, 0.75)", "student_t(5, 1.75, 0.75)"))
  expect_equal(pr[pr$class == "Intercept", ]$prior, c("constant(0)"))

  # 2 fixed effects + intercept
  formula <- bmf(kappa ~ set_size + session, c ~ set_size + session)
  pr <- default_prior(formula, data, model)
  expect_equal(pr[pr$class == "Intercept", ]$prior, c("student_t(5, 2, 0.75)", "student_t(5, 1.75, 0.75)", "constant(0)"))
  expect_equal(pr[pr$coef == "" & pr$class == "b", ]$prior, c("normal(0, 1)", "normal(0, 1)"))

  # 2 fixed effects + intercept suppressed
  formula <- bmf(kappa ~ 0 + set_size + session, c ~ 0 + set_size + session)
  pr <- default_prior(formula, data, model)
  expect_equal(pr[pr$class == "Intercept", ]$prior, "constant(0)")
  expect_equal(pr[pr$coef == "" & pr$class == "b", ]$prior, c("normal(0, 1)", "normal(0, 1)"))
  expect_equal(pr[pr$coef == "set_size1" & pr$class == "b", ]$prior, c("student_t(5, 2, 0.75)", "student_t(5, 1.75, 0.75)"))
  expect_equal(pr[pr$coef == "session2" & pr$class == "b", ]$prior, c("", ""))

  # 2 fixed effects + interaction + intercept
  formula <- bmf(kappa ~ set_size * session, c ~ set_size * session)
  pr <- default_prior(formula, data, model)
  expect_equal(pr[pr$class == "Intercept", ]$prior, c("student_t(5, 2, 0.75)", "student_t(5, 1.75, 0.75)", "constant(0)"))
  expect_equal(pr[pr$coef == "" & pr$class == "b", ]$prior, c("normal(0, 1)", "normal(0, 1)"))

  # 2 fixed effects + interaction + intercept suppressed
  formula <- bmf(kappa ~ 0 + set_size * session, c ~ 0 + set_size * session)
  pr <- default_prior(formula, data, model)
  expect_equal(pr[pr$class == "Intercept", ]$prior, c("constant(0)"))
  expect_equal(pr[pr$coef == "" & pr$class == "b", ]$prior, c("normal(0, 1)", "normal(0, 1)"))
  expect_equal(pr[pr$coef == "set_size1" & pr$class == "b", ]$prior, c("student_t(5, 2, 0.75)", "student_t(5, 1.75, 0.75)"))
  expect_equal(pr[pr$coef == "session2" & pr$class == "b", ]$prior, c("", ""))
  expect_equal(pr[pr$coef == "set_size2:session2" & pr$class == "b", ]$prior, c("", ""))

  # interaction only between 2 fixed effects
  formula <- bmf(kappa ~ 0 + set_size:session, c ~ 0 + set_size:session)
  pr <- default_prior(formula, data, model)
  expect_equal(pr[pr$class == "Intercept", ]$prior, "constant(0)")
  expect_equal(pr[pr$coef == "" & pr$class == "b", ]$prior, c("student_t(5, 2, 0.75)", "student_t(5, 1.75, 0.75)"))
})


test_that("default priors work when there are no fixed parameters", {
  formula <- bmf(
    mu ~ 1,
    c ~ 1,
    kappa ~ 1
  )

  pr <- default_prior(formula, oberauer_lin_2017, sdm("dev_rad"))
  expect_s3_class(pr, "brmsprior")
})

test_that("default priors work when there are non-linear transformations of default parameters", {
  withr::local_options(bmm.silent = 2)
  expect_warning(
    dp <- default_prior(
      object = bmmformula(c ~ exp(nlc), nlc ~ 1),
      data = oberauer_lin_2017,
      model = sdm(resp_error = "dev_rad")
    ),
    "Non-linear transformations"
  )
  expect_true(!("c" %in% dp$dpar))
  expect_true("nlc" %in% dp$nlpar)
})
