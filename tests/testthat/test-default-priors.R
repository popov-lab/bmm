test_that("default priors are set correctly with fixed effects only", {
  data <- oberauer_lin_2017
  model <- mixture2p("dev_rad")

  # Intercept only
  formula <- bmf(kappa ~ 1, thetat ~ 1)
  pr <- default_prior(formula, data, model)
  expect_equal(pr[pr$coef == "Intercept", ]$prior, c("normal(2, 1)", "logistic(0, 1)"))
  expect_equal(pr[pr$coef == "" & pr$class == "b", ]$prior, c("", ""))

  # 1 fixed effect + intercept
  formula <- bmf(kappa ~ set_size, thetat ~ set_size)
  pr <- default_prior(formula, data, model)
  expect_equal(pr[pr$coef == "" & pr$class == "b", ]$prior, c("", "normal(0, 1)"))
  expect_equal(pr[pr$coef == "Intercept", ]$prior, c("normal(2, 1)", "logistic(0, 1)"))
  expect_true(all(grepl("constant", pr[pr$dpar %in% c("mu1", "mu2", "kappa2"), ]$prior)))

  # 1 fixed effect intercept suppressed
  formula <- bmf(kappa ~ 0 + set_size, thetat ~ 0 + set_size)
  pr <- default_prior(formula, data, model)
  expect_equal(pr[pr$coef == "" & pr$class == "b", ]$prior, c("normal(2, 1)", "logistic(0, 1)"))
  expect_equal(pr[pr$coef == "Intercept", ]$prior, character(0))
  expect_true(all(grepl("constant", pr[pr$dpar %in% c("mu1", "mu2", "kappa2"), ]$prior)))

  # 2 fixed effects + intercept
  formula <- bmf(kappa ~ set_size + session, thetat ~ set_size + session)
  pr <- default_prior(formula, data, model)
  expect_equal(pr[pr$coef == "Intercept", ]$prior, c("normal(2, 1)", "logistic(0, 1)"))
  expect_equal(pr[pr$coef == "" & pr$class == "b", ]$prior, c("", "normal(0, 1)"))

  # 2 fixed effects + intercept suppressed
  formula <- bmf(kappa ~ 0 + set_size + session, thetat ~ 0 + set_size + session)
  pr <- default_prior(formula, data, model)
  expect_equal(pr[pr$coef == "Intercept", ]$prior, character(0))
  expect_equal(pr[pr$coef == "" & pr$class == "b", ]$prior, c("", "normal(0, 1)"))
  expect_equal(pr[pr$coef == "set_size1" & pr$class == "b", ]$prior, c("normal(2, 1)", "logistic(0, 1)"))
  expect_equal(pr[pr$coef == "session2" & pr$class == "b", ]$prior, c("", ""))

  # 2 fixed effects + interaction + intercept
  formula <- bmf(kappa ~ set_size * session, thetat ~ set_size * session)
  pr <- default_prior(formula, data, model)
  expect_equal(pr[pr$coef == "Intercept", ]$prior, c("normal(2, 1)", "logistic(0, 1)"))
  expect_equal(pr[pr$coef == "" & pr$class == "b", ]$prior, c("", "normal(0, 1)"))

  # 2 fixed effects + interaction + intercept suppressed
  formula <- bmf(kappa ~ 0 + set_size * session, thetat ~ 0 + set_size * session)
  pr <- default_prior(formula, data, model)
  expect_equal(pr[pr$coef == "Intercept", ]$prior, character(0))
  expect_equal(pr[pr$coef == "" & pr$class == "b", ]$prior, c("", "normal(0, 1)"))
  expect_equal(pr[pr$coef == "set_size1" & pr$class == "b", ]$prior, c("normal(2, 1)", "logistic(0, 1)"))
  expect_equal(pr[pr$coef == "session2" & pr$class == "b", ]$prior, c("", ""))
  expect_equal(pr[pr$coef == "set_size2:session2" & pr$class == "b", ]$prior, c("", ""))

  # interaction only between 2 fixed effects
  formula <- bmf(kappa ~ 0 + set_size:session, thetat ~ 0 + set_size:session)
  pr <- default_prior(formula, data, model)
  expect_equal(pr[pr$coef == "Intercept", ]$prior, character(0))
  expect_equal(pr[pr$coef == "" & pr$class == "b", ]$prior, c("normal(2, 1)", "logistic(0, 1)"))
})


test_that("default priors are set correctly with random effects", {
  data <- oberauer_lin_2017
  model <- mixture2p("dev_rad")

  # Intercept only
  formula <- bmf(kappa ~ 1 + (1 | ID), thetat ~ 1 + (1 | ID))
  pr <- default_prior(formula, data, model)
  expect_equal(pr[pr$coef == "Intercept" & pr$class == "b", ]$prior, c("normal(2, 1)", "logistic(0, 1)"))
  expect_equal(pr[pr$coef == "" & pr$class == "b", ]$prior, c("", ""))
  expect_equal(unique(pr[pr$class == "sd", ]$prior), c("", "exponential(1)"))

  # 1 fixed effect + intercept
  formula <- bmf(kappa ~ set_size + (1 | ID), thetat ~ set_size + (1 | ID))
  pr <- default_prior(formula, data, model)
  expect_equal(pr[pr$coef == "" & pr$class == "b", ]$prior, c("", "normal(0, 1)"))
  expect_equal(pr[pr$coef == "Intercept" & pr$class == "b", ]$prior, c("normal(2, 1)", "logistic(0, 1)"))
  expect_true(all(grepl("constant", pr[pr$dpar %in% c("mu1", "mu2", "kappa2"), ]$prior)))
  expect_equal(unique(pr[pr$class == "sd", ]$prior), c("", "exponential(1)"))

  # 1 fixed effect intercept suppressed
  formula <- bmf(kappa ~ 0 + set_size + (1 | ID), thetat ~ 0 + set_size + (1 | ID))
  pr <- default_prior(formula, data, model)
  expect_equal(pr[pr$coef == "" & pr$class == "b", ]$prior, c("normal(2, 1)", "logistic(0, 1)"))
  expect_equal(pr[pr$coef == "Intercept" & pr$class == "b", ]$prior, character(0))
  expect_true(all(grepl("constant", pr[pr$dpar %in% c("mu1", "mu2", "kappa2"), ]$prior)))
  expect_equal(unique(pr[pr$class == "sd", ]$prior), c("", "exponential(1)"))

  # 2 fixed effects + intercept
  formula <- bmf(kappa ~ set_size + session + (1 | ID), thetat ~ set_size + session + (1 | ID))
  pr <- default_prior(formula, data, model)
  expect_equal(pr[pr$coef == "Intercept" & pr$class == "b", ]$prior, c("normal(2, 1)", "logistic(0, 1)"))
  expect_equal(pr[pr$coef == "" & pr$class == "b", ]$prior, c("", "normal(0, 1)"))
  expect_equal(unique(pr[pr$class == "sd", ]$prior), c("", "exponential(1)"))

  # 2 fixed effects + intercept suppressed
  formula <- bmf(kappa ~ 0 + set_size + session + (1 | ID), thetat ~ 0 + set_size + session + (1 | ID))
  pr <- default_prior(formula, data, model)
  expect_equal(pr[pr$coef == "Intercept" & pr$class == "b", ]$prior, character(0))
  expect_equal(pr[pr$coef == "" & pr$class == "b", ]$prior, c("", "normal(0, 1)"))
  expect_equal(pr[pr$coef == "set_size1" & pr$class == "b", ]$prior, c("normal(2, 1)", "logistic(0, 1)"))
  expect_equal(pr[pr$coef == "session2" & pr$class == "b", ]$prior, c("", ""))
  expect_equal(unique(pr[pr$class == "sd", ]$prior), c("", "exponential(1)"))

  # 2 fixed effects + interaction + intercept
  formula <- bmf(kappa ~ set_size * session + (1 | ID), thetat ~ set_size * session + (1 | ID))
  pr <- default_prior(formula, data, model)
  expect_equal(pr[pr$coef == "Intercept" & pr$class == "b", ]$prior, c("normal(2, 1)", "logistic(0, 1)"))
  expect_equal(pr[pr$coef == "" & pr$class == "b", ]$prior, c("", "normal(0, 1)"))
  expect_equal(unique(pr[pr$class == "sd", ]$prior), c("", "exponential(1)"))

  # 2 fixed effects + interaction + intercept suppressed
  formula <- bmf(kappa ~ 0 + set_size * session + (1 | ID), thetat ~ 0 + set_size * session + (1 | ID))
  pr <- default_prior(formula, data, model)
  expect_equal(pr[pr$coef == "Intercept" & pr$class == "b", ]$prior, character(0))
  expect_equal(pr[pr$coef == "" & pr$class == "b", ]$prior, c("", "normal(0, 1)"))
  expect_equal(pr[pr$coef == "set_size1" & pr$class == "b", ]$prior, c("normal(2, 1)", "logistic(0, 1)"))
  expect_equal(pr[pr$coef == "session2" & pr$class == "b", ]$prior, c("", ""))
  expect_equal(pr[pr$coef == "set_size2:session2" & pr$class == "b", ]$prior, c("", ""))
  expect_equal(unique(pr[pr$class == "sd", ]$prior), c("", "exponential(1)"))

  # interaction only between 2 fixed effects
  formula <- bmf(kappa ~ 0 + set_size:session + (1 | ID), thetat ~ 0 + set_size:session + (1 | ID))
  pr <- default_prior(formula, data, model)
  expect_equal(pr[pr$coef == "Intercept" & pr$class == "b", ]$prior, character(0))
  expect_equal(pr[pr$coef == "" & pr$class == "b", ]$prior, c("normal(2, 1)", "logistic(0, 1)"))
  expect_equal(unique(pr[pr$class == "sd", ]$prior), c("", "exponential(1)"))
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
