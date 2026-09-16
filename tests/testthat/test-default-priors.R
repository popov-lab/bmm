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
  expect_equal(unique(pr[pr$class == "sd", ]$prior), c("student_t(3, 0, 2.5)", ""))

  # 1 fixed effect + intercept
  formula <- bmf(kappa ~ set_size + (1 | ID), thetat ~ set_size + (1 | ID))
  pr <- default_prior(formula, data, model)
  expect_equal(pr[pr$coef == "" & pr$class == "b", ]$prior, c("", "normal(0, 1)"))
  expect_equal(pr[pr$coef == "Intercept" & pr$class == "b", ]$prior, c("normal(2, 1)", "logistic(0, 1)"))
  expect_true(all(grepl("constant", pr[pr$dpar %in% c("mu1", "mu2", "kappa2"), ]$prior)))
  expect_equal(unique(pr[pr$class == "sd", ]$prior), c("student_t(3, 0, 2.5)", ""))

  # 1 fixed effect intercept suppressed
  formula <- bmf(kappa ~ 0 + set_size + (1 | ID), thetat ~ 0 + set_size + (1 | ID))
  pr <- default_prior(formula, data, model)
  expect_equal(pr[pr$coef == "" & pr$class == "b", ]$prior, c("normal(2, 1)", "logistic(0, 1)"))
  expect_equal(pr[pr$coef == "Intercept" & pr$class == "b", ]$prior, character(0))
  expect_true(all(grepl("constant", pr[pr$dpar %in% c("mu1", "mu2", "kappa2"), ]$prior)))
  expect_equal(unique(pr[pr$class == "sd", ]$prior), c("student_t(3, 0, 2.5)", ""))

  # 2 fixed effects + intercept
  formula <- bmf(kappa ~ set_size + session + (1 | ID), thetat ~ set_size + session + (1 | ID))
  pr <- default_prior(formula, data, model)
  expect_equal(pr[pr$coef == "Intercept" & pr$class == "b", ]$prior, c("normal(2, 1)", "logistic(0, 1)"))
  expect_equal(pr[pr$coef == "" & pr$class == "b", ]$prior, c("", "normal(0, 1)"))
  expect_equal(unique(pr[pr$class == "sd", ]$prior), c("student_t(3, 0, 2.5)", ""))

  # 2 fixed effects + intercept suppressed
  formula <- bmf(kappa ~ 0 + set_size + session + (1 | ID), thetat ~ 0 + set_size + session + (1 | ID))
  pr <- default_prior(formula, data, model)
  expect_equal(pr[pr$coef == "Intercept" & pr$class == "b", ]$prior, character(0))
  expect_equal(pr[pr$coef == "" & pr$class == "b", ]$prior, c("", "normal(0, 1)"))
  expect_equal(pr[pr$coef == "set_size1" & pr$class == "b", ]$prior, c("normal(2, 1)", "logistic(0, 1)"))
  expect_equal(pr[pr$coef == "session2" & pr$class == "b", ]$prior, c("", ""))
  expect_equal(unique(pr[pr$class == "sd", ]$prior), c("student_t(3, 0, 2.5)", ""))

  # 2 fixed effects + interaction + intercept
  formula <- bmf(kappa ~ set_size * session + (1 | ID), thetat ~ set_size * session + (1 | ID))
  pr <- default_prior(formula, data, model)
  expect_equal(pr[pr$coef == "Intercept" & pr$class == "b", ]$prior, c("normal(2, 1)", "logistic(0, 1)"))
  expect_equal(pr[pr$coef == "" & pr$class == "b", ]$prior, c("", "normal(0, 1)"))
  expect_equal(unique(pr[pr$class == "sd", ]$prior), c("student_t(3, 0, 2.5)", ""))

  # 2 fixed effects + interaction + intercept suppressed
  formula <- bmf(kappa ~ 0 + set_size * session + (1 | ID), thetat ~ 0 + set_size * session + (1 | ID))
  pr <- default_prior(formula, data, model)
  expect_equal(pr[pr$coef == "Intercept" & pr$class == "b", ]$prior, character(0))
  expect_equal(pr[pr$coef == "" & pr$class == "b", ]$prior, c("", "normal(0, 1)"))
  expect_equal(pr[pr$coef == "set_size1" & pr$class == "b", ]$prior, c("normal(2, 1)", "logistic(0, 1)"))
  expect_equal(pr[pr$coef == "session2" & pr$class == "b", ]$prior, c("", ""))
  expect_equal(pr[pr$coef == "set_size2:session2" & pr$class == "b", ]$prior, c("", ""))
  expect_equal(unique(pr[pr$class == "sd", ]$prior), c("student_t(3, 0, 2.5)", ""))

  # interaction only between 2 fixed effects
  formula <- bmf(kappa ~ 0 + set_size:session + (1 | ID), thetat ~ 0 + set_size:session + (1 | ID))
  pr <- default_prior(formula, data, model)
  expect_equal(pr[pr$coef == "Intercept" & pr$class == "b", ]$prior, character(0))
  expect_equal(pr[pr$coef == "" & pr$class == "b", ]$prior, c("normal(2, 1)", "logistic(0, 1)"))
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

sd_rows <- function(pr) {
  pr[pr$class == "sd" & pr$coef == "" & pr$group == "" & pr$prior != "", ]
}

test_that("an sd default becomes one blanket prior per parameter with random effects", {
  data <- oberauer_lin_2017
  model <- mixture2p("dev_rad")
  model$default_priors$kappa$sd <- "exponential(1)"
  model$default_priors$thetat$sd <- "exponential(2)"

  # random intercepts on both parameters
  pr <- default_prior(bmf(kappa ~ 1 + (1 | ID), thetat ~ 1 + (1 | ID)), data, model)
  sd_pr <- sd_rows(pr)
  expect_equal(nrow(sd_pr), 2)
  expect_equal(sd_pr[sd_pr$nlpar == "kappa", ]$prior, "exponential(1)")
  expect_equal(sd_pr[sd_pr$nlpar == "thetat", ]$prior, "exponential(2)")

  # random intercept + slope: still a single blanket row
  pr <- default_prior(bmf(kappa ~ set_size + (set_size | ID), thetat ~ 1), data, model)
  sd_pr <- sd_rows(pr)
  expect_equal(nrow(sd_pr), 1)
  expect_equal(sd_pr$nlpar, "kappa")
  expect_false(any(pr$class == "sd" & pr$nlpar == "thetat"))

  # random slope without random intercept
  pr <- default_prior(bmf(kappa ~ set_size + (0 + set_size | ID), thetat ~ 1), data, model)
  expect_equal(sd_rows(pr)$prior, "exponential(1)")

  # crossed grouping factors share the one blanket row
  pr <- default_prior(bmf(kappa ~ 1 + (1 | ID) + (1 | session), thetat ~ 1), data, model)
  expect_equal(nrow(sd_rows(pr)), 1)
  expect_setequal(unique(pr[pr$class == "sd" & pr$nlpar == "kappa", ]$group), c("", "ID", "session"))
})

test_that("no sd prior is emitted for a parameter without random effects", {
  data <- oberauer_lin_2017
  model <- mixture2p("dev_rad")
  model$default_priors$kappa$sd <- "exponential(1)"

  pr <- default_prior(bmf(kappa ~ 1, thetat ~ 1), data, model)
  expect_false(any(pr$class == "sd"))

  # brms rejects priors on non-existent parameters only when the model is built
  fit <- bmm(bmf(kappa ~ 1, thetat ~ 1), data, model,
    backend = "mock", mock_fit = 1, rename = FALSE
  )
  expect_false(any(fit$prior$class == "sd"))
})

test_that("sd defaults are set for distributional parameters", {
  data <- oberauer_lin_2017
  model <- sdm("dev_rad")
  model$default_priors$c$sd <- "exponential(1)"

  pr <- default_prior(bmf(c ~ 1 + (1 | ID), kappa ~ 1 + (1 | ID)), data, model)
  sd_pr <- sd_rows(pr)
  expect_equal(sd_pr[sd_pr$dpar == "c", ]$prior, "exponential(1)")
  # brms keeps its own dpar-level default where bmm sets none
  expect_equal(sd_pr[sd_pr$dpar == "kappa", ]$prior, "student_t(3, 0, 2.5)")
})

test_that("user sd priors take precedence over the sd default", {
  data <- oberauer_lin_2017
  model <- mixture2p("dev_rad")
  model$default_priors$kappa$sd <- "exponential(1)"
  formula <- bmf(kappa ~ 1 + (1 | ID), thetat ~ 1)

  # same key as the default: the default row is replaced
  fit <- bmm(formula, data, model,
    prior = brms::prior_("normal(0, 0.5)", class = "sd", nlpar = "kappa"),
    backend = "mock", mock_fit = 1, rename = FALSE
  )
  expect_equal(sd_rows(fit$prior)$prior, "normal(0, 0.5)")
  expect_false("exponential(1)" %in% fit$prior$prior)

  # more specific key: both rows survive and brms applies the user's
  fit <- bmm(formula, data, model,
    prior = brms::prior_("normal(0, 0.5)", class = "sd", group = "ID", nlpar = "kappa"),
    backend = "mock", mock_fit = 1, rename = FALSE
  )
  expect_true("exponential(1)" %in% fit$prior$prior)
  expect_match(brms::stancode(fit), "normal_lpdf(sd_1 | 0, 0.5)", fixed = TRUE)
  expect_no_match(brms::stancode(fit), "exponential_lpdf(sd_1", fixed = TRUE)
})

test_that("bmm.default_priors = FALSE also disables the sd defaults", {
  withr::local_options(bmm.default_priors = FALSE)
  model <- mixture2p("dev_rad")
  model$default_priors$kappa$sd <- "exponential(1)"

  pr <- default_prior(bmf(kappa ~ 1 + (1 | ID), thetat ~ 1), oberauer_lin_2017, model)
  expect_false("exponential(1)" %in% pr$prior)
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
