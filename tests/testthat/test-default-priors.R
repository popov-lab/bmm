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
  expect_equal(pr[pr$coef == "" & pr$class == "b", ]$prior, c("normal(0, 1)", "normal(0, 0.5)"))
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
  expect_equal(pr[pr$coef == "" & pr$class == "b", ]$prior, c("normal(0, 1)", "normal(0, 0.5)"))

  # 2 fixed effects + intercept suppressed
  formula <- bmf(kappa ~ 0 + set_size + session, thetat ~ 0 + set_size + session)
  pr <- default_prior(formula, data, model)
  expect_equal(pr[pr$coef == "Intercept", ]$prior, character(0))
  expect_equal(pr[pr$coef == "" & pr$class == "b", ]$prior, c("normal(0, 1)", "normal(0, 0.5)"))
  expect_equal(pr[pr$coef == "set_size1" & pr$class == "b", ]$prior, c("normal(2, 1)", "logistic(0, 1)"))
  expect_equal(pr[pr$coef == "session2" & pr$class == "b", ]$prior, c("normal(0, 1)", "normal(0, 0.5)"))

  # 2 fixed effects + interaction + intercept
  formula <- bmf(kappa ~ set_size * session, thetat ~ set_size * session)
  pr <- default_prior(formula, data, model)
  expect_equal(pr[pr$coef == "Intercept", ]$prior, c("normal(2, 1)", "logistic(0, 1)"))
  expect_equal(pr[pr$coef == "" & pr$class == "b", ]$prior, c("normal(0, 1)", "normal(0, 0.5)"))

  # 2 fixed effects + interaction + intercept suppressed
  formula <- bmf(kappa ~ 0 + set_size * session, thetat ~ 0 + set_size * session)
  pr <- default_prior(formula, data, model)
  expect_equal(pr[pr$coef == "Intercept", ]$prior, character(0))
  expect_equal(pr[pr$coef == "" & pr$class == "b", ]$prior, c("normal(0, 1)", "normal(0, 0.5)"))
  expect_equal(pr[pr$coef == "set_size1" & pr$class == "b", ]$prior, c("normal(2, 1)", "logistic(0, 1)"))
  expect_equal(pr[pr$coef == "session2" & pr$class == "b", ]$prior, c("normal(0, 1)", "normal(0, 0.5)"))
  expect_equal(pr[pr$coef == "set_size2:session2" & pr$class == "b", ]$prior, c("", ""))

  # interaction only between 2 fixed effects
  formula <- bmf(kappa ~ 0 + set_size:session, thetat ~ 0 + set_size:session)
  pr <- default_prior(formula, data, model)
  expect_equal(pr[pr$coef == "Intercept", ]$prior, character(0))
  expect_equal(pr[pr$coef == "" & pr$class == "b", ]$prior, c("normal(2, 1)", "logistic(0, 1)"))

  # Intercept + interaction only (no main effects)
  formula <- bmf(kappa ~ 1 + set_size:session, thetat ~ 1 + set_size:session)
  pr <- default_prior(formula, data, model)
  expect_equal(pr[pr$coef == "Intercept", ]$prior, c("normal(2, 1)", "logistic(0, 1)"))
  # Class-level effects prior should be set for interactions
  expect_equal(pr[pr$coef == "" & pr$class == "b", ]$prior, c("normal(0, 1)", "normal(0, 0.5)"))
})


test_that("default priors handle all predictor levels when intercept suppressed", {
  model <- mixture2p("dev_rad")
  data <- oberauer_lin_2017
  
  # With intercept suppressed and single predictor, class-level main prior is set
  formula <- bmf(kappa ~ 0 + set_size, thetat ~ 0 + set_size)
  pr <- default_prior(formula, data, model)
  
  # Class-level b prior should use main prior (not effects)
  expect_equal(pr[pr$coef == "" & pr$class == "b", ]$prior, c("normal(2, 1)", "logistic(0, 1)"))
  
  # With two predictors, first predictor level gets explicit main prior
  formula <- bmf(kappa ~ 0 + set_size + session, thetat ~ 0 + set_size + session)
  pr <- default_prior(formula, data, model)
  
  # First predictor: first level gets explicit main prior
  expect_equal(pr[pr$coef == "set_size1" & pr$class == "b", ]$prior, c("normal(2, 1)", "logistic(0, 1)"))
  
  # Second predictor: second level gets explicit effects prior  
  expect_equal(pr[pr$coef == "session2" & pr$class == "b", ]$prior, c("normal(0, 1)", "normal(0, 0.5)"))
  
  # Class-level effects prior is set for non-first predictor coefficients
  expect_equal(pr[pr$coef == "" & pr$class == "b", ]$prior, c("normal(0, 1)", "normal(0, 0.5)"))
})


test_that("default priors distinguish main vs effects priors correctly", {
  model <- mixture2p("dev_rad")
  data <- oberauer_lin_2017
  
  # Pattern 1: Intercept suppressed, single predictor → class-level MAIN prior
  formula <- bmf(kappa ~ 0 + session, thetat ~ 0 + session)
  pr <- default_prior(formula, data, model)
  
  # Class-level b prior should use main prior (not effects)
  expect_equal(pr[pr$coef == "" & pr$class == "b", ]$prior, c("normal(2, 1)", "logistic(0, 1)"))
  # Individual coefficients inherit from class-level (show as empty)
  expect_equal(pr[pr$coef == "session1" & pr$class == "b", ]$prior, c("", ""))
  expect_equal(pr[pr$coef == "session2" & pr$class == "b", ]$prior, c("", ""))
  
  # Pattern 2: Intercept included, single predictor → class-level EFFECTS prior set
  formula <- bmf(kappa ~ 1 + session, thetat ~ 1 + session)
  pr <- default_prior(formula, data, model)
  
  expect_equal(pr[pr$coef == "Intercept", ]$prior, c("normal(2, 1)", "logistic(0, 1)"))
  expect_equal(pr[pr$coef == "" & pr$class == "b", ]$prior, c("normal(0, 1)", "normal(0, 0.5)"))
  # session2 inherits from class-level (shows as empty)
  expect_equal(pr[pr$coef == "session2" & pr$class == "b", ]$prior, c("", ""))
  
  # Pattern 3: Two predictors with intercept → class-level EFFECTS prior set
  formula <- bmf(kappa ~ 1 + set_size + session, thetat ~ 1 + set_size + session)
  pr <- default_prior(formula, data, model)
  
  expect_equal(pr[pr$coef == "Intercept", ]$prior, c("normal(2, 1)", "logistic(0, 1)"))
  expect_equal(pr[pr$coef == "" & pr$class == "b", ]$prior, c("normal(0, 1)", "normal(0, 0.5)"))
  # Both non-reference levels inherit from class-level
  expect_equal(pr[pr$coef == "set_size2" & pr$class == "b", ]$prior, c("", ""))
  expect_equal(pr[pr$coef == "session2" & pr$class == "b", ]$prior, c("", ""))
})


test_that("default priors handle interaction-only formulas correctly", {
  model <- mixture2p("dev_rad")
  data <- oberauer_lin_2017
  
  # Intercept + interaction only: class-level effects prior should be set
  formula <- bmf(kappa ~ 1 + set_size:session, thetat ~ 1 + set_size:session)
  pr <- default_prior(formula, data, model)
  
  expect_equal(pr[pr$coef == "Intercept", ]$prior, c("normal(2, 1)", "logistic(0, 1)"))
  # Critical: class-level b prior should be effects prior (not flat)
  expect_equal(pr[pr$coef == "" & pr$class == "b", ]$prior, c("normal(0, 1)", "normal(0, 0.5)"))
  
  # Verify individual interaction terms inherit (show as empty but class-level is set)
  interaction_priors <- pr[grepl(":", pr$coef) & pr$class == "b", ]
  expect_true(all(interaction_priors$prior == ""))
  
  # No intercept + interaction only: class-level should be MAIN prior
  formula <- bmf(kappa ~ 0 + set_size:session, thetat ~ 0 + set_size:session)
  pr <- default_prior(formula, data, model)
  
  expect_equal(pr[pr$coef == "" & pr$class == "b", ]$prior, c("normal(2, 1)", "logistic(0, 1)"))
})


test_that("default priors are consistent across different contrast codings", {
  model <- mixture2p("dev_rad")
  data_orig <- oberauer_lin_2017
  
  # Test with treatment coding (default)
  data_treatment <- data_orig
  contrasts(data_treatment$set_size) <- contr.treatment(nlevels(data_treatment$set_size))
  contrasts(data_treatment$session) <- contr.treatment(nlevels(data_treatment$session))
  
  # Test with sum coding
  data_sum <- data_orig
  contrasts(data_sum$set_size) <- contr.sum(nlevels(data_sum$set_size))
  contrasts(data_sum$session) <- contr.sum(nlevels(data_sum$session))
  
  # Test with helmert coding
  data_helmert <- data_orig
  contrasts(data_helmert$set_size) <- contr.helmert(nlevels(data_helmert$set_size))
  contrasts(data_helmert$session) <- contr.helmert(nlevels(data_helmert$session))
  
  # Optionally test with bayestestR::contr.equalprior if available
  data_equalprior <- NULL
  if (requireNamespace("bayestestR", quietly = TRUE)) {
    data_equalprior <- data_orig
    contr_equalprior <- getFromNamespace("contr.equalprior", "bayestestR")
    contrasts(data_equalprior$set_size) <- contr_equalprior(nlevels(data_equalprior$set_size))
    contrasts(data_equalprior$session) <- contr_equalprior(nlevels(data_equalprior$session))
  }
  
  # Test Pattern 1: Single predictor with intercept
  formula <- bmf(kappa ~ 1 + set_size, thetat ~ 1 + set_size)
  
  pr_treatment <- default_prior(formula, data_treatment, model)
  pr_sum <- default_prior(formula, data_sum, model)
  pr_helmert <- default_prior(formula, data_helmert, model)
  
  # All should have same class-level priors
  expect_equal(
    pr_treatment[pr_treatment$coef == "" & pr_treatment$class == "b", ]$prior,
    pr_sum[pr_sum$coef == "" & pr_sum$class == "b", ]$prior
  )
  expect_equal(
    pr_treatment[pr_treatment$coef == "" & pr_treatment$class == "b", ]$prior,
    pr_helmert[pr_helmert$coef == "" & pr_helmert$class == "b", ]$prior
  )
  expect_equal(
    pr_treatment[pr_treatment$coef == "Intercept", ]$prior,
    pr_sum[pr_sum$coef == "Intercept", ]$prior
  )
  
  if (!is.null(data_equalprior)) {
    pr_equalprior <- default_prior(formula, data_equalprior, model)
    expect_equal(
      pr_treatment[pr_treatment$coef == "" & pr_treatment$class == "b", ]$prior,
      pr_equalprior[pr_equalprior$coef == "" & pr_equalprior$class == "b", ]$prior
    )
  }
  
  # Test Pattern 2: Multiple predictors with intercept
  formula <- bmf(kappa ~ 1 + set_size + session, thetat ~ 1 + set_size + session)
  
  pr_treatment <- default_prior(formula, data_treatment, model)
  pr_sum <- default_prior(formula, data_sum, model)
  pr_helmert <- default_prior(formula, data_helmert, model)
  
  expect_equal(
    pr_treatment[pr_treatment$coef == "" & pr_treatment$class == "b", ]$prior,
    pr_sum[pr_sum$coef == "" & pr_sum$class == "b", ]$prior
  )
  expect_equal(
    pr_treatment[pr_treatment$coef == "" & pr_treatment$class == "b", ]$prior,
    pr_helmert[pr_helmert$coef == "" & pr_helmert$class == "b", ]$prior
  )
  expect_equal(
    pr_treatment[pr_treatment$coef == "Intercept", ]$prior,
    pr_sum[pr_sum$coef == "Intercept", ]$prior
  )
  
  # Test Pattern 3: Single predictor, no intercept
  formula <- bmf(kappa ~ 0 + set_size, thetat ~ 0 + set_size)
  
  pr_treatment <- default_prior(formula, data_treatment, model)
  pr_sum <- default_prior(formula, data_sum, model)
  pr_helmert <- default_prior(formula, data_helmert, model)
  
  # Class-level main prior should be consistent
  expect_equal(
    pr_treatment[pr_treatment$coef == "" & pr_treatment$class == "b", ]$prior,
    pr_sum[pr_sum$coef == "" & pr_sum$class == "b", ]$prior
  )
  expect_equal(
    pr_treatment[pr_treatment$coef == "" & pr_treatment$class == "b", ]$prior,
    pr_helmert[pr_helmert$coef == "" & pr_helmert$class == "b", ]$prior
  )
  expect_equal(
    pr_treatment[pr_treatment$coef == "" & pr_treatment$class == "b", ]$prior,
    c("normal(2, 1)", "logistic(0, 1)")
  )
  
  # Test Pattern 4: Multiple predictors, no intercept
  formula <- bmf(kappa ~ 0 + set_size + session, thetat ~ 0 + set_size + session)
  
  pr_treatment <- default_prior(formula, data_treatment, model)
  pr_sum <- default_prior(formula, data_sum, model)
  pr_helmert <- default_prior(formula, data_helmert, model)
  
  # Class-level effects prior should be consistent
  expect_equal(
    pr_treatment[pr_treatment$coef == "" & pr_treatment$class == "b", ]$prior,
    pr_sum[pr_sum$coef == "" & pr_sum$class == "b", ]$prior
  )
  expect_equal(
    pr_treatment[pr_treatment$coef == "" & pr_treatment$class == "b", ]$prior,
    pr_helmert[pr_helmert$coef == "" & pr_helmert$class == "b", ]$prior
  )
  expect_equal(
    pr_treatment[pr_treatment$coef == "" & pr_treatment$class == "b", ]$prior,
    c("normal(0, 1)", "normal(0, 0.5)")
  )
  
  # First predictor should have explicit main prior regardless of contrast coding
  expect_equal(
    pr_treatment[pr_treatment$coef == "set_size1" & pr_treatment$class == "b", ]$prior,
    pr_sum[pr_sum$coef == "set_size1" & pr_sum$class == "b", ]$prior
  )
  expect_equal(
    pr_treatment[pr_treatment$coef == "set_size1" & pr_treatment$class == "b", ]$prior,
    pr_helmert[pr_helmert$coef == "set_size1" & pr_helmert$class == "b", ]$prior
  )
  expect_equal(
    pr_treatment[pr_treatment$coef == "set_size1" & pr_treatment$class == "b", ]$prior,
    c("normal(2, 1)", "logistic(0, 1)")
  )
  
  # Test Pattern 5: Interactions with intercept
  formula <- bmf(kappa ~ 1 + set_size:session, thetat ~ 1 + set_size:session)
  
  pr_treatment <- default_prior(formula, data_treatment, model)
  pr_sum <- default_prior(formula, data_sum, model)
  pr_helmert <- default_prior(formula, data_helmert, model)
  
  # Class-level effects prior should be set consistently
  expect_equal(
    pr_treatment[pr_treatment$coef == "" & pr_treatment$class == "b", ]$prior,
    pr_sum[pr_sum$coef == "" & pr_sum$class == "b", ]$prior
  )
  expect_equal(
    pr_treatment[pr_treatment$coef == "" & pr_treatment$class == "b", ]$prior,
    pr_helmert[pr_helmert$coef == "" & pr_helmert$class == "b", ]$prior
  )
  expect_equal(
    pr_treatment[pr_treatment$coef == "" & pr_treatment$class == "b", ]$prior,
    c("normal(0, 1)", "normal(0, 0.5)")
  )
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
  expect_equal(pr[pr$coef == "" & pr$class == "b", ]$prior, c("normal(0, 1)", "normal(0, 0.5)"))
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
  expect_equal(pr[pr$coef == "" & pr$class == "b", ]$prior, c("normal(0, 1)", "normal(0, 0.5)"))
  expect_equal(unique(pr[pr$class == "sd", ]$prior), c("student_t(3, 0, 2.5)", ""))

  # 2 fixed effects + intercept suppressed
  formula <- bmf(kappa ~ 0 + set_size + session + (1 | ID), thetat ~ 0 + set_size + session + (1 | ID))
  pr <- default_prior(formula, data, model)
  expect_equal(pr[pr$coef == "Intercept" & pr$class == "b", ]$prior, character(0))
  expect_equal(pr[pr$coef == "" & pr$class == "b", ]$prior, c("normal(0, 1)", "normal(0, 0.5)"))
  expect_equal(pr[pr$coef == "set_size1" & pr$class == "b", ]$prior, c("normal(2, 1)", "logistic(0, 1)"))
  expect_equal(pr[pr$coef == "session2" & pr$class == "b", ]$prior, c("normal(0, 1)", "normal(0, 0.5)"))
  expect_equal(unique(pr[pr$class == "sd", ]$prior), c("student_t(3, 0, 2.5)", ""))

  # 2 fixed effects + interaction + intercept
  formula <- bmf(kappa ~ set_size * session + (1 | ID), thetat ~ set_size * session + (1 | ID))
  pr <- default_prior(formula, data, model)
  expect_equal(pr[pr$coef == "Intercept" & pr$class == "b", ]$prior, c("normal(2, 1)", "logistic(0, 1)"))
  expect_equal(pr[pr$coef == "" & pr$class == "b", ]$prior, c("normal(0, 1)", "normal(0, 0.5)"))
  expect_equal(unique(pr[pr$class == "sd", ]$prior), c("student_t(3, 0, 2.5)", ""))

  # 2 fixed effects + interaction + intercept suppressed
  formula <- bmf(kappa ~ 0 + set_size * session + (1 | ID), thetat ~ 0 + set_size * session + (1 | ID))
  pr <- default_prior(formula, data, model)
  expect_equal(pr[pr$coef == "Intercept" & pr$class == "b", ]$prior, character(0))
  expect_equal(pr[pr$coef == "" & pr$class == "b", ]$prior, c("normal(0, 1)", "normal(0, 0.5)"))
  expect_equal(pr[pr$coef == "set_size1" & pr$class == "b", ]$prior, c("normal(2, 1)", "logistic(0, 1)"))
  expect_equal(pr[pr$coef == "session2" & pr$class == "b", ]$prior, c("normal(0, 1)", "normal(0, 0.5)"))
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
  expect_equal(pr[pr$coef == "session2" & pr$class == "b", ]$prior, c("normal(0, 1)", "normal(0, 1)"))

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
  expect_equal(pr[pr$coef == "session2" & pr$class == "b", ]$prior, c("normal(0, 1)", "normal(0, 1)"))
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
