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

test_that("SD priors are set for random intercept only", {
  data <- oberauer_lin_2017
  model <- mixture2p("dev_rad")

  formula <- bmf(kappa ~ 1 + (1 | ID), thetat ~ 1 + (1 | ID))
  pr <- default_prior(formula, data, model)
  sd_pr <- pr[pr$class == "sd" & pr$prior != "", ]

  # sd_main on intercept for kappa and thetat
  int_sd <- sd_pr[sd_pr$coef == "Intercept", ]
  expect_true(all(int_sd$prior == "exponential(1)"))
  expect_setequal(int_sd$nlpar, c("kappa", "thetat"))

  # sd_effects as blanket for kappa and thetat
  blanket_sd <- sd_pr[sd_pr$coef == "", ]
  expect_true(all(blanket_sd$prior == "exponential(1)"))
  expect_setequal(blanket_sd$nlpar, c("kappa", "thetat"))
})


test_that("SD priors are set for random intercept + slope", {
  data <- oberauer_lin_2017
  model <- mixture2p("dev_rad")

  formula <- bmf(kappa ~ set_size + (set_size | ID), thetat ~ 1 + (1 | ID))
  pr <- default_prior(formula, data, model)
  sd_pr <- pr[pr$class == "sd" & pr$prior != "", ]

  # kappa: sd_main on Intercept, sd_effects as blanket
  kappa_sd <- sd_pr[sd_pr$nlpar == "kappa", ]
  expect_true("Intercept" %in% kappa_sd$coef)
  expect_true("" %in% kappa_sd$coef)
  expect_true(all(kappa_sd$prior == "exponential(1)"))

  expect_true(all(sd_pr[sd_pr$nlpar == "thetat", "prior"] == "exponential(1)"))
})


test_that("SD priors are set for random slope only (no intercept)", {
  data <- oberauer_lin_2017
  model <- mixture2p("dev_rad")

  formula <- bmf(kappa ~ set_size + (0 + set_size | ID), thetat ~ 1)
  pr <- default_prior(formula, data, model)
  sd_pr <- pr[pr$class == "sd" & pr$prior != "", ]

  # sd_effects blanket present for kappa
  blanket <- sd_pr[sd_pr$nlpar == "kappa" & sd_pr$coef == "", ]
  expect_equal(nrow(blanket), 1)
  expect_equal(blanket$prior, "exponential(1)")

  expect_equal(nrow(sd_pr[sd_pr$nlpar == "kappa" & sd_pr$coef == "Intercept", ]), 0)
})


test_that("SD priors are set for crossed random effects", {
  data <- oberauer_lin_2017
  model <- mixture2p("dev_rad")

  formula <- bmf(kappa ~ 1 + (1 | ID) + (1 | session), thetat ~ 1 + (1 | ID))
  pr <- default_prior(formula, data, model)
  sd_pr <- pr[pr$class == "sd" & pr$prior != "", ]

  kappa_int <- sd_pr[sd_pr$nlpar == "kappa" & sd_pr$coef == "Intercept", ]
  expect_equal(nrow(kappa_int), 2)
  expect_setequal(kappa_int$group, c("ID", "session"))
  expect_true(all(kappa_int$prior == "exponential(1)"))

  thetat_int <- sd_pr[sd_pr$nlpar == "thetat" & sd_pr$coef == "Intercept", ]
  expect_equal(nrow(thetat_int), 1)
  expect_equal(thetat_int$group, "ID")
})


test_that("no SD priors without random effects", {
  data <- oberauer_lin_2017
  model <- mixture2p("dev_rad")

  formula <- bmf(kappa ~ 1, thetat ~ 1)
  pr <- default_prior(formula, data, model)

  expect_false(any(pr$class == "sd"))
})


test_that("SD priors work with SDM model (dpars)", {
  data <- oberauer_lin_2017
  model <- sdm("dev_rad")

  formula <- bmf(c ~ 1 + (1 | ID), kappa ~ 1 + (1 | ID))
  pr <- default_prior(formula, data, model)
  sd_pr <- pr[pr$class == "sd" & pr$prior != "", ]

  expect_true(all(sd_pr$prior == "exponential(1)"))
  expect_setequal(sd_pr$dpar, c("c", "kappa"))
})


test_that("user-specified SD priors override defaults", {
  data <- oberauer_lin_2017
  model <- mixture2p("dev_rad")

  formula <- bmf(kappa ~ 1 + (1 | ID), thetat ~ 1)
  user_prior <- brms::prior_("normal(0, 0.5)", class = "sd", nlpar = "kappa")
  pr <- default_prior(formula, data, model, prior = user_prior)
  sd_pr <- pr[pr$class == "sd" & pr$prior != "" & pr$nlpar == "kappa", ]

  expect_true("normal(0, 0.5)" %in% sd_pr$prior)
})


test_that("SD priors use parameter-specific rates for ezdm", {
  sim_data <- data.frame(
    mean_rt = rnorm(20, 0.5, 0.1), var_rt = runif(20, 0.01, 0.05),
    n_upper = sample(30:70, 20, TRUE), n_trials = rep(100, 20),
    id = rep(1:10, 2), cond = rep(c("A", "B"), each = 10)
  )
  model <- ezdm("mean_rt", "var_rt", "n_upper", "n_trials", version = "3par")

  formula <- bmf(
    drift ~ 1 + (1 | id),
    bound ~ 1 + (1 | id),
    ndt ~ 1 + (1 | id)
  )
  pr <- default_prior(formula, sim_data, model)
  sd_pr <- pr[pr$class == "sd" & pr$prior != "", ]

  expect_equal(sd_pr[sd_pr$dpar == "drift" & sd_pr$coef == "Intercept", "prior"], "exponential(1)")
  expect_equal(sd_pr[sd_pr$dpar == "drift" & sd_pr$coef == "", "prior"], "exponential(2)")

  expect_equal(sd_pr[sd_pr$dpar == "bound" & sd_pr$coef == "Intercept", "prior"], "exponential(2)")
  expect_equal(sd_pr[sd_pr$dpar == "bound" & sd_pr$coef == "", "prior"], "exponential(2)")

  expect_equal(sd_pr[sd_pr$dpar == "ndt" & sd_pr$coef == "Intercept", "prior"], "exponential(2)")
  expect_equal(sd_pr[sd_pr$dpar == "ndt" & sd_pr$coef == "", "prior"], "exponential(3)")
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
