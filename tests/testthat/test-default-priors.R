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
  expect_equal(pr[pr$coef == "session2" & pr$class == "b", ]$prior, c("", ""))

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
  expect_equal(pr[pr$coef == "" & pr$class == "b", ]$prior, c("normal(0, 1)", "normal(0, 0.5)"))
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
  expect_equal(pr[pr$coef == "" & pr$class == "b", ]$prior, c("normal(0, 1)", "normal(0, 0.5)"))
  expect_equal(unique(pr[pr$class == "sd", ]$prior), c("", "exponential(1)"))

  # 2 fixed effects + intercept suppressed
  formula <- bmf(kappa ~ 0 + set_size + session + (1 | ID), thetat ~ 0 + set_size + session + (1 | ID))
  pr <- default_prior(formula, data, model)
  expect_equal(pr[pr$coef == "Intercept" & pr$class == "b", ]$prior, character(0))
  expect_equal(pr[pr$coef == "" & pr$class == "b", ]$prior, c("normal(0, 1)", "normal(0, 0.5)"))
  expect_equal(pr[pr$coef == "set_size1" & pr$class == "b", ]$prior, c("normal(2, 1)", "logistic(0, 1)"))
  expect_equal(pr[pr$coef == "session2" & pr$class == "b", ]$prior, c("", ""))
  expect_equal(unique(pr[pr$class == "sd", ]$prior), c("", "exponential(1)"))

  # 2 fixed effects + interaction + intercept
  formula <- bmf(kappa ~ set_size * session + (1 | ID), thetat ~ set_size * session + (1 | ID))
  pr <- default_prior(formula, data, model)
  expect_equal(pr[pr$coef == "Intercept" & pr$class == "b", ]$prior, c("normal(2, 1)", "logistic(0, 1)"))
  expect_equal(pr[pr$coef == "" & pr$class == "b", ]$prior, c("normal(0, 1)", "normal(0, 0.5)"))
  expect_equal(unique(pr[pr$class == "sd", ]$prior), c("", "exponential(1)"))

  # 2 fixed effects + interaction + intercept suppressed
  formula <- bmf(kappa ~ 0 + set_size * session + (1 | ID), thetat ~ 0 + set_size * session + (1 | ID))
  pr <- default_prior(formula, data, model)
  expect_equal(pr[pr$coef == "Intercept" & pr$class == "b", ]$prior, character(0))
  expect_equal(pr[pr$coef == "" & pr$class == "b", ]$prior, c("normal(0, 1)", "normal(0, 0.5)"))
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


test_that("mixture3p mixing weights get an effects prior on the softmax scale", {
  data <- oberauer_lin_2017
  model <- mixture3p("dev_rad", nt_features = paste0("col_nt", 1:7), set_size = "set_size")

  # set_size as predictor requires a suppressed intercept in this model
  formula <- bmf(kappa ~ 1, thetat ~ session, thetant ~ 0 + set_size)
  pr <- default_prior(formula, data, model)
  b_rows <- pr[pr$coef == "" & pr$class == "b", ]
  expect_equal(b_rows[b_rows$nlpar == "thetat", ]$prior, "normal(0, 0.5)")
  expect_equal(b_rows[b_rows$nlpar == "thetant", ]$prior, "logistic(0, 1)")
  expect_equal(pr[pr$coef == "set_size1" & pr$nlpar == "thetant", ]$prior, "constant(-100)")

  formula <- bmf(kappa ~ 1, thetat ~ 1, thetant ~ 0 + set_size + session)
  pr <- default_prior(formula, data, model)
  b_rows <- pr[pr$coef == "" & pr$class == "b", ]
  expect_equal(b_rows[b_rows$nlpar == "thetant", ]$prior, "normal(0, 0.5)")
  expect_equal(pr[pr$coef == "set_size1" & pr$nlpar == "thetant", ]$prior, "constant(-100)")
  expect_equal(pr[pr$coef == "set_size2" & pr$nlpar == "thetant", ]$prior, "logistic(0, 1)")
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
  model$default_priors$c$sd <- "exponential(3)"
  model$default_priors$kappa$sd <- NULL

  pr <- default_prior(bmf(c ~ 1 + (1 | ID), kappa ~ 1 + (1 | ID)), data, model)
  sd_pr <- sd_rows(pr)
  expect_equal(sd_pr[sd_pr$dpar == "c", ]$prior, "exponential(3)")
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

cor_default <- function(pr) {
  pr[pr$class == "cor" & pr$group == "", ]$prior
}

test_that("correlated random effects get an lkj(2) default on their correlations", {
  data <- oberauer_lin_2017
  model <- mixture2p("dev_rad")

  pr <- default_prior(bmf(kappa ~ set_size + (set_size | ID), thetat ~ 1), data, model)
  expect_equal(cor_default(pr), "lkj(2)")

  pr <- default_prior(bmf(kappa ~ 0 + set_size + (0 + set_size | ID), thetat ~ 1), data, model)
  expect_equal(cor_default(pr), "lkj(2)")

  # an ID shared across parameters builds one correlation matrix out of two intercepts
  pr <- default_prior(bmf(kappa ~ 1 + (1 | q | ID), thetat ~ 1 + (1 | q | ID)), data, model)
  expect_equal(cor_default(pr), "lkj(2)")
})

test_that("no cor prior is emitted for a model without a correlation matrix", {
  data <- oberauer_lin_2017
  model <- mixture2p("dev_rad")
  no_cor <- list(
    intercept_only = bmf(kappa ~ 1 + (1 | ID), thetat ~ 1),
    uncorrelated = bmf(kappa ~ set_size + (set_size || ID), thetat ~ 1),
    separate_parameters = bmf(kappa ~ 1 + (1 | ID), thetat ~ 1 + (1 | ID)),
    no_ranef = bmf(kappa ~ 1, thetat ~ 1)
  )

  for (formula in no_cor) {
    expect_false(any(default_prior(formula, data, model)$class == "cor"))
  }

  # brms rejects priors on non-existent parameters only when the model is built
  for (formula in no_cor[c("intercept_only", "uncorrelated")]) {
    expect_no_error(bmm(formula, data, model, backend = "mock", mock_fit = 1, rename = FALSE))
  }
})

test_that("the lkj default reaches the Stan code and yields to a user prior", {
  data <- oberauer_lin_2017
  model <- mixture2p("dev_rad")
  formula <- bmf(kappa ~ set_size + (set_size | ID), thetat ~ 1)

  fit <- bmm(formula, data, model, backend = "mock", mock_fit = 1, rename = FALSE)
  expect_match(brms::stancode(fit), "lkj_corr_cholesky_lpdf(L_1 | 2)", fixed = TRUE)

  fit <- bmm(formula, data, model,
    prior = brms::prior_("lkj(4)", class = "cor"),
    backend = "mock", mock_fit = 1, rename = FALSE
  )
  # a fitted object stores the row as class "L", so the Stan code is the evidence
  expect_match(brms::stancode(fit), "lkj_corr_cholesky_lpdf(L_1 | 4)", fixed = TRUE)
  expect_no_match(brms::stancode(fit), "lkj_corr_cholesky_lpdf(L_1 | 2)", fixed = TRUE)
})

test_that("a fit's stored correlation prior replaces the cor default on a refit", {
  data <- oberauer_lin_2017
  model <- mixture2p("dev_rad")
  formula <- bmf(kappa ~ set_size + (set_size | ID), thetat ~ 1)
  fit <- bmm(formula, data, model,
    prior = brms::prior_("lkj(4)", class = "cor"),
    backend = "mock", mock_fit = 1, rename = FALSE
  )

  # update() feeds the stored prior back as the user prior, and a fit stores the
  # row as class "L": next to a "cor" default brms sees a duplicated prior
  expect_true("L" %in% fit$prior$class)
  refit <- bmm(formula, data, model,
    prior = fit$prior,
    backend = "mock", mock_fit = 1, rename = FALSE
  )
  expect_match(brms::stancode(refit), "lkj_corr_cholesky_lpdf(L_1 | 4)", fixed = TRUE)
})

test_that("bmm.default_priors = FALSE also disables the cor default", {
  withr::local_options(bmm.default_priors = FALSE)
  pr <- default_prior(
    bmf(kappa ~ set_size + (set_size | ID), thetat ~ 1), oberauer_lin_2017, mixture2p("dev_rad")
  )
  expect_equal(cor_default(pr), "lkj(1)")
})

sd_default <- function(pr, par) {
  rows <- sd_rows(pr)
  rows[rows$nlpar == par | rows$dpar == par, ]$prior
}

test_that("every model ships an sd default on the link scale of each parameter", {
  data <- oberauer_lin_2017

  pr <- default_prior(
    bmf(kappa ~ 1 + (1 | ID), thetat ~ 1 + (1 | ID), thetant ~ 1 + (1 | ID)), data,
    mixture3p("dev_rad", nt_features = paste0("col_nt", 1:7), set_size = "set_size")
  )
  for (par in c("kappa", "thetat", "thetant")) expect_equal(sd_default(pr, par), "exponential(1)")

  pr <- default_prior(
    bmf(kappa ~ 1 + (1 | ID), a ~ 1 + (1 | ID), c ~ 1 + (1 | ID), s ~ 1 + (1 | ID)), data,
    imm("dev_rad", nt_features = paste0("col_nt", 1:7), nt_distances = paste0("dist_nt", 1:7), set_size = "set_size")
  )
  for (par in c("kappa", "a", "c", "s")) expect_equal(sd_default(pr, par), "exponential(1)")

  pr <- default_prior(bmf(kappa ~ 1 + (1 | ID), c ~ 1 + (1 | ID)), data, sdm("dev_rad"))
  for (par in c("kappa", "c")) expect_equal(sd_default(pr, par), "exponential(1)")

  m3_formula <- bmf(c ~ 1 + (1 | ID), a ~ 1 + (1 | ID))
  for (rule in c("simple", "softmax")) {
    model <- m3(resp_cats = c("corr", "other", "npl"), num_options = c("n_corr", "n_other", "n_npl"),
      choice_rule = rule, version = "ss")
    pr <- default_prior(m3_formula, oberauer_lewandowsky_2019_e1, model)
    for (par in c("c", "a")) expect_equal(sd_default(pr, par), "exponential(1)")
  }

  # the custom version derives its defaults from the link of each parameter
  withr::local_options(bmm.silent = 2)
  model <- m3(resp_cats = c("corr", "other", "npl"), num_options = c("n_corr", "n_other", "n_npl"),
    choice_rule = "softmax", version = "custom", links = list(c = "log", a = "logit"))
  pr <- suppressWarnings(default_prior(bmf(
    c ~ 1 + (1 | ID), a ~ 1 + (1 | ID),
    corr ~ b + a + c, other ~ b + a, npl ~ b
  ), oberauer_lewandowsky_2019_e1, model))
  for (par in c("c", "a")) expect_equal(sd_default(pr, par), "exponential(1)")

  rt_data <- data.frame(rt = rep(c(0.6, 0.8, 1.1, 0.7), 5), response = rep(c(1, 0), 10), id = factor(rep(1:5, each = 4)))
  rt_formula <- bmf(drift ~ 1 + (1 | id), bound ~ 1 + (1 | id), ndt ~ 1 + (1 | id))
  pr <- default_prior(rt_formula, rt_data, ddm("rt", "response"))
  expect_equal(sd_default(pr, "drift"), "exponential(1)")
  for (par in c("bound", "ndt")) expect_equal(sd_default(pr, par), "exponential(2)")

  # the toy data has a 50% error rate, which cswald "simple" warns about
  pr <- suppressWarnings(default_prior(rt_formula, rt_data, cswald("rt", "response", version = "simple")))
  for (par in c("drift", "bound", "ndt")) expect_equal(sd_default(pr, par), "exponential(2)")
  pr <- default_prior(rt_formula + bmf(zr ~ 1 + (1 | id)), rt_data, cswald("rt", "response", version = "crisk"))
  expect_equal(sd_default(pr, "drift"), "exponential(1)")
  for (par in c("bound", "ndt", "zr")) expect_equal(sd_default(pr, par), "exponential(2)")

  ez_data <- data.frame(
    mean_rt = rep(c(0.5, 0.6), 10), var_rt = rep(c(0.02, 0.03), 10),
    n_upper = rep(c(60, 70), 10), n_trials = 100, id = factor(rep(1:10, each = 2))
  )
  pr <- default_prior(rt_formula, ez_data, ezdm("mean_rt", "var_rt", "n_upper", "n_trials", version = "3par"))
  expect_equal(sd_default(pr, "drift"), "exponential(1)")
  for (par in c("bound", "ndt")) expect_equal(sd_default(pr, par), "exponential(2)")

  pr <- default_prior(
    bmf(d ~ 1 + (1 | id), criterion ~ 0 + condition + (1 | id), sdratio ~ 1 + (1 | id)),
    broeder_schuetz_2009_e3,
    sdt_yn(response = "n_old", stimulus = "stimulus", n_trials = "n_trials")
  )
  expect_equal(sd_default(pr, "d"), "exponential(1)")
  for (par in c("criterion", "sdratio")) expect_equal(sd_default(pr, par), "exponential(2)")

  mafc_data <- data.frame(
    n_correct = rep(c(30, 24), 10), n_trials = 40, id = factor(rep(1:10, each = 2))
  )
  pr <- default_prior(
    bmf(d ~ 1 + (1 | id)), mafc_data,
    sdt_mafc(response = "n_correct", n_trials = "n_trials", m = 4)
  )
  expect_equal(sd_default(pr, "d"), "exponential(1)")
})

test_that("a freed mu / mu1 gets regularizing main, effects and sd priors on the tan_half scale", {
  data <- oberauer_lin_2017

  pr <- default_prior(bmf(mu ~ 1 + set_size + (1 | ID), c ~ 1, kappa ~ 1), data, sdm("dev_rad"))
  expect_equal(pr[pr$class == "Intercept" & pr$dpar == "", ]$prior, "normal(0, 0.5)")
  expect_equal(pr[pr$class == "b" & pr$coef == "" & pr$dpar == "", ]$prior, "normal(0, 0.25)")
  # brms treats mu as its own parameter, so its sd row is the global one
  expect_equal(pr[pr$class == "sd" & pr$coef == "" & pr$group == "" & pr$dpar == "" & pr$nlpar == "", ]$prior, "exponential(4)")

  pr <- default_prior(bmf(mu1 ~ 0 + session + (1 | ID), kappa ~ 1, thetat ~ 1), data, mixture2p("dev_rad"))
  expect_equal(pr[pr$class == "b" & pr$coef == "" & pr$dpar == "mu1", ]$prior, "normal(0, 0.5)")
  expect_equal(sd_default(pr, "mu1"), "exponential(4)")
})

test_that("the set-size-1 sd constraint survives next to the blanket sd prior", {
  data <- oberauer_lin_2017
  model <- mixture3p("dev_rad", nt_features = paste0("col_nt", 1:7), set_size = "set_size")
  formula <- bmf(kappa ~ 1, thetat ~ 1, thetant ~ 0 + set_size + (0 + set_size | ID))

  pr <- default_prior(formula, data, model)
  constraint <- pr[pr$class == "sd" & pr$coef == "set_size1" & pr$nlpar == "thetant", ]
  expect_equal(constraint$prior, "constant(1e-8)")
  expect_equal(constraint$group, "ID")
  expect_equal(sd_default(pr, "thetant"), "exponential(1)")

  fit <- bmm(formula, data, model, backend = "mock", mock_fit = 1, rename = FALSE)
  code <- brms::stancode(fit)
  expect_match(code, "1e-08", fixed = TRUE)
  expect_match(code, "exponential_lpdf(sd_1", fixed = TRUE)

  model <- imm("dev_rad", nt_features = paste0("col_nt", 1:7), nt_distances = paste0("dist_nt", 1:7), set_size = "set_size")
  formula <- bmf(kappa ~ 1, c ~ 1, a ~ 0 + set_size + (0 + set_size | ID), s ~ 0 + set_size + (0 + set_size | ID))
  pr <- default_prior(formula, data, model)
  constraint <- pr[pr$class == "sd" & pr$coef == "set_size1", ]
  expect_setequal(constraint$nlpar, c("a", "s"))
  expect_true(all(constraint$prior == "constant(1e-8)"))
  for (par in c("a", "s")) expect_equal(sd_default(pr, par), "exponential(1)")
})

test_that("mock fits carry class-level b and sd priors for every formula shape and contrast coding", {
  withr::local_options(bmm.sort_data = FALSE)
  expect_defaults_applied <- function(fit, pars, re_pars) {
    pr <- fit$prior
    class_rows <- pr[pr$coef == "" & pr$group == "" & (pr$nlpar %in% pars | pr$dpar %in% pars), ]
    b_rows <- class_rows[class_rows$class == "b", ]
    expect_true(all(!b_rows$prior %in% c("", "(flat)")), label = "class-level b priors are not flat")
    sd_rows <- class_rows[class_rows$class == "sd", ]
    expect_setequal(c(sd_rows$nlpar, sd_rows$dpar)[c(sd_rows$nlpar, sd_rows$dpar) != ""], re_pars)
    expect_true(all(grepl("^exponential", sd_rows$prior)))
  }

  data <- oberauer_lin_2017
  model <- mixture2p("dev_rad")
  formulas <- list(
    bmf(kappa ~ set_size + (1 | ID), thetat ~ session + (1 | ID)),
    bmf(kappa ~ 0 + set_size + (0 + set_size || ID), thetat ~ 0 + set_size + session),
    bmf(kappa ~ set_size * session + (1 | ID), thetat ~ session)
  )
  re_pars <- list(c("kappa", "thetat"), "kappa", "kappa")
  for (contrasts in list(c("contr.treatment", "contr.poly"), c("contr.sum", "contr.poly"))) {
    withr::local_options(contrasts = contrasts)
    for (i in seq_along(formulas)) {
      fit <- bmm(formulas[[i]], data, model, backend = "mock", mock_fit = 1, rename = FALSE)
      expect_defaults_applied(fit, c("kappa", "thetat"), re_pars[[i]])
    }
  }

  fit <- bmm(bmf(kappa ~ set_size + (1 | ID), c ~ 0 + set_size + (0 + set_size || ID)), data, sdm("dev_rad"),
    backend = "mock", mock_fit = 1, rename = FALSE)
  expect_defaults_applied(fit, c("kappa", "c"), c("kappa", "c"))

  model <- m3(resp_cats = c("corr", "other", "npl"), num_options = c("n_corr", "n_other", "n_npl"),
    choice_rule = "softmax", version = "ss")
  fit <- bmm(bmf(c ~ 1 + cond + (1 + cond || ID), a ~ 1 + cond + (1 | ID)), oberauer_lewandowsky_2019_e1, model,
    backend = "mock", mock_fit = 1, rename = FALSE)
  expect_defaults_applied(fit, c("c", "a"), c("c", "a"))
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

test_that("sdt_yn emits an sd prior only for parameters with random effects", {
  data <- broeder_schuetz_2009_e3
  model <- sdt_yn(response = "n_old", stimulus = "stimulus", n_trials = "n_trials")

  # default_prior() alone does not validate a prior against the model, so the
  # no-random-effects case is only provable through a fit
  formula_fixed <- bmf(d ~ 1, criterion ~ 0 + condition, sdratio ~ 1)
  fit <- bmm(formula_fixed, data, model, backend = "mock", mock_fit = 1, rename = FALSE)
  expect_false(any(fit$prior$class == "sd"))
})
