test_that("create_initfun returns function for sdm", {
  # prepare info for tests
  ff <- bmmformula(kappa ~ 1, c ~ 1)
  dat <- oberauer_lin_2017
  mod <- sdm(resp_error = "dev_rad")
  config_args <- configure_model(mod, data = dat, formula = ff)

  # create initfun
  init_fun <- create_initfun(mod, dat, config_args$formula)

  # run tests
  expect_equal(class(init_fun), "function")
  expect_equal(class(unlist(init_fun())), "numeric")
})


test_that("create_initfun returns 1 for mixture2p models", {
  # prepare info for tests
  dat <- oberauer_lin_2017
  model_mix2p <- mixture2p(resp_error = "dev_rad")

  ff_mix2p <- bmf(thetat ~ 1, kappa ~ 1)

  config_args_mix2p <- configure_model(model_mix2p, data = dat, formula = ff_mix2p)

  # create initfun
  init_fun <- create_initfun(model_mix2p, dat, config_args_mix2p$formula)

  # run tests
  expect_equal(class(init_fun), "numeric")
  expect_equal(init_fun, 1)
})

# =============================================================================
# BASIC FUNCTIONALITY TESTS
# =============================================================================

test_that("initfun generates valid numeric initial values", {
  ff <- bmmformula(kappa ~ 1, c ~ 1)
  dat <- oberauer_lin_2017

  mod <- sdm(resp_error = "dev_rad")
  config_args <- configure_model(mod, data = dat, formula = ff)

  init_fun <- create_initfun(mod, dat, config_args$formula)
  inits <- init_fun()

  expect_type(inits, "list")
  expect_true(all(sapply(inits, function(x) is.numeric(x) || is.matrix(x) || is.array(x))))
  expect_true(all(sapply(inits, function(x) all(is.finite(x)))))
})

# =============================================================================
# INTERCEPT-ONLY MODELS (real type parameters)
# =============================================================================

test_that("initfun generates correct intercept values for sdm", {
  ff <- bmmformula(kappa ~ 1, c ~ 1)
  dat <- oberauer_lin_2017
  mod <- sdm(resp_error = "dev_rad")
  config_args <- configure_model(mod, data = dat, formula = ff)

  init_fun <- create_initfun(mod, dat, config_args$formula)
  inits <- init_fun()

  # Check that Intercept parameters exist
  intercept_names <- grep("Intercept", names(inits), value = TRUE)
  expect_true(length(intercept_names) > 0)

  # Check values are scalars (length 1)
  for (nm in intercept_names) {
    expect_equal(length(inits[[nm]]), 1)
  }
})

# =============================================================================
# MODELS WITH PREDICTOR EFFECTS (vector parameters)
# =============================================================================

test_that("initfun handles single predictor without intercept", {
  dat <- oberauer_lin_2017
  dat$condition <- factor(rep(c("A", "B"), length.out = nrow(dat)))

  ff <- bmmformula(kappa ~ 0 + condition, c ~ 1)
  mod <- sdm(resp_error = "dev_rad")
  config_args <- configure_model(mod, data = dat, formula = ff)

  init_fun <- create_initfun(mod, dat, config_args$formula)
  inits <- init_fun()

  # b_kappa should have 2 values (one per level)
  b_kappa <- inits[["b_kappa"]]
  expect_equal(length(b_kappa), 2)
  expect_true(all(is.finite(b_kappa)))
})

test_that("initfun handles predictor with intercept", {
  dat <- oberauer_lin_2017
  dat$condition <- factor(rep(c("A", "B"), length.out = nrow(dat)))

  ff <- bmmformula(kappa ~ 1 + condition, c ~ 1)
  mod <- sdm(resp_error = "dev_rad")
  config_args <- configure_model(mod, data = dat, formula = ff)

  init_fun <- create_initfun(mod, dat, config_args$formula)
  inits <- init_fun()

  # Should have intercept + effect coded predictor
  expect_true("Intercept_kappa" %in% names(inits) || any(grepl("b_kappa", names(inits))))
})

test_that("initfun handles multiple predictors", {
  dat <- oberauer_lin_2017
  dat$cond1 <- factor(rep(c("A", "B"), length.out = nrow(dat)))
  dat$cond2 <- factor(rep(c("X", "Y", "Z"), length.out = nrow(dat)))

  ff <- bmmformula(kappa ~ 0 + cond1 + cond2, c ~ 1)
  mod <- sdm(resp_error = "dev_rad")
  config_args <- configure_model(mod, data = dat, formula = ff)

  init_fun <- create_initfun(mod, dat, config_args$formula)
  inits <- init_fun()

  # b_kappa should have values for first term transformed, rest small
  b_kappa <- inits[["b_kappa"]]
  expect_true(length(b_kappa) >= 2)
  expect_true(all(is.finite(b_kappa)))
})

test_that("initfun handles interaction terms", {
  dat <- oberauer_lin_2017
  dat$cond1 <- factor(rep(c("A", "B"), length.out = nrow(dat)))
  dat$cond2 <- factor(rep(c("X", "Y"), length.out = nrow(dat)))

  ff <- bmmformula(kappa ~ 0 + cond1:cond2, c ~ 1)
  mod <- sdm(resp_error = "dev_rad")
  config_args <- configure_model(mod, data = dat, formula = ff)

  init_fun <- create_initfun(mod, dat, config_args$formula)
  inits <- init_fun()

  # Should handle interaction term correctly
  b_kappa <- inits[["b_kappa"]]
  expect_equal(length(b_kappa), 4)
  expect_true(all(is.finite(b_kappa)))
})

test_that("initfun handles interaction terms with other terms", {
  dat <- oberauer_lin_2017
  dat$cond1 <- factor(rep(c("A", "B"), length.out = nrow(dat)))
  dat$cond2 <- factor(rep(c("X", "Y"), length.out = nrow(dat)))
  dat$cond3 <- factor(rep(c("S", "T"), length.out = nrow(dat)))

  ff <- bmmformula(kappa ~ 0 + cond1:cond2 + cond1:cond3, c ~ 1)
  mod <- sdm(resp_error = "dev_rad")
  config_args <- configure_model(mod, data = dat, formula = ff)

  init_fun <- create_initfun(mod, dat, config_args$formula)
  inits <- init_fun()

  # Should handle interaction term correctly
  b_kappa <- inits[["b_kappa"]]
  expect_equal(length(b_kappa), 6)
  expect_true(all(is.finite(b_kappa)))
})

# =============================================================================
# RANDOM EFFECTS TESTS
# =============================================================================

test_that("initfun generates sd parameters for random effects", {
  dat <- oberauer_lin_2017

  ff <- bmmformula(kappa ~ 1 + (1 | ID), c ~ 1)
  mod <- sdm(resp_error = "dev_rad")
  config_args <- configure_model(mod, data = dat, formula = ff)

  init_fun <- create_initfun(mod, dat, config_args$formula)
  inits <- init_fun()

  # Should have sd_ parameters
  sd_names <- grep("^sd_", names(inits), value = TRUE)
  expect_true(length(sd_names) > 0)

  # sd parameters should be positive and small
  for (nm in sd_names) {
    expect_true(all(inits[[nm]] > 0))
    expect_true(all(inits[[nm]] < 1))
  }
})

test_that("initfun generates z values for random effects", {
  dat <- oberauer_lin_2017

  ff <- bmmformula(kappa ~ 1 + (1 | ID), c ~ 1)
  mod <- sdm(resp_error = "dev_rad")
  config_args <- configure_model(mod, data = dat, formula = ff)

  init_fun <- create_initfun(mod, dat, config_args$formula)
  inits <- init_fun()

  # Should have z_ parameters (arrays)
  z_names <- grep("^z_", names(inits), value = TRUE)
  expect_true(length(z_names) > 0)

  # z values should be small (around 0)
  for (nm in z_names) {
    expect_true(all(abs(inits[[nm]]) <= 0.5))
  }
})

test_that("initfun handles correlated random effects", {
  dat <- oberauer_lin_2017
  dat$condition <- factor(rep(c("A", "B"), length.out = nrow(dat)))

  ff <- bmmformula(kappa ~ 1 + condition + (1 + condition | ID), c ~ 1)
  mod <- sdm(resp_error = "dev_rad")
  config_args <- configure_model(mod, data = dat, formula = ff)

  init_fun <- create_initfun(mod, dat, config_args$formula)
  inits <- init_fun()

  # Should have correlation matrix (L_ or cor_ parameters)
  cor_names <- grep("^(L_|cor_)", names(inits), value = TRUE)
  expect_true(length(cor_names) > 0)
})

# =============================================================================
# LINK FUNCTION TESTS
# =============================================================================

test_that("initfun applies log link correctly for kappa", {
  ff <- bmmformula(kappa ~ 1, c ~ 1)
  dat <- oberauer_lin_2017
  mod <- sdm(resp_error = "dev_rad")
  config_args <- configure_model(mod, data = dat, formula = ff)

  init_fun <- create_initfun(mod, dat, config_args$formula)
  inits <- init_fun()

  # kappa intercept should be on log scale
  # init_ranges for kappa are c(2.5, 3.5), log transformed should be ~ log(2.5) to log(3.5)
  kappa_int <- inits[["Intercept_kappa"]]
  expect_true(kappa_int > log(2) && kappa_int < log(4))
})

test_that("initfun handles NULL/missing links as identity", {
  dat <- oberauer_lin_2017
  mod <- sdm(resp_error = "dev_rad")

  # Manually remove a link to simulate NULL case
  mod$links$kappa <- NULL

  ff <- bmmformula(kappa ~ 1, c ~ 1)
  config_args <- configure_model(mod, data = dat, formula = ff)

  # This should not error due to our fix
  init_fun <- create_initfun(mod, dat, config_args$formula)
  expect_true(is.function(init_fun))

  inits <- init_fun()
  expect_true(is.list(inits))
  expect_true(all(sapply(inits, function(x) all(is.finite(x)))))
})

# =============================================================================
# REPRODUCIBILITY AND RANDOMNESS TESTS
# =============================================================================

test_that("initfun generates different values on repeated calls", {
  ff <- bmmformula(kappa ~ 1, c ~ 1)
  dat <- oberauer_lin_2017
  mod <- sdm(resp_error = "dev_rad")
  config_args <- configure_model(mod, data = dat, formula = ff)

  init_fun <- create_initfun(mod, dat, config_args$formula)

  inits1 <- init_fun()
  inits2 <- init_fun()

  # At least one parameter should differ (randomness)
  all_equal <- all(mapply(function(a, b) identical(a, b), inits1, inits2))
  expect_false(all_equal)
})

test_that("initfun values are within expected ranges", {
  ff <- bmmformula(kappa ~ 1, c ~ 1)
  dat <- oberauer_lin_2017
  mod <- sdm(resp_error = "dev_rad")
  config_args <- configure_model(mod, data = dat, formula = ff)

  init_fun <- create_initfun(mod, dat, config_args$formula)

  # Run multiple times to check consistency
  for (i in 1:10) {
    inits <- init_fun()

    # All values should be finite
    expect_true(all(sapply(inits, function(x) all(is.finite(x)))))

    # No extreme values
    numeric_vals <- unlist(lapply(inits, as.numeric))
    expect_true(all(abs(numeric_vals) < 100))
  }
})

# =============================================================================
# EDGE CASES
# =============================================================================

test_that("initfun handles single random effect group correctly", {
  dat <- oberauer_lin_2017

  # Use a formula that results in single sd parameter per group
  ff <- bmmformula(kappa ~ 1 + (1 | ID), c ~ 1)
  mod <- sdm(resp_error = "dev_rad")
  config_args <- configure_model(mod, data = dat, formula = ff)

  init_fun <- create_initfun(mod, dat, config_args$formula)
  inits <- init_fun()

  # sd parameters should be properly formatted even when single-dimensional
  sd_names <- grep("^sd_", names(inits), value = TRUE)
  for (nm in sd_names) {
    expect_true(is.numeric(inits[[nm]]) || is.array(inits[[nm]]))
    expect_true(all(is.finite(inits[[nm]])))
  }
})

test_that("initfun handles numeric predictors", {
  dat <- oberauer_lin_2017
  dat$continuous_pred <- rnorm(nrow(dat))

  ff <- bmmformula(kappa ~ 1 + continuous_pred, c ~ 1)
  mod <- sdm(resp_error = "dev_rad")
  config_args <- configure_model(mod, data = dat, formula = ff)

  init_fun <- create_initfun(mod, dat, config_args$formula)
  inits <- init_fun()

  expect_true(is.list(inits))
  expect_true(all(sapply(inits, function(x) all(is.finite(x)))))
})

# =============================================================================
# STRUCTURE VALIDATION
# =============================================================================

test_that("initfun output matches standata dimensions", {
  # Use a model with predictors to ensure b_ parameters exist
  dat <- oberauer_lin_2017
  
  ff <- bmmformula(kappa ~ 1 + set_size, c ~ 1)
  mod <- sdm(resp_error = "dev_rad")
  config_args <- configure_model(mod, data = dat, formula = ff)

  init_fun <- create_initfun(mod, dat, config_args$formula)
  inits <- init_fun()

  standata <- brms::standata(config_args$formula, dat, config_args$formula$family)

  # Verify that we have b_ parameters to test
  b_names <- grep("^b_", names(inits), value = TRUE)
  expect_true(length(b_names) > 0, info = "Should have at least one b_ parameter")

  # Verify dimensions match for b_ (non-intercept) parameters
  # Note: b_ parameters correspond to Kc_ (centered, excluding intercept) in standata
  for (nm in b_names) {
    param <- sub("^b_", "", nm)
    # For models with intercepts, brms uses Kc_ for centered predictors
    dim_name_c <- paste0("Kc_", param)
    # For models without intercepts, brms uses K_
    dim_name <- paste0("K_", param)
    
    if (dim_name_c %in% names(standata)) {
      expect_equal(
        length(inits[[nm]]), 
        standata[[dim_name_c]],
        info = paste("Dimension mismatch for parameter:", nm)
      )
    } else if (dim_name %in% names(standata)) {
      expect_equal(
        length(inits[[nm]]), 
        standata[[dim_name]],
        info = paste("Dimension mismatch for parameter:", nm)
      )
    }
  }
})

test_that("initfun output matches standata dimensions for no-intercept models", {
  # Use a model without intercept to test K_ dimension matching
  dat <- oberauer_lin_2017
  
  ff <- bmmformula(kappa ~ 0 + set_size, c ~ 1)
  mod <- sdm(resp_error = "dev_rad")
  config_args <- configure_model(mod, data = dat, formula = ff)

  init_fun <- create_initfun(mod, dat, config_args$formula)
  inits <- init_fun()

  standata <- brms::standata(config_args$formula, dat, config_args$formula$family)

  # Verify that we have b_ parameters to test
  b_names <- grep("^b_", names(inits), value = TRUE)
  expect_true(length(b_names) > 0, info = "Should have at least one b_ parameter")

  # For models without intercepts, kappa should NOT have an Intercept_kappa parameter
  expect_false("Intercept_kappa" %in% names(inits), 
               info = "No-intercept model should not have Intercept_kappa parameter")

  # Verify dimensions match using K_ (not Kc_) for no-intercept models
  for (nm in b_names) {
    param <- sub("^b_", "", nm)
    dim_name <- paste0("K_", param)
    
    expect_true(dim_name %in% names(standata), 
                info = paste("K_ dimension should exist for no-intercept model:", dim_name))
    
    expect_equal(
      length(inits[[nm]]),
      standata[[dim_name]],
      info = paste("Dimension mismatch for no-intercept parameter:", nm)
    )
  }
})


# -----------------------------------------------------------------------------
# match_stan_to_model_par tests (substring collision scenarios)
# -----------------------------------------------------------------------------

test_that("match_stan_to_model_par handles exact matches", {
  model_pars <- c("kappa", "c", "mu")
  expect_equal(match_stan_to_model_par("Intercept_kappa", model_pars), "kappa")
  expect_equal(match_stan_to_model_par("Intercept_c", model_pars), "c")
  expect_equal(match_stan_to_model_par("Intercept", model_pars), "mu")
})

test_that("match_stan_to_model_par avoids substring collisions", {
  model_pars <- c("s", "sim", "ndt", "bound")
  expect_equal(match_stan_to_model_par("Intercept_sim", model_pars), "sim")
  expect_equal(match_stan_to_model_par("Intercept_s", model_pars), "s")
  expect_equal(match_stan_to_model_par("b_sim", model_pars), "sim")
  expect_equal(match_stan_to_model_par("sd_1", model_pars), "sd")
})

test_that("match_stan_to_model_par avoids single-letter collisions", {
  model_pars <- c("c", "correct", "a", "activation")
  expect_equal(match_stan_to_model_par("Intercept_correct", model_pars), "correct")
  expect_equal(match_stan_to_model_par("Intercept_c", model_pars), "c")
  expect_equal(match_stan_to_model_par("Intercept_activation", model_pars), "activation")
  expect_equal(match_stan_to_model_par("Intercept_a", model_pars), "a")
})

test_that("match_stan_to_model_par prefers longest match for prefix params", {
  model_pars <- c("mu", "mu1", "mu2", "kappa", "kappa2")
  expect_equal(match_stan_to_model_par("Intercept_mu1", model_pars), "mu1")
  expect_equal(match_stan_to_model_par("Intercept_mu2", model_pars), "mu2")
  expect_equal(match_stan_to_model_par("Intercept_kappa2", model_pars), "kappa2")
  expect_equal(match_stan_to_model_par("Intercept_kappa", model_pars), "kappa")
})

test_that("match_stan_to_model_par falls back for structural params", {
  model_pars <- c("kappa", "c")
  expect_equal(match_stan_to_model_par("sd_1", model_pars), "sd")
  expect_equal(match_stan_to_model_par("z_1", model_pars), "z")
  expect_equal(match_stan_to_model_par("cor_1", model_pars), "cor")
})


# -----------------------------------------------------------------------------
# nlpar parameter resolution (#362)
# -----------------------------------------------------------------------------

test_that("init terms resolve via nlpars for non-linear models (#362)", {
  # native-multinomial / non-linear models (e.g. sdt_rating) carry their
  # parameters in bterms$nlpars rather than bterms$dpars; an intercept-only
  # nlpar is represented as a `b_` vector and so routes through init_vector_param
  bterms <- brms::brmsterms(brms::bf(y ~ eta, eta ~ 1, nl = TRUE))
  expect_null(bterms$dpars[["eta"]])
  expect_false(is.null(bterms$nlpars[["eta"]]))

  par_terms <- bterms$dpars[["eta"]] %||% bterms$nlpars[["eta"]]
  expect_identical(par_terms, bterms$nlpars[["eta"]])

  dat <- data.frame(y = rnorm(10))
  inits <- init_vector_param(
    "b_eta", dim = 1, init_range = c(0, 1), link = "identity",
    bterms = par_terms, data = dat
  )
  expect_length(inits, 1)
  expect_true(all(is.finite(inits)))

  # the pre-fix dpars-only lookup yields NULL and errors in has_intercept()
  expect_error(init_vector_param(
    "b_eta", dim = 1, init_range = c(0, 1), link = "identity",
    bterms = bterms$dpars[["eta"]], data = dat
  ))
})
