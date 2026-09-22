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


# the init function bmm() builds: the prior decides which parameters exist
configured_initfun <- function(model, formula, data) {
  model <- check_model(model, data, formula)
  data <- check_data(model, data, formula)
  formula <- check_formula(model, data, formula)
  config_args <- configure_model(model, data, formula)
  prior <- configure_prior(model, data, config_args$formula, NULL)
  create_initfun(model, data, config_args$formula, prior)
}

stan_parameter_names <- function(model, formula, data) {
  names(extract_parameter_dimensions(
    extract_stan_blocks(stancode(formula, data, model))$parameters
  ))
}

expect_re_inits <- function(inits, correlated = TRUE) {
  sd_names <- grep("^sd_", names(inits), value = TRUE)
  expect_gt(length(sd_names), 0)
  for (nm in sd_names) {
    expect_true(all(inits[[nm]] >= 0.05 & inits[[nm]] <= 0.1))
  }
  z_names <- grep("^z_", names(inits), value = TRUE)
  expect_gt(length(z_names), 0)
  for (nm in z_names) {
    expect_true(all(abs(inits[[nm]]) <= 0.5))
  }
  L_names <- grep("^L_", names(inits), value = TRUE)
  if (!correlated) {
    expect_length(L_names, 0)
    return(invisible())
  }
  expect_gt(length(L_names), 0)
  for (nm in L_names) {
    expect_equal(inits[[nm]], diag(nrow = nrow(inits[[nm]])))
  }
}

test_that("create_initfun returns a function with random-effects inits for mixture2p", {
  dat <- oberauer_lin_2017
  model <- mixture2p(resp_error = "dev_rad")
  init_fun <- configured_initfun(model, bmf(thetat ~ 1, kappa ~ 1), dat)
  expect_true(is.function(init_fun))
  inits <- init_fun()
  expect_true(all(vapply(inits, function(x) all(is.finite(x)), logical(1))))

  formula <- bmf(
    kappa ~ 0 + set_size + (0 + set_size | ID),
    thetat ~ 0 + set_size + (0 + set_size | ID)
  )
  expect_re_inits(configured_initfun(model, formula, dat)())
})

test_that("create_initfun returns 0 for m3 with simple choice rule and identity link", {
  dat <- oberauer_lewandowsky_2019_e1
  ff <- bmf(c ~ 1, a ~ 1)

  model <- m3(
    resp_cats = c("corr", "other", "npl"),
    num_options = c("n_corr", "n_other", "n_npl"),
    choice_rule = "simple",
    version = "ss"
  )

  # default simple links are log -> falls through to the default method
  expect_true(is.function(configured_initfun(model, ff, dat)))

  # an identity link on any parameter requires zeros for stable sampling
  model$links$c <- "identity"
  expect_equal(configured_initfun(model, ff, dat), 0)
})

test_that("m3 with a softmax choice rule gets random-effects inits", {
  model <- m3(
    resp_cats = c("corr", "other", "npl"),
    num_options = c("n_corr", "n_other", "n_npl"),
    choice_rule = "softmax",
    version = "ss"
  )
  formula <- bmf(c ~ 1 + cond + (1 + cond | ID), a ~ 1 + cond + (1 + cond | ID))
  inits <- configured_initfun(model, formula, oberauer_lewandowsky_2019_e1)()
  expect_re_inits(inits)
  # m3 declares no init_ranges: its coefficients keep the radius-1 start
  expect_true(all(abs(c(inits$b_c, inits$b_a)) <= 1))
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

# =============================================================================
# SOFTPLUS LINK (opt-in alternative to log for positive parameters)
# =============================================================================

test_that("create_initfun handles a softplus link override on a positive parameter", {
  ff <- bmmformula(kappa ~ 1, c ~ 1)
  dat <- oberauer_lin_2017
  mod <- sdm(resp_error = "dev_rad", links = list(kappa = "softplus"))
  config_args <- configure_model(mod, data = dat, formula = ff)

  init_fun <- create_initfun(mod, dat, config_args$formula)
  inits <- init_fun()

  expect_type(inits, "list")
  expect_true(all(sapply(inits, function(x) all(is.finite(x)))))
})

# =============================================================================
# RANDOM-EFFECTS INITS FOR MODELS WITHOUT init_ranges BEFORE THIS CHANGE
# =============================================================================

test_that("imm and mixture3p get random-effects inits", {
  dat <- oberauer_lin_2017
  imm_model <- imm(
    resp_error = "dev_rad", nt_features = paste0("col_nt", 1:7),
    nt_distances = paste0("dist_nt", 1:7), set_size = "set_size"
  )
  imm_formula <- bmf(
    c ~ 0 + set_size + (0 + set_size | ID), a ~ 1 + (1 | ID),
    s ~ 1, kappa ~ 1
  )
  expect_re_inits(configured_initfun(imm_model, imm_formula, dat)())

  mix3p_model <- mixture3p("dev_rad", nt_features = paste0("col_nt", 1:7), set_size = "set_size")
  intercepts <- bmf(kappa ~ 1 + (1 | ID), thetat ~ 1 + (1 | ID), thetant ~ 1 + (1 | ID))
  expect_re_inits(configured_initfun(mix3p_model, intercepts, dat)(), correlated = FALSE)
  slopes <- bmf(
    kappa ~ 0 + set_size + (0 + set_size | ID), thetat ~ 1,
    thetant ~ 0 + set_size + (0 + set_size | ID)
  )
  expect_re_inits(configured_initfun(mix3p_model, slopes, dat)())
})

test_that("the init list names every parameter of the Stan model and nothing else", {
  dat <- oberauer_lin_2017
  cases <- list(
    list(
      model = mixture2p("dev_rad"),
      formula = bmf(
        kappa ~ 0 + set_size + (0 + set_size | ID),
        thetat ~ 0 + set_size + (0 + set_size | ID)
      )
    ),
    # set size 1 pins thetant and its sd to constants, which brms implements by
    # replacing the vectors b_thetant and sd_3 with per-coefficient scalars
    list(
      model = mixture3p("dev_rad", nt_features = paste0("col_nt", 1:7), set_size = "set_size"),
      formula = bmf(
        kappa ~ 0 + set_size + (0 + set_size | ID),
        thetat ~ 0 + set_size + (0 + set_size | ID),
        thetant ~ 0 + set_size + (0 + set_size | ID)
      )
    )
  )
  for (case in cases) {
    spars <- stan_parameter_names(case$model, case$formula, dat)
    inits <- configured_initfun(case$model, case$formula, dat)()
    expect_setequal(names(inits), spars)
  }
  expect_true(any(grepl("^par_b_thetant_", spars)))
  par_sd <- unlist(inits[grep("^par_sd_", names(inits))])
  expect_gt(length(par_sd), 0)
  expect_true(all(par_sd >= 0.05 & par_sd <= 0.1))
  par_b <- unlist(inits[grep("^par_b_thetant_", names(inits))])
  expect_true(all(abs(par_b) <= 1.1))
  sdata <- standata(cases[[2]]$formula, dat, cases[[2]]$model)
  expect_length(inits$b_kappa, sdata$K_kappa)
  expect_equal(dim(inits$z_1), c(sdata$M_1, sdata$N_1))
})

test_that("uncorrelated random effects get z inits and no correlation matrix", {
  dat <- oberauer_lin_2017
  model <- mixture2p("dev_rad")
  formula <- bmf(kappa ~ 0 + set_size + (0 + set_size || ID), thetat ~ 1)
  inits <- configured_initfun(model, formula, dat)()
  sdata <- standata(formula, dat, model)
  expect_re_inits(inits, correlated = FALSE)
  expect_equal(dim(inits$z_1), c(sdata$M_1, sdata$N_1))
})

test_that("a user init replaces the package init", {
  dat <- oberauer_lin_2017
  model <- mixture2p("dev_rad")
  formula <- bmf(kappa ~ 1 + (1 | ID), thetat ~ 1)
  fit <- bmm(formula, dat, model, backend = "mock", mock_fit = 1, rename = FALSE)
  expect_true(is.function(fit$stan_args$init))
  fit <- bmm(formula, dat, model, init = 0, backend = "mock", mock_fit = 1, rename = FALSE)
  expect_equal(fit$stan_args$init, 0)
})

# =============================================================================
# POPULATION-LEVEL RANGES OF THE CIRCULAR MIXTURE MODELS
# =============================================================================

# central 50% of the main default prior, mapped to the native scale
central_range <- function(location, scale, inverse) {
  signif(inverse(location + c(-1, 1) * scale * stats::qnorm(0.75)), 2)
}

test_that("mixture2p starts inside the central 50% of its default priors", {
  dat <- oberauer_lin_2017
  model <- mixture2p("dev_rad")
  ranges <- model$init_ranges
  expect_equal(ranges$kappa, central_range(2, 1, exp))
  expect_equal(ranges$mu1, central_range(0, 0.5, function(x) 2 * atan(x)))
  expect_equal(ranges$thetat, signif(stats::plogis(c(-1, 1) * stats::qlogis(0.75)), 2))

  formula <- bmf(kappa ~ 1 + (1 | ID), thetat ~ 1 + set_size, mu1 ~ 1)
  init_fun <- configured_initfun(model, formula, dat)
  for (i in 1:20) {
    inits <- init_fun()
    # nlpars carry their intercept in the first coefficient
    expect_true(inits$b_kappa[1] >= log(ranges$kappa[1]) && inits$b_kappa[1] <= log(ranges$kappa[2]))
    expect_true(inits$b_thetat[1] >= stats::qlogis(ranges$thetat[1]) &&
      inits$b_thetat[1] <= stats::qlogis(ranges$thetat[2]))
    expect_true(all(abs(inits$b_thetat[-1]) <= 0.1))
    expect_true(abs(inits$Intercept_mu1) <= tan(ranges$mu1[2] / 2))
  }
})

test_that("mixture3p softmax weights start inside their range on the sampling scale", {
  model <- mixture3p("dev_rad", nt_features = paste0("col_nt", 1:7), set_size = "set_size")
  expect_equal(model$init_ranges$thetat, c(-1.1, 1.1))
  expect_equal(model$init_ranges$thetant, c(-1.1, 1.1))
  formula <- bmf(kappa ~ 1, thetat ~ 1, thetant ~ 1)
  init_fun <- configured_initfun(model, formula, oberauer_lin_2017)
  for (i in 1:20) {
    inits <- init_fun()
    expect_true(abs(inits$b_thetat[1]) <= 1.1)
    expect_true(abs(inits$b_thetant[1]) <= 1.1)
  }
})

test_that("imm versions carry ranges for exactly their parameters", {
  ranges <- imm(
    resp_error = "dev_rad", nt_features = paste0("col_nt", 1:7),
    nt_distances = paste0("dist_nt", 1:7), set_size = "set_size"
  )$init_ranges
  expect_setequal(names(ranges), c("mu1", "kappa", "a", "c", "s"))
  expect_equal(ranges$a, central_range(0, 1, exp))
  expect_equal(ranges$c, ranges$a)
  expect_equal(ranges$s, ranges$a)
  abc <- imm(
    resp_error = "dev_rad", nt_features = paste0("col_nt", 1:7),
    set_size = "set_size", version = "abc"
  )$init_ranges
  expect_setequal(names(abc), c("mu1", "kappa", "a", "c"))
  bsc <- imm(
    resp_error = "dev_rad", nt_features = paste0("col_nt", 1:7),
    nt_distances = paste0("dist_nt", 1:7), set_size = "set_size", version = "bsc"
  )$init_ranges
  expect_setequal(names(bsc), c("mu1", "kappa", "c", "s"))
})

# =============================================================================
# RADIUS-1 FILL FOR PARAMETERS WITHOUT A RANGE
# =============================================================================

test_that("init_default_param() maps the draw through the declared bounds", {
  spec <- function(type, bounds, dims = "K", types = type) nlist(type, types, bounds, dims)
  unbounded <- init_default_param(spec("vector", NULL), 20, list())
  lower <- init_default_param(spec("vector", list(lower = "2")), 20, list())
  upper <- init_default_param(spec("vector", list(upper = "-3")), 20, list())
  both <- init_default_param(
    spec("vector", list(lower = "0", upper = "min_Y")), 20, list(min_Y = 4)
  )

  expect_true(all(abs(unbounded) <= 1))
  expect_true(all(lower > 2))
  expect_true(all(upper < -3))
  expect_true(all(both > 0 & both < 4))
  expect_equal(dim(init_default_param(spec("matrix", NULL), c(2, 3), list())), c(2, 3))
  # a size-1 vector must stay an array so that it is written as [x], not x
  expect_equal(dim(init_default_param(spec("vector", NULL), 1, list())), 1)
  expect_null(dim(init_default_param(spec("real", NULL, dims = "1"), 1, list())))
  # the same holds for an array of reals, which the parser gives a trailing 1
  real_array <- spec("real", NULL, dims = c("K", "1"), types = c("array", "real"))
  expect_equal(dim(init_default_param(real_array, c(1, 1), list())), 1)
})

test_that("init_default_param() leaves unsupported declarations to the sampler", {
  spec <- function(type, bounds) nlist(type, types = type, bounds, dims = "K")
  expect_null(init_default_param(spec("simplex", NULL), 3, list()))
  expect_null(init_default_param(spec("real", list(lower = "unknown_var")), 1, list()))
})

test_that("Stan dimensions and bounds resolve from literals and data", {
  sdata <- list(K = 3L, N = 10L, lb = 0.5)
  expect_equal(resolve_stan_dim(c("K", "N"), sdata), c(3, 10))
  expect_equal(resolve_stan_dim("1", sdata), 1)
  expect_equal(resolve_stan_bound(NULL, sdata, -Inf), -Inf)
  expect_equal(resolve_stan_bound("2.5", sdata, -Inf), 2.5)
  expect_equal(resolve_stan_bound("lb", sdata, -Inf), 0.5)
  expect_true(is.na(resolve_stan_bound("missing", sdata, -Inf)))
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

test_that("the intercept of a non-linear parameter starts in range (#362)", {
  # native-multinomial / non-linear models (e.g. sdt_rating) carry their
  # parameters as nlpars, whose intercept brms folds into the first coefficient
  # of the b_ vector instead of declaring an Intercept_ parameter
  model <- list(
    parameters = list(eta = ""), init_ranges = list(eta = c(2, 3)), links = list(eta = "identity")
  )
  X <- stats::model.matrix(~ 1, data.frame(y = 1:3))
  inits <- init_fixef_param("b_eta", "vector", 1, model, list(X_eta = X), stan_names = "b_eta")
  expect_equal(dim(inits), 1)
  expect_true(inits >= 2 && inits <= 3)
})

test_that("a per-coefficient scalar starts where its position in the vector would", {
  model <- list(
    parameters = list(eta = ""), init_ranges = list(eta = c(2, 3)), links = list(eta = "identity")
  )
  X <- stats::model.matrix(~ 1 + x, data.frame(x = factor(1:3)))
  sdata <- list(X_eta = X)
  intercept <- init_fixef_param("par_b_eta_1", "real", 1, model, sdata, stan_names = "par_b_eta_1")
  effect <- init_fixef_param("par_b_eta_3", "real", 1, model, sdata, stan_names = "par_b_eta_3")
  expect_null(dim(intercept))
  expect_true(intercept >= 2 && intercept <= 3)
  expect_true(abs(effect) <= 0.1)
})

# =============================================================================
# SHAPES AND RANGES ACROSS FORMULA STRUCTURES
# =============================================================================

# rstan reads a real as a bare number and everything else as an array with the
# declared dimensions; a matrix type declares its size once
expect_stan_shapes <- function(model, formula, data) {
  pars <- extract_parameter_dimensions(
    extract_stan_blocks(stancode(formula, data, model))$parameters
  )
  sdata <- standata(formula, data, model)
  inits <- configured_initfun(model, formula, data)()
  expect_setequal(names(inits), names(pars))
  for (nm in names(inits)) {
    dims <- resolve_stan_dim(pars[[nm]]$dims, sdata)
    expected <- switch(pars[[nm]]$type,
      real = dims[-length(dims)],
      cholesky_factor_corr = c(dims, dims[length(dims)]),
      dims
    )
    if (length(expected) == 0) {
      expect_null(dim(inits[[nm]]), label = nm)
    } else {
      expect_equal(dim(inits[[nm]]), expected, label = nm)
    }
  }
}

test_that("every init has the shape its Stan declaration asks for", {
  dat <- oberauer_lin_2017
  dat$grp <- factor(as.integer(dat$ID) %% 2)
  mix3p <- mixture3p("dev_rad", nt_features = paste0("col_nt", 1:7), set_size = "set_size")
  # size-1 coefficient vectors, which rstan cannot read from a bare number
  expect_stan_shapes(mixture2p("dev_rad"), bmf(kappa ~ 1, thetat ~ 1 + set_size), dat)
  # per-coefficient scalars and uncorrelated random effects
  expect_stan_shapes(mix3p, bmf(
    kappa ~ 1, thetat ~ 0 + set_size + (0 + set_size || ID), thetant ~ 0 + set_size
  ), dat)
  # the main dpar's coefficients and a Stan intercept
  expect_stan_shapes(sdm("dev_rad"), bmf(mu ~ 0 + set_size, c ~ 1 + set_size, kappa ~ 1), dat)
  # one correlation matrix per level of a gr(by = ) variable
  expect_stan_shapes(mixture2p("dev_rad"), bmf(
    kappa ~ 1 + (1 + set_size | gr(ID, by = grp)), thetat ~ 1
  ), dat)
})

test_that("the first term of a no-intercept formula starts in range whatever its form", {
  dat <- oberauer_lin_2017
  dat$ss_num <- as.numeric(as.character(dat$set_size))
  model <- mixture2p("dev_rad")
  log_range <- log(model$init_ranges$kappa)
  n_in_range <- function(formula, data = dat) {
    b <- configured_initfun(model, formula, data)()$b_kappa
    in_range <- b >= log_range[1] & b <= log_range[2]
    expect_true(all(in_range | abs(b) <= 0.1))
    c(length(b), sum(in_range))
  }
  expect_equal(n_in_range(bmf(kappa ~ 0 + factor(ss_num), thetat ~ 1)), c(8, 8))
  expect_equal(n_in_range(bmf(kappa ~ 0 + set_size + session, thetat ~ 1)), c(9, 8))
  expect_equal(n_in_range(bmf(kappa ~ 0 + set_size:session, thetat ~ 1)), c(16, 16))
  expect_equal(n_in_range(bmf(kappa ~ 0 + ss_num + session, thetat ~ 1)), c(3, 1))
  # brms drops the unused level, and so does the count
  without_8 <- dat[dat$set_size != 8, ]
  expect_equal(n_in_range(bmf(kappa ~ 0 + set_size, thetat ~ 1), without_8), c(7, 7))
})

test_that("the coefficients of the main dpar start from its range", {
  dat <- oberauer_lin_2017
  model <- sdm("dev_rad")
  upper <- tan(model$init_ranges$mu[2] / 2)
  cells <- configured_initfun(model, bmf(mu ~ 0 + set_size, c ~ 1, kappa ~ 1), dat)()
  expect_length(cells$b, 8)
  expect_true(all(abs(cells$b) <= upper))
  expect_gt(max(abs(cells$b)), 0.1)
  effects <- configured_initfun(model, bmf(mu ~ 1 + set_size, c ~ 1, kappa ~ 1), dat)()
  expect_true(all(abs(effects$b) <= 0.1))
  expect_true(abs(effects$Intercept) <= upper)
})
