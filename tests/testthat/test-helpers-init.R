skip_on_cran() # cross-model contract file, skipped on CRAN as a whole (ADR 0004, R10)

test_that("create_initfun returns function for sdm", {
  # prepare info for tests
  ff <- bmmformula(kappa ~ 1, c ~ 1)
  dat <- small_data()
  mod <- sdm(resp_error = "dev_rad")
  config_args <- configure_model(mod, data = dat, formula = ff)

  # create initfun
  init_fun <- create_initfun(mod, dat, config_args$formula)

  # run tests
  expect_equal(class(init_fun), "function")
  expect_equal(class(unlist(init_fun())), "numeric")
})


# the init function bmm() builds: the prior decides which parameters exist
configured_initfun <- function(model, formula, data, user_prior = NULL) {
  model <- check_model(model, data, formula)
  data <- check_data(model, data, formula)
  formula <- check_formula(model, data, formula)
  config_args <- configure_model(model, data, formula)
  prior <- configure_prior(model, data, config_args$formula, user_prior)
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
  dat <- small_data()
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

  # the init ranges are on the native scale and reach Stan on the link scale
  inits <- init_fun()
  expect_true(exp(inits$Intercept_kappa) >= 3.8 && exp(inits$Intercept_kappa) <= 15)
  expect_true(
    stats::plogis(inits$Intercept_thetat) >= 0.25 &&
      stats::plogis(inits$Intercept_thetat) <= 0.75
  )
})

test_that("create_initfun keeps the variable-precision gamma shape integrable", {
  dat <- oberauer_lin_2017
  model <- mixture2p(resp_error = "dev_rad", variable_precision = TRUE)
  init_fun <- configured_initfun(model, bmf(thetat ~ 1, kappa ~ 1, tau ~ 1), dat)

  shapes <- replicate(50, {
    inits <- init_fun()
    .circmix_J(exp(inits$Intercept_kappa)) / exp(inits$Intercept_tau)
  })
  expect_true(all(shapes > 60 / model$vp_nodes))
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
  expect_warning(
    init <- configured_initfun(model, ff, dat),
    "allow values that the model's default"
  )
  expect_equal(init, 0)
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
  dat <- small_data()

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
  dat <- small_data()
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
  dat <- small_data()
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
  dat <- small_data()
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
  dat <- small_data()
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
  dat <- small_data()
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
  dat <- small_data()
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
  dat <- small_data()

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
  dat <- small_data()

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
  dat <- small_data()
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
  dat <- small_data()
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
  dat <- small_data()
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
  dat <- small_data()
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
  dat <- small_data()
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
  dat <- small_data()

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
  withr::local_seed(445)
  dat <- small_data()
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
  dat <- small_data()
  
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
  dat <- small_data()
  
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
  ff <- bmmformula(drift ~ 1, bound ~ 1, ndt ~ 1)
  dat <- data.frame(
    mean_rt = c(0.50, 0.60, 0.70, 0.55),
    var_rt = c(0.020, 0.030, 0.025, 0.020),
    n_upper = c(80, 85, 75, 78),
    n_trials = 100
  )
  mod <- ezdm("mean_rt", "var_rt", "n_upper", "n_trials",
    links = list(bound = "softplus")
  )
  dat <- check_data(mod, dat, ff)
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
  dat <- small_data()
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
  # set size 1 has no non-target, which mixture3p warns about for thetant
  expect_re_inits(suppressWarnings(configured_initfun(mix3p_model, slopes, dat))())
})

test_that("the init list names every parameter of the Stan model and nothing else", {
  dat <- small_data()
  cases <- list(
    list(
      model = mixture2p("dev_rad"),
      formula = bmf(
        kappa ~ 0 + set_size + (0 + set_size | ID),
        thetat ~ 0 + set_size + (0 + set_size | ID)
      )
    ),
    list(
      model = mixture3p("dev_rad",
        nt_features = paste0("col_nt", 1:7), set_size = "set_size",
        version = "slot_averaging", variable_precision = TRUE
      ),
      formula = bmf(
        kappa ~ 1 + (1 | ID), tau ~ 1, K ~ 1 + (1 | ID), ptarget ~ 1 + (1 | ID)
      )
    ),
    # nothing pins the set-size-1 level of thetant any more, so b_thetant and its
    # sd stay vectors rather than per-coefficient scalars
    list(
      model = mixture3p("dev_rad", nt_features = paste0("col_nt", 1:7), set_size = "set_size"),
      formula = bmf(
        kappa ~ 0 + set_size + (0 + set_size | ID),
        thetat ~ 0 + set_size + (0 + set_size | ID),
        thetant ~ 0 + set_size + (0 + set_size | ID)
      )
    )
  )
  # set size 1 has no non-target, which mixture3p warns about for thetant
  for (case in cases) suppressWarnings({
    spars <- stan_parameter_names(case$model, case$formula, dat)
    inits <- configured_initfun(case$model, case$formula, dat)()
    expect_setequal(names(inits), spars)
  })
  expect_false(any(grepl("^par_b_thetant_", spars)))
  expect_length(inits$b_thetant, 8)
  # the softmax weights take their init range on the sampling scale
  expect_true(all(abs(inits$b_thetant) <= 1.1))
  sdata <- suppressWarnings(standata(cases[[3]]$formula, dat, cases[[3]]$model))
  expect_length(inits$b_kappa, sdata$K_kappa)
  expect_equal(dim(inits$z_1), c(sdata$M_1, sdata$N_1))
})

test_that("uncorrelated random effects get z inits and no correlation matrix", {
  dat <- small_data()
  model <- mixture2p("dev_rad")
  formula <- bmf(kappa ~ 0 + set_size + (0 + set_size || ID), thetat ~ 1)
  inits <- configured_initfun(model, formula, dat)()
  sdata <- standata(formula, dat, model)
  expect_re_inits(inits, correlated = FALSE)
  expect_equal(dim(inits$z_1), c(sdata$M_1, sdata$N_1))
})

test_that("a user init replaces the package init", {
  dat <- small_data()
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
  dat <- small_data()
  model <- mixture2p("dev_rad")
  ranges <- model$init_ranges
  expect_equal(ranges$kappa, central_range(2, 1, exp))
  expect_equal(ranges$mu, central_range(0, 0.5, function(x) 2 * atan(x)))
  expect_equal(ranges$thetat, signif(stats::plogis(c(-1, 1) * stats::qlogis(0.75)), 2))

  formula <- bmf(kappa ~ 1 + (1 | ID), thetat ~ 1 + set_size, mu ~ 1)
  init_fun <- configured_initfun(model, formula, dat)
  for (i in 1:20) {
    inits <- init_fun()
    # dpars carry their intercept separately from the coefficients
    expect_true(inits$Intercept_kappa >= log(ranges$kappa[1]) &&
      inits$Intercept_kappa <= log(ranges$kappa[2]))
    expect_true(inits$Intercept_thetat >= stats::qlogis(ranges$thetat[1]) &&
      inits$Intercept_thetat <= stats::qlogis(ranges$thetat[2]))
    expect_true(all(abs(inits$b_thetat) <= 0.1))
    # brms names the intercept of the custom family's first parameter, mu, plainly
    expect_true(abs(inits$Intercept) <= tan(ranges$mu[2] / 2))
  }
})

test_that("mixture3p softmax weights start inside their range on the sampling scale", {
  model <- mixture3p("dev_rad", nt_features = paste0("col_nt", 1:7), set_size = "set_size")
  expect_equal(model$init_ranges$thetat, c(-1.1, 1.1))
  expect_equal(model$init_ranges$thetant, c(-1.1, 1.1))
  formula <- bmf(kappa ~ 1, thetat ~ 1, thetant ~ 1)
  init_fun <- configured_initfun(model, formula, small_data())
  for (i in 1:20) {
    inits <- init_fun()
    expect_true(abs(inits$Intercept_thetat) <= 1.1)
    expect_true(abs(inits$Intercept_thetant) <= 1.1)
  }
})

test_that("imm versions carry ranges for exactly their parameters", {
  ranges <- imm(
    resp_error = "dev_rad", nt_features = paste0("col_nt", 1:7),
    nt_distances = paste0("dist_nt", 1:7), set_size = "set_size"
  )$init_ranges
  expect_setequal(names(ranges), c("mu", "kappa", "a", "c", "s", "b"))
  expect_equal(ranges$a, central_range(0, 1, exp))
  expect_equal(ranges$c, ranges$a)
  expect_equal(ranges$s, ranges$a)
  abc <- imm(
    resp_error = "dev_rad", nt_features = paste0("col_nt", 1:7),
    set_size = "set_size", version = "abc"
  )$init_ranges
  expect_setequal(names(abc), c("mu", "kappa", "a", "c", "b"))
  bsc <- imm(
    resp_error = "dev_rad", nt_features = paste0("col_nt", 1:7),
    nt_distances = paste0("dist_nt", 1:7), set_size = "set_size", version = "bsc"
  )$init_ranges
  expect_setequal(names(bsc), c("mu", "kappa", "c", "s", "b"))
})

# =============================================================================
# RADIUS-1 FILL FOR PARAMETERS WITHOUT A RANGE
# =============================================================================

test_that("init_default_param() maps the draw through the declared bounds", {
  spec <- function(type, dims = "K", types = type) nlist(type, types, dims)
  unbounded <- init_default_param(spec("vector"), 20, c(-Inf, Inf))
  lower <- init_default_param(spec("vector"), 20, c(2, Inf))
  upper <- init_default_param(spec("vector"), 20, c(-Inf, -3))
  both <- init_default_param(spec("vector"), 20, c(0, 4))

  expect_true(all(abs(unbounded) <= 1))
  expect_true(all(lower > 2))
  expect_true(all(upper < -3))
  expect_true(all(both > 0 & both < 4))
  expect_equal(dim(init_default_param(spec("matrix"), c(2, 3), c(-Inf, Inf))), c(2, 3))
  # a size-1 vector must stay an array so that it is written as [x], not x
  expect_equal(dim(init_default_param(spec("vector"), 1, c(-Inf, Inf))), 1)
  expect_null(dim(init_default_param(spec("real", dims = "1"), 1, c(-Inf, Inf))))
  # the same holds for an array of reals, which the parser gives a trailing 1
  real_array <- spec("real", dims = c("K", "1"), types = c("array", "real"))
  expect_equal(dim(init_default_param(real_array, c(1, 1), c(-Inf, Inf))), 1)
})

test_that("init_default_param() leaves unsupported declarations to the sampler", {
  spec <- list(type = "ordered", types = "ordered", dims = "K")
  expect_null(init_default_param(spec, 3, c(-Inf, Inf)))
})

test_that("a simplex starts from Stan's stick-breaking draw", {
  spec <- list(type = "simplex", types = "simplex", bounds = NULL, dims = "K")
  for (i in 1:20) {
    s <- init_default_param(spec, 4, c(-Inf, Inf))
    expect_equal(dim(s), 4)
    expect_true(all(s > 0))
    expect_equal(sum(s), 1)
  }
  # a zero draw on every break gives equal shares
  expect_equal(init_simplex(rep(0, 4), 5), rep(0.2, 5))
})

test_that("a parameter whose size does not resolve is left to the sampler", {
  model <- list(parameters = list(kappa = ""), init_ranges = list(), links = list())
  spec <- list(type = "vector", types = "vector", bounds = NULL, dims = "missing[1]")
  expect_null(init_stan_param("zs_kappa_1_1", spec, model, list()))
  spec$dims <- "K"
  expect_length(init_stan_param("zs_kappa_1_1", spec, model, list(K = 3L)), 3)
})

test_that("a parameter whose bounds do not resolve is left to the sampler", {
  withr::local_seed(694)
  model <- list(parameters = list(kappa = ""), init_ranges = list(kappa = c(1, 2)), links = list(kappa = "log"))
  sdata <- list(K = 2L, X_kappa = matrix(1, 3, 2), min_Y = 4, Y = c(1, 4, 2))
  sigma <- list(type = "real", types = "real", bounds = list(lower = "0", upper = "Y"), dims = "1")
  b <- list(type = "vector", types = "vector", bounds = list(lower = "unknown_var"), dims = "K")
  expect_null(init_stan_param("sigma", sigma, model, sdata))
  expect_null(init_stan_param("b_kappa", b, model, sdata))
  sigma$bounds$upper <- "min_Y"
  sigmas <- replicate(50, init_stan_param("sigma", sigma, model, sdata))
  expect_true(all(sigmas > 0 & sigmas < 4))
})

test_that("a declaration the parser cannot type is left to the sampler", {
  model <- list(
    parameters = list(kappa = ""), init_ranges = list(kappa = c(1, 2)),
    links = list(kappa = "log")
  )
  # a name the population-level rule would otherwise claim
  spec <- parse_parameters_line("sum_to_zero_vector[K] b_kappa;")
  expect_null(init_stan_param("b_kappa", spec, model, list(K = 3L, X_kappa = matrix(1, 3, 3))))
})

test_that("Stan dimensions and bounds resolve from literals and data", {
  sdata <- list(K = 3L, N = 10L, lb = 0.5)
  expect_equal(resolve_stan_dim(c("K", "N"), sdata), c(3, 10))
  expect_equal(resolve_stan_dim("1", sdata), 1)
  # mo() and s() size their parameters by an element of a data array
  expect_equal(resolve_stan_dim(c("Jmo[1]", "Jmo[2]"), list(Jmo = c(4L, 2L))), c(4, 2))
  expect_true(is.na(resolve_stan_dim("missing[1]", sdata)))
  expect_true(is.na(resolve_stan_dim("Jmo[i]", list(Jmo = c(4L, 2L)))))
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
  inits <- init_fixef_param("b_eta", "vector", 1, model, list(X_eta = X), c(-Inf, Inf))
  expect_equal(dim(inits), 1)
  expect_true(inits >= 2 && inits <= 3)
})

test_that("a per-coefficient scalar starts where its position in the vector would", {
  model <- list(
    parameters = list(eta = ""), init_ranges = list(eta = c(2, 3)), links = list(eta = "identity")
  )
  X <- stats::model.matrix(~ 1 + x, data.frame(x = factor(1:3)))
  sdata <- list(X_eta = X)
  intercept <- init_fixef_param("par_b_eta_1", "real", 1, model, sdata, c(-Inf, Inf))
  effect <- init_fixef_param("par_b_eta_3", "real", 1, model, sdata, c(-Inf, Inf))
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
  dat <- small_data()
  dat$grp <- factor(as.integer(dat$ID) %% 2)
  mix3p <- mixture3p("dev_rad", nt_features = paste0("col_nt", 1:7), set_size = "set_size")
  # size-1 coefficient vectors, which rstan cannot read from a bare number
  expect_stan_shapes(mixture2p("dev_rad"), bmf(kappa ~ 1, thetat ~ 1 + set_size), dat)
  # uncorrelated random effects; set size 1 has no non-target, which
  # mixture3p warns about for thetant
  suppressWarnings(expect_stan_shapes(mix3p, bmf(
    kappa ~ 1, thetat ~ 0 + set_size + (0 + set_size || ID), thetant ~ 0 + set_size
  ), dat))
  # the main dpar's coefficients and a Stan intercept
  expect_stan_shapes(sdm("dev_rad"), bmf(mu ~ 0 + set_size, c ~ 1 + set_size, kappa ~ 1), dat)
  # one correlation matrix per level of a gr(by = ) variable
  expect_stan_shapes(mixture2p("dev_rad"), bmf(
    kappa ~ 1 + (1 + set_size | gr(ID, by = grp)), thetat ~ 1
  ), dat)
})

# one test per case: each builds the Stan code and data, and together they
# would take most of the per-test time budget
for (name in names(model_test_cases("init"))) {
  test_that(paste0("the init function of ", name, " gives each Stan parameter its declared shape"), {
    withr::local_options(bmm.silent = 2)
    case <- model_test_cases("init")[[name]]
    # the toy RT data have a 50% error rate, which cswald "simple" warns about
    suppressWarnings(expect_stan_shapes(case$model, case$formula, case$data))
  })
}

test_that("mo() and s() terms get a complete init list", {
  # their parameters are sized by an element of a data array, which the parser
  # used to split at the inner bracket; cmdstanr then reported a partial list
  dat <- small_data()
  dat$x <- as.numeric(as.character(dat$set_size))
  dat$xo <- factor(pmin(dat$x, 4), ordered = TRUE)
  expect_stan_shapes(sdm("dev_rad"), bmf(c ~ 1 + mo(xo), kappa ~ 1), dat)
  expect_stan_shapes(mixture2p("dev_rad"), bmf(kappa ~ 1 + s(x, k = 4), thetat ~ 1), dat)
})

test_that("the first term of a no-intercept formula starts in range whatever its form", {
  dat <- small_data()
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
  # the numeric term carries no constant, so the session cell means hold it
  expect_equal(n_in_range(bmf(kappa ~ 0 + ss_num + session, thetat ~ 1)), c(3, 2))
  # brms drops the unused level, and so does the count
  without_8 <- dat[dat$set_size != 8, ]
  expect_equal(n_in_range(bmf(kappa ~ 0 + set_size, thetat ~ 1), without_8), c(7, 7))
})

test_that("range_coefficients() puts the linear predictor at the target on every row", {
  target <- 1.7
  dat <- data.frame(f = factor(rep(1:3, 4)), x = rep(1:3, 4))
  # brms holds the intercept of a centred design in a parameter of its own, so
  # the centred columns have to leave whatever value it draws alone
  predictor <- function(X, centered) {
    b <- range_coefficients(X, target, centered = centered)
    if (!centered) {
      return(as.vector(X %*% b))
    }
    as.vector(target + scale(X[, -1, drop = FALSE], scale = FALSE) %*% b)
  }

  # the order of the terms must not matter
  expect_equal(predictor(stats::model.matrix(~ 0 + f + x, dat), FALSE), rep(target, 12))
  expect_equal(predictor(stats::model.matrix(~ 0 + x + f, dat), FALSE), rep(target, 12))
  expect_equal(predictor(stats::model.matrix(~ 0 + f:x, dat), FALSE), rep(target, 12))
  expect_equal(predictor(stats::model.matrix(~ 1 + f, dat), TRUE), rep(target, 12))

  # f:x above is a relabelled cell-means design because f and x are perfectly
  # aligned in dat; genuinely crossed factors cannot reach a constant, the same
  # orthogonality property as the poly() case below
  crossed <- expand.grid(f = factor(1:3), x = 1:4)
  X <- stats::model.matrix(~ 0 + f:x, crossed)
  b <- range_coefficients(X, target, centered = FALSE)
  expect_false(isTRUE(all.equal(as.vector(X %*% b), rep(target, nrow(crossed)))))
  expect_equal(as.vector(t(X) %*% (X %*% b - target)), rep(0, ncol(X)))

  # a 0 + Intercept design matrix has no assign attribute
  X <- stats::model.matrix(~ 1 + f, dat)
  attr(X, "assign") <- NULL
  expect_equal(predictor(X, FALSE), rep(target, 12))

  # a design that cannot reach a constant gets the closest fit: the residual is
  # orthogonal to every column of the design
  X <- stats::model.matrix(~ 0 + poly(x, 2), dat)
  b <- range_coefficients(X, target, centered = FALSE)
  expect_false(isTRUE(all.equal(as.vector(X %*% b), rep(target, 12))))
  expect_equal(as.vector(t(X) %*% (X %*% b - target)), c(0, 0))

  # an aliased column solves to NA, which would poison every value Stan reads
  X <- cbind(Intercept = 1, x = dat$x, copy = dat$x)
  b <- range_coefficients(X, target, centered = FALSE)
  expect_false(anyNA(b))
  expect_equal(as.vector(X %*% b), rep(target, 12))

  # log(0) in a predictor reaches brms's design matrix, and qr() rejects it; the
  # coefficients start at zero rather than taking bmm() down before Stan sees it
  X <- cbind(Intercept = 1, x = c(1, 2, -Inf, 4))
  expect_equal(range_coefficients(X, target, centered = FALSE), c(0, 0))
  expect_equal(range_coefficients(X, target, centered = TRUE), 0)
})

test_that("the coefficients of the main dpar start from its range", {
  withr::local_seed(20260922)
  dat <- small_data()
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

test_that("a declaration whose size the parser cannot resolve is left to the sampler", {
  dat <- small_data()
  dat$x <- as.numeric(as.character(dat$set_size))
  init_fun <- configured_initfun(mixture2p("dev_rad"), bmf(kappa ~ 1 + s(x, k = 4), thetat ~ 1), dat)
  expect_no_error(init_fun())
})

test_that("a 0 + Intercept formula gets finite initial values", {
  # brms's X_<par> carries no assign attribute for this term, and the intercept
  # is not centred, so only its own column starts in range
  dat <- small_data()
  dat$x <- as.numeric(as.character(dat$set_size))
  model <- mixture2p("dev_rad")
  b <- configured_initfun(model, bmf(kappa ~ 0 + Intercept + x, thetat ~ 1), dat)()$b_kappa
  expect_true(all(is.finite(b)))
  expect_length(b, 2)
  expect_true(b[1] >= log(model$init_ranges$kappa[1]) && b[1] <= log(model$init_ranges$kappa[2]))
  expect_true(abs(b[2]) <= 0.1)
})

test_that("a constant intercept leaves every coefficient of a centred parameter an effect", {
  # Stan has no Intercept_c, but X_c still holds the intercept column and b_c
  # has one entry fewer: none of them may take the intercept's range
  dat <- small_data()
  prior <- brms::set_prior("constant(1)", class = "Intercept", dpar = "c")
  fit <- bmm(bmf(c ~ 1 + set_size, kappa ~ 1), dat, sdm("dev_rad"),
    prior = prior, backend = "mock", mock_fit = 1, rename = FALSE
  )
  b <- fit$stan_args$init()$b_c
  expect_length(b, 7)
  expect_true(all(abs(b) <= 0.1))
})

# =============================================================================
# INIT RANGES THE LINK CAN ONLY PARTLY REPRESENT (#460)
# =============================================================================

# criterion's range straddles 0, which softplus cannot represent
softplus_criterion <- list(
  parameters = list(criterion = ""), init_ranges = list(criterion = c(-0.5, 0.5)),
  links = list(criterion = "softplus")
)

test_that("an intercept starts inside the part of its range the link can represent", {
  withr::local_seed(460)
  inits <- unlist(lapply(1:200, function(i) {
    init_fixef_param("Intercept_criterion", "real", 1, softplus_criterion, list(), c(-Inf, Inf))
  }))
  expect_true(is.numeric(inits) && all(is.finite(inits)))
  native <- link_transform(inits, "softplus", inverse = TRUE)
  expect_true(all(native > 0 & native <= 0.5))
})

test_that("coefficients start inside the part of the range the link can represent", {
  withr::local_seed(460)
  X <- stats::model.matrix(~ 1 + x, data.frame(x = factor(c(1, 2, 2))))
  expect_no_warning(
    b <- lapply(1:200, function(i) {
      init_fixef_param("b_criterion", "vector", 2, softplus_criterion, list(X_criterion = X), c(-Inf, Inf))
    })
  )
  native <- link_transform(unlist(lapply(b, function(x) X %*% x)), "softplus", inverse = TRUE)
  expect_true(all(native > 0 & native <= 0.5))
})

test_that("a range the link cannot represent at all is left to the default draw", {
  withr::local_seed(460)
  model <- list(
    parameters = list(bound = ""), init_ranges = list(bound = c(1.5, 2)),
    links = list(bound = "logit")
  )
  spec <- list(type = "real", types = "real", bounds = NULL, dims = "1")
  expect_no_warning(
    inits <- vapply(1:200, function(i) init_stan_param("Intercept_bound", spec, model, list()), numeric(1))
  )
  expect_true(all(abs(inits) <= 1))
})

test_that("sdt_yn with a softplus criterion gets a finite init list", {
  withr::local_seed(460)
  dat <- data.frame(y = c(12, 15, 5, 7), stimulus = c(1, 1, 0, 0), n = 20)
  model <- sdt_yn("y", "stimulus", "n", links = list(criterion = "softplus"))
  init_fun <- configured_initfun(model, bmf(d ~ 1, criterion ~ 1), dat)
  expect_no_warning(inits <- lapply(1:50, function(i) init_fun()))
  expect_true(all(vapply(inits, function(x) "Intercept_criterion" %in% names(x), logical(1))))
  expect_true(all(is.finite(unlist(inits))))
})

test_that("a softmax range stays on the sampling scale and a single point is not drawn", {
  # a softmax weight's range is not a (0, 1) probability, so it is not clipped
  model <- list(
    parameters = list(thetat = ""), init_ranges = list(thetat = c(-2, -1)),
    links = list(thetat = "softmax")
  )
  init <- init_fixef_param("Intercept_thetat", "real", 1, model, list(), c(-Inf, Inf))
  expect_true(init >= -2 && init <= -1)
  # log(0) is -Inf: the one point log cannot represent leaves nothing to draw from
  expect_null(representable_range(c(-1, 0), "log"))
})

test_that("an init function carries nothing from the frame that built it", {
  withr::local_seed(445)
  dat <- small_data()
  model <- sdm(resp_error = "dev_rad")
  formula <- bmf(c ~ 1, kappa ~ 1)
  # brms stores the init function in the fit, and bmm()'s frame holds the data
  # and the fit itself, which saveRDS() would then write into the file again
  init_beside <- function(n) {
    ballast <- runif(n)
    configure_fit(formula, dat, model)$config_args$init
  }
  # each function is sized on its own, so that the other one's ballast is not counted
  serialized_size <- function(n) length(serialize(init_beside(n), NULL))
  expect_equal(serialized_size(1e6), serialized_size(1))
  inits <- init_beside(1e6)()
  expect_type(inits, "list")
  expect_named(inits)
})

test_that("an init function carries no data column the model does not use", {
  model <- sdm(resp_error = "dev_rad")
  formula <- bmf(c ~ 1, kappa ~ 1)
  serialized_size <- function(dat) {
    length(serialize(configure_fit(formula, dat, model)$config_args$init, NULL))
  }
  # check_data() keeps the 18 columns this model does not read, and a copy of
  # the data in the closure would make every saved fit grow with them
  expect_equal(serialized_size(oberauer_lin_2017), serialized_size(oberauer_lin_2017["dev_rad"]))
})

# =============================================================================
# PARAMETERS THE DATA BOUND (#621)
# =============================================================================

# 80 trials of two participants whose fastest response is 0.03 s, so that the
# ndt of cswald starts from c(0.012, 0.018): log(0.018 / 0.012) = 0.4055 wide on
# the log link, 0.006 on the identity link. trial runs to 40
bounded_ndt_data <- function() {
  data <- expand.grid(trial = 1:40, id = factor(1:2))
  data$rt <- 0.03 + rep(c(0, 0.2, 0.45, 0.7), length.out = 80)
  data$response <- rep(c(1, 0, 1, 1, 0), length.out = 80)
  data$setsize <- factor(rep(c(2, 4, 6), length.out = 80), ordered = TRUE)
  data$id2 <- factor(rep(c(2, 1, 1, 2, 2), length.out = 80))
  data$g <- factor(ifelse(data$id == 1, "a", "b"))
  data
}

bounded_ndt_inits <- function(formula, links = NULL, n = 50) {
  # check_data() warns about responses below 0.1 s and about the identity link
  init_fun <- suppressWarnings(
    configured_initfun(cswald("rt", "response", links = links), formula, bounded_ndt_data())
  )
  replicate(n, init_fun(), simplify = FALSE)
}

draws_of <- function(inits, name, index = NULL) {
  unlist(lapply(inits, function(x) if (is.null(index)) x[[name]] else x[[name]][index]))
}

test_that("a random slope on trial number starts with an SD scaled to the predictor and the range (#621)", {
  withr::local_seed(621)
  inits <- bounded_ndt_inits(bmf(drift ~ 1, bound ~ 1, ndt ~ trial + (1 + trial | id)))
  slope_sd <- draws_of(inits, "sd_1", 2)
  expect_lte(max(slope_sd), 0.1 * 0.4055 / 40)
  expect_gte(min(slope_sd), 0.05 * 0.4054 / 40)
  intercept_sd <- draws_of(inits, "sd_1", 1)
  expect_lte(max(intercept_sd), 0.1 * 0.4055)
  expect_gte(min(intercept_sd), 0.05 * 0.4054)
})

test_that("a random intercept on an identity link starts with an SD scaled to the range (#621)", {
  withr::local_seed(621)
  inits <- bounded_ndt_inits(bmf(drift ~ 1, bound ~ 1, ndt ~ 1 + (1 | id)), links = list(ndt = "identity"))
  sd <- draws_of(inits, "sd_1")
  expect_lte(max(sd), 0.1 * 0.006 + 1e-12)
  expect_gte(min(sd), 0.05 * 0.006 - 1e-12)
  expect_true(all(abs(draws_of(inits, "z_1")) <= 0.5))
})

test_that("the SDs of one group term across parameters each follow their own parameter and column (#621)", {
  withr::local_seed(621)
  # drift (log link, c(1, 2)) is log(2) = 0.6931 wide, bound (log, c(1.5, 2))
  # log(4 / 3) = 0.2877; the columns are 1 except trial, which reaches 40
  inits <- bounded_ndt_inits(bmf(
    drift ~ 1 + (1 + trial | p | id), bound ~ 1 + (1 | id), ndt ~ trial + (1 + trial | p | id)
  ))
  upper <- c(0.6931, 0.6931 / 40, 0.4055, 0.4055 / 40)
  for (k in 1:4) {
    sd <- draws_of(inits, "sd_1", k)
    expect_lte(max(sd), 0.1 * upper[k] + 1e-6, label = paste("sd_1 entry", k))
    expect_gte(min(sd), 0.05 * upper[k] - 1e-6, label = paste("sd_1 entry", k))
  }
  expect_lte(max(draws_of(inits, "sd_2")), 0.1 * 0.2877 + 1e-6)
  expect_gte(min(draws_of(inits, "sd_2")), 0.05 * 0.2877 - 1e-6)
})

test_that("an mm() term reads its membership suffix as a member, not as an entry (#621)", {
  withr::local_seed(621)
  # brms names its columns Z_1_drift_<entry>_<member>
  inits <- bounded_ndt_inits(bmf(drift ~ 1 + (1 + trial | mm(id, id2)), bound ~ 1, ndt ~ 1))
  intercept_sd <- draws_of(inits, "sd_1", 1)
  expect_gte(min(intercept_sd), 0.05 * 0.6931 - 1e-6)
  expect_lte(max(intercept_sd), 0.1 * 0.6931 + 1e-6)
  expect_lte(max(draws_of(inits, "sd_1", 2)), 0.1 * 0.6931 / 40 + 1e-6)
})

test_that("a gr(by = ) term shrinks each row of its SD matrix by that row's column (#621)", {
  withr::local_seed(621)
  inits <- bounded_ndt_inits(bmf(drift ~ 1, bound ~ 1, ndt ~ 1 + (1 + trial | gr(id, by = g))))
  expect_equal(dim(inits[[1]]$sd_1), c(2, 2))
  slope_sd <- unlist(lapply(inits, function(x) x$sd_1[2, ]))
  expect_lte(max(slope_sd), 0.1 * 0.4055 / 40)
  intercept_sd <- unlist(lapply(inits, function(x) x$sd_1[1, ]))
  expect_gte(min(intercept_sd), 0.05 * 0.4054)
})

test_that("an s() term starts its spline near zero (#621)", {
  withr::local_seed(621)
  inits <- bounded_ndt_inits(bmf(drift ~ 1, bound ~ 1, ndt ~ 1 + s(trial, k = 5)))
  # the s(trial) columns brms builds stay within -1 and 1, so only the range shrinks
  sds <- draws_of(inits, "sds_ndt_1")
  expect_lte(max(sds), 0.1 * 0.4055)
  expect_gte(min(sds), 0.05 * 0.4054)
  expect_true(all(abs(draws_of(inits, "zs_ndt_1_1")) <= 0.5))
  bs <- draws_of(inits, "bs_ndt")
  expect_lte(max(abs(bs)), 0.1 * 0.4055)
  expect_gt(max(abs(bs)), 0)
})

test_that("a mo() coefficient starts near zero, divided by the categories it spans (#621)", {
  withr::local_seed(621)
  # setsize has three levels, so mo() multiplies its coefficient by up to 2
  inits <- bounded_ndt_inits(bmf(drift ~ 1, bound ~ 1, ndt ~ 1 + mo(setsize)))
  bsp <- draws_of(inits, "bsp_ndt")
  expect_lte(max(abs(bsp)), 0.1 * 0.4055 / 2)
  expect_gt(max(abs(bsp)), 0.05 * 0.4054 / 2)
  simo <- draws_of(inits, "simo_ndt_1")
  expect_true(all(simo > 0 & simo < 1))
})

test_that("unit columns on a range wider than 1 keep the SD draws they had before #621", {
  dat <- small_data()
  # written from the base of #621 with the same seed: the shrink factor is 1
  # for every entry here, so the draws must not move
  cases <- list(
    mixture2p = list(
      model = mixture2p("dev_rad"),
      formula = bmf(kappa ~ 0 + set_size + (0 + set_size | ID), thetat ~ 0 + set_size + (0 + set_size | ID)),
      data = dat,
      sd = list(
        sd_1 = c(0.05957906, 0.06431895, 0.09658438, 0.0689548, 0.06884546, 0.0503031, 0.05774589, 0.09981795),
        sd_2 = c(0.0915609, 0.09422463, 0.09371997, 0.05057796, 0.05859956, 0.08822028, 0.09902488, 0.08486954)
      )
    ),
    m3 = list(
      model = m3(
        resp_cats = c("corr", "other", "npl"), num_options = c("n_corr", "n_other", "n_npl"),
        choice_rule = "softmax", version = "ss"
      ),
      formula = bmf(c ~ 1 + cond + (1 + cond | ID), a ~ 1 + cond + (1 + cond | ID)),
      data = oberauer_lewandowsky_2019_e1,
      sd = list(sd_1 = c(0.06884546, 0.0503031, 0.05774589), sd_2 = c(0.06797572, 0.05016843, 0.09157917))
    ),
    imm = list(
      model = imm(
        resp_error = "dev_rad", nt_features = paste0("col_nt", 1:7),
        nt_distances = paste0("dist_nt", 1:7), set_size = "set_size"
      ),
      formula = bmf(c ~ 0 + set_size + (0 + set_size | ID), a ~ 1 + (1 | ID), s ~ 1, kappa ~ 1),
      data = dat,
      sd = list(
        sd_1 = 0.09658438,
        sd_2 = c(0.0503031, 0.05774589, 0.09981795, 0.06557165, 0.07012099, 0.05040553, 0.08198434, 0.080004)
      )
    ),
    mixture3p = list(
      model = mixture3p("dev_rad", nt_features = paste0("col_nt", 1:7), set_size = "set_size"),
      formula = bmf(
        kappa ~ 0 + set_size + (0 + set_size | ID), thetat ~ 0 + set_size + (0 + set_size | ID),
        thetant ~ 0 + set_size + (0 + set_size | ID)
      ),
      data = dat,
      sd = list(
        sd_1 = c(0.06431895, 0.09658438, 0.0689548, 0.06884546, 0.0503031, 0.05774589, 0.09981795, 0.06557165),
        sd_2 = c(0.09422463, 0.09371997, 0.05057796, 0.05859956, 0.08822028, 0.09902488, 0.08486954, 0.0583179),
        sd_3 = c(0.07510217, 0.07860099, 0.08637935, 0.06904934, 0.05220266, 0.07224673, 0.06396304, 0.06847676)
      )
    )
  )
  for (name in names(cases)) {
    case <- cases[[name]]
    # set size 1 has no non-target, which mixture3p warns about for thetant
    init_fun <- suppressWarnings(configured_initfun(case$model, case$formula, case$data))
    inits <- withr::with_seed(621, init_fun())
    for (sd in names(case$sd)) {
      expect_equal(as.vector(inits[[sd]]), case$sd[[sd]], tolerance = 1e-6, label = paste(name, sd))
    }
  }
})

# =============================================================================
# BOUNDS THE USER SETS ON COEFFICIENTS (#694)
# =============================================================================

# Stan rejects a start on or past a bound, and puts one a rounding error inside
# it at the far end of the unconstrained scale. A start keeps 0.01 times its
# coefficient's shrink from the bound, or a twentieth of the interval between
# two bounds where that is less. The shrink is 1 for a column reaching at most 1
# under an init range at least 1 wide on the link scale, as in eta_model. A start
# at -Inf or Inf clears a margin from an infinite bound but is no start at all
inside <- function(x, lower = -Inf, upper = Inf, shrink = 1) {
  margin <- min(0.01 * shrink, (upper - lower) / 20)
  all(is.finite(x) & x >= lower + margin & x <= upper - margin)
}

eta_model <- list(
  parameters = list(eta = ""), init_ranges = list(eta = c(2, 3)), links = list(eta = "identity")
)
eta_data <- list(X_eta = stats::model.matrix(~ 1 + x, data.frame(x = factor(1:3))))

test_that("solved coefficients start inside a bound set on the class (#694)", {
  withr::local_seed(694)
  draw <- function(bounds) init_fixef_param("b_eta", "vector", 3, eta_model, eta_data, bounds)
  lower <- replicate(200, draw(c(0, Inf)))
  upper <- replicate(200, draw(c(-Inf, 0)))
  narrow <- replicate(200, draw(c(-1e-4, 1e-4)))
  # the effects the constant target does not need solve to zero
  expect_true(inside(lower[2:3, ], lower = 0))
  expect_true(all(lower[2:3, ] <= 0.1))
  expect_true(inside(upper, upper = 0))
  expect_true(all(narrow > -1e-4 & narrow < 1e-4))
  # an intercept above the bound keeps its start in range
  expect_true(all(lower[1, ] >= 2 & lower[1, ] <= 3))
})

test_that("a bounded coefficient on a wide column moves by a step its column shrinks (#694)", {
  withr::local_seed(694)
  # a centred design (Kc_eta): coefficient k belongs to column k + 1 of X
  sdata <- list(X_eta = cbind(Intercept = 1, small = rep(0:1, 20), wide = 1:40), Kc_eta = 2L)
  b <- replicate(200, init_fixef_param("b_eta", "vector", 2, eta_model, sdata, c(0, Inf)))
  expect_true(inside(b[1, ], lower = 0) && all(b[1, ] <= 0.1))
  expect_true(inside(b[2, ], lower = 0, shrink = 1 / 40) && all(b[2, ] <= 0.1 / 40))
})

test_that("a start that moves takes its own shrink when another start stays (#694)", {
  withr::local_seed(694)
  starts <- replicate(200, into_bounds(c(0.5, -1e-17), c(0, Inf), c(1, 1 / 40)))
  expect_true(all(starts[1, ] == 0.5))
  expect_true(inside(starts[2, ], lower = 0, shrink = 1 / 40) && all(starts[2, ] <= 0.1 / 40))
})

test_that("starts between close bounds spread strictly inside them (#694)", {
  withr::local_seed(694)
  close <- replicate(200, into_bounds(c(-1, 1), c(0, 0.02), 1))
  expect_true(inside(close, lower = 0, upper = 0.02))
  expect_gt(length(unique(close[1, ])), 1)
  expect_gt(length(unique(close[2, ])), 1)
  # four doubles apart, some steps fall below half the spacing of the doubles
  # and round onto the bound
  tiny <- c(1, 1 + 4 * .Machine$double.eps)
  starts <- replicate(200, into_bounds(c(0, 2), tiny, 1))
  expect_true(all(starts > tiny[1] & starts < tiny[2]))
})

test_that("a start well inside the bounds is the start without them (#694)", {
  # the intercept starts in 2-3 and the effects solve to zero
  draw <- function(bounds) withr::with_seed(1, init_fixef_param("b_eta", "vector", 3, eta_model, eta_data, bounds))
  expect_equal(draw(c(-1, 4)), draw(c(-Inf, Inf)))
})

test_that("a bound that never binds leaves every seeded start as it was (#694)", {
  inits <- function(user_prior) {
    init_fun <- configured_initfun(
      mixture2p(resp_error = "response_error"),
      bmf(kappa ~ 1 + setsize + (1 | subID), thetat ~ 1 + setsize),
      zhang_luck_2008,
      user_prior = user_prior
    )
    withr::with_seed(5, init_fun())
  }
  expect_identical(
    inits(brms::prior(normal(0, 1), class = b, dpar = kappa, lb = -5)),
    inits(brms::prior(normal(0, 1), class = b, dpar = kappa))
  )
})

test_that("an intercept and a per-coefficient scalar start inside their bounds (#694)", {
  withr::local_seed(694)
  intercept <- replicate(200, init_fixef_param("Intercept_eta", "real", 1, eta_model, list(), c(5, Inf)))
  effect <- replicate(200, init_fixef_param("par_b_eta_3", "real", 1, eta_model, eta_data, c(0, Inf)))
  expect_true(inside(intercept, lower = 5))
  # replicate() would fold a vector returned for the scalar into a matrix
  expect_length(effect, 200)
  expect_true(inside(effect, lower = 0))
  below <- replicate(200, init_fixef_param("Intercept_eta", "real", 1, eta_model, list(), c(-Inf, 0)))
  expect_true(inside(below, upper = 0) && all(below >= -0.1))
})

test_that("an intercept on a narrow init range moves by a step that range shrinks (#694)", {
  withr::local_seed(694)
  # the shrink is the 0.006 width of the range, which no column widens
  model <- list(parameters = list(ndt = ""), init_ranges = list(ndt = c(0.012, 0.018)), links = list(ndt = "identity"))
  intercept <- replicate(200, init_fixef_param("Intercept_ndt", "real", 1, model, list(), c(0.05, Inf)))
  expect_true(inside(intercept, lower = 0.05, shrink = 0.006) && all(intercept <= 0.05 + 0.1 * 0.006))
})

test_that("mo() coefficients start inside a bound set on the class (#694)", {
  withr::local_seed(694)
  model <- list(parameters = list(kappa = ""), init_ranges = list(kappa = c(1, 2)), links = list(kappa = "log"))
  sdata <- list(Jmo_kappa = 3L)
  bsp <- replicate(200, init_special_coef_param("bsp_kappa", "vector", 2, model, sdata, c(0, Inf)))
  # the range is log(2) wide on the log scale, and mo() reaches Jmo = 3
  shrink <- log(2) / 3
  expect_true(inside(bsp, lower = 0, shrink = shrink) && all(bsp <= 0.1 * shrink))
  below <- replicate(200, init_special_coef_param("bsp_kappa", "vector", 2, model, sdata, c(-Inf, 0)))
  expect_true(inside(below, upper = 0, shrink = shrink) && all(below >= -0.1 * shrink))
})

test_that("the inits of a fit with lb = 0 on b respect the bound (#694)", {
  withr::local_seed(694)
  init_fun <- configured_initfun(
    mixture2p(resp_error = "response_error"),
    bmf(kappa ~ 1 + setsize, thetat ~ 1 + setsize),
    zhang_luck_2008,
    user_prior = brms::prior(normal(0, 1), class = b, dpar = kappa, lb = 0)
  )
  inits <- lapply(1:50, function(i) init_fun())
  env <- environment(init_fun)
  # the coefficients belong to every column of the centred design but the first
  shrink <- apply(env$standata_list$X_kappa[, -1, drop = FALSE], 2, function(column) {
    coef_shrink("kappa", column, env$model)
  })
  expect_true(inside(unlist(lapply(inits, `[[`, "b_kappa")), lower = 0, shrink = min(shrink)))
  expect_true(all(is.finite(unlist(inits))))
})

# SDs start at 0.05 to 0.1 times their shrink, so a bound on either side of
# that interval binds every draw
sd_starts <- function(formula, user_prior, sd_name, data = zhang_luck_2008) {
  init_fun <- configured_initfun(mixture2p(resp_error = "response_error"), formula, data, user_prior)
  env <- environment(init_fun)
  list(
    starts = unlist(lapply(1:20, function(i) init_fun()[[sd_name]])),
    shrink = min(sd_entry_shrink(sd_name, env$model, env$standata_list))
  )
}

test_that("group-level SDs start apart inside an upper bound set on the class (#694)", {
  withr::local_seed(694)
  for (ub in c(0.02, 0.04)) {
    sd <- sd_starts(
      bmf(kappa ~ 1 + (1 | subID), thetat ~ 1),
      brms::set_prior("normal(0, 1)", class = "sd", dpar = "kappa", ub = ub),
      "sd_1"
    )
    expect_true(inside(sd$starts, lower = 0, upper = ub, shrink = sd$shrink))
    expect_gt(length(unique(sd$starts)), 1)
  }
})

test_that("group-level SDs start inside a lower bound set on the class (#694)", {
  withr::local_seed(694)
  sd <- sd_starts(
    bmf(kappa ~ 1 + (1 | subID), thetat ~ 1),
    brms::prior(normal(0, 1), class = sd, dpar = kappa, lb = 0.5),
    "sd_1"
  )
  expect_true(inside(sd$starts, lower = 0.5, shrink = sd$shrink))
  expect_true(all(sd$starts <= 0.5 + 0.1 * sd$shrink))
})

test_that("the SDs of an s() term start inside a bound set on the class (#694)", {
  withr::local_seed(694)
  dat <- zhang_luck_2008
  dat$z <- dat$trial / 100
  sds <- sd_starts(
    bmf(kappa ~ 1 + s(z, k = 4), thetat ~ 1),
    brms::prior(normal(0, 1), class = sds, dpar = kappa, ub = 0.02),
    "sds_kappa_1",
    dat
  )
  expect_true(inside(sds$starts, lower = 0, upper = 0.02, shrink = sds$shrink))
  expect_gt(length(unique(sds$starts)), 1)
})
