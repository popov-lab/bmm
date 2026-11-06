save_pars <- brms::save_pars

test_that("update.bmmfit works", {
  skip_on_cran()
  withr::local_seed(445)
  fit1 <- restructure(load_fixture_fit("bmmfit_example1.rds"))
  data <- fit1$data

  # formula is replaced
  up <- suppressMessages(update(fit1, formula. = bmf(c ~ 1, kappa ~ 1), testmode = TRUE))
  expect_true(is(up, "bmmfit"))
  expect_equal(up$bmm$user_formula$c, c ~ 1, ignore_attr = TRUE)

  # data is replaced, old formula is kept
  new_data <- data
  new_data$dev_rad <- rnorm(nrow(new_data), 0, 0.5)
  up <- suppressMessages(
    update(fit1, newdata = new_data, save_pars = save_pars(group = FALSE), testmode = TRUE)
  )
  expect_true(is(up, "bmmfit"))
  expect_equal(attr(up$data, "data_name"), "new_data")
  expect_equal(up$bmm$user_formula$c, c ~ 0 + set_size, ignore_formula_env = T, ignore_attr = TRUE)

  # prior is replaced
  up <- suppressMessages(
    update(
      fit1,
      formula. = bmf(c ~ 1, kappa ~ 1), testmode = TRUE,
      prior = brms::set_prior("normal(0,0.1)", class = "Intercept", dpar = "kappa")
    )
  )
  expect_true(is(up, "bmmfit"))

  # refuse to change model
  expect_error(
    update(fit1, model = mixture2p(resp_error = "dev_rad")),
    "You cannot update with a different model"
  )

  up <- suppressMessages(update(fit1, save_pars = save_pars(group = FALSE), testmode = TRUE))
  expect_true(is(up, "bmmfit"))
  up <- suppressMessages(update(fit1, save_pars = save_pars(latent = FALSE), testmode = TRUE))
  expect_true(is(up, "bmmfit"))
  expect_error(update(fit1, data = data), "use argument 'newdata'")
})

sdm_fixture <- function() {
  restructure(load_fixture_fit("bmmfit_example1.rds"))
}

update_mock <- function(object, ...) {
  suppressMessages(update(object, ..., backend = "mock", mock_fit = 1, rename = FALSE))
}

constant_priors <- function(fit) {
  sum(grepl("^constant\\(", as.data.frame(fit$prior)$prior))
}

test_that("update() stores the formula that bmm() stores for the same inputs", {
  # brms::combine_models() compares the formula elements by position, so an
  # updated fit merges with the original only if the element order matches too
  skip_on_cran()
  fit1 <- sdm_fixture()
  bmm_mock <- function(formula) {
    suppressMessages(bmm(formula, fit1$data, fit1$bmm$model,
                         backend = "mock", mock_fit = 1, rename = FALSE))
  }

  expect_equal(formula(update_mock(fit1)),
               formula(bmm_mock(fit1$bmm$user_formula)),
               ignore_formula_env = TRUE)

  new_formula <- bmf(c ~ 1, kappa ~ 1)
  expect_equal(formula(update_mock(fit1, formula. = new_formula)),
               formula(bmm_mock(new_formula)),
               ignore_formula_env = TRUE)
})

test_that("update() stores the formula that bmm() stores, without a fixture", {
  skip_on_cran() # over 0.4 s in CRAN mode (ADR 0004, R10)
  skip_if_not_installed("rstan")
  withr::local_seed(445)
  stub <- methods::new("stanfit", sim = list(
    iter = 10L, warmup = 5L, chains = 1L, thin = 1L,
    samples = list(structure(list(), args = list(control = list())))
  ))
  dat <- data.frame(y = rsdm(60, kappa = 5))
  bmm_mock <- function(formula) {
    suppressMessages(bmm(formula, dat, sdm("y"), backend = "mock", mock_fit = stub, rename = FALSE))
  }
  fit <- bmm_mock(bmf(c ~ 1, kappa ~ 1))
  # a mock fit cannot be reused, so update() takes the recompile path
  update_mock <- function(...) {
    suppressMessages(update(fit, ..., backend = "mock", mock_fit = stub, rename = FALSE, recompile = TRUE))
  }

  expect_equal(formula(update_mock()), formula(fit), ignore_formula_env = TRUE)

  new_formula <- bmf(c ~ 1, kappa ~ 1, mu ~ 1)
  expect_equal(formula(update_mock(formula. = new_formula)),
               formula(bmm_mock(new_formula)),
               ignore_formula_env = TRUE)
})

test_that("update() that keeps the Stan code names the bmm version once (#551)", {
  skip_on_cran()
  fit1 <- sdm_fixture()
  # the mock backend recompiles and so writes a fresh header; testmode keeps
  # the fit's code, which already carries the bmm version
  expected <- paste0("and bmm ", utils::packageVersion("bmm"))
  header <- function(fit) sub("^.*? (and bmm[^\n]*)\n.*$", "\\1", unclass(fit$model))
  up <- suppressMessages(update(fit1, testmode = TRUE))
  expect_identical(header(up), expected)
  expect_identical(header(suppressMessages(update(up, testmode = TRUE))), expected)
  up <- suppressMessages(update(fit1, formula. = bmf(c ~ 1, kappa ~ 1), testmode = TRUE))
  expect_identical(header(up), expected)
})

test_that("update.bmmfit frees parameters that the new formula predicts", {
  skip_on_cran()
  withr::local_seed(445)
  fit1 <- sdm_fixture()

  # mu is fixed to 0 by default in sdm; predicting it must drop the constant()
  # prior the original fit stored for it
  up <- update_mock(fit1, formula. = bmf(c ~ 1, kappa ~ 1, mu ~ 1))
  expect_length(up$bmm$model$fixed_parameters, 0)
  expect_equal(constant_priors(up), 0)
  expect_match(brms::stancode(up), "real Intercept;", fixed = TRUE)

  # freeing a parameter and replacing the data at once
  new_data <- data.frame(
    dev_rad = rsdm(80, c = 4, kappa = 3),
    set_size = rep(1:2, each = 40)
  )
  up <- update_mock(fit1, formula. = bmf(c ~ 1, kappa ~ 1, mu ~ 1), newdata = new_data)
  expect_length(up$bmm$model$fixed_parameters, 0)
  expect_equal(constant_priors(up), 0)
  expect_equal(nrow(up$data), 80)

  # a parameter left alone stays fixed
  expect_equal(constant_priors(update_mock(fit1)), 1)
})

test_that("update.bmmfit fixes parameters that the new formula sets to a constant", {
  skip_on_cran()
  fit1 <- sdm_fixture()

  # the reverse direction: the old fit's free prior on kappa must not survive and
  # override the constant() that the new formula asks for
  up <- update_mock(fit1, formula. = bmf(c ~ 0 + set_size, kappa = 5))
  expect_equal(up$bmm$model$fixed_parameters$kappa, 5)
  expect_true(any(grepl("^constant\\(5\\)", as.data.frame(up$prior)$prior)))
  expect_match(brms::stancode(up), "Intercept_kappa = 5;", fixed = TRUE)
  expect_false(grepl("student_t_lpdf(Intercept_kappa", brms::stancode(up), fixed = TRUE))

  # changing the value of an already fixed parameter is the same defect: the old
  # constant() row would otherwise survive and pin mu at the original 0
  up <- update_mock(fit1, formula. = bmf(c ~ 0 + set_size, mu = 0.5))
  expect_equal(up$bmm$model$fixed_parameters$mu, 0.5)
  expect_match(brms::stancode(up), "Intercept = 0.5;", fixed = TRUE)
})

test_that("a predictor added by update() gets the prior bmm() gives it", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  withr::local_seed(445)
  stub <- methods::new("stanfit", sim = list(
    iter = 10L, warmup = 5L, chains = 1L, thin = 1L,
    samples = list(structure(list(), args = list(control = list())))
  ))
  stub_fit <- function(formula, model) suppressWarnings(suppressMessages(
    bmm(formula, dat, model, backend = "mock", mock_fit = stub, rename = FALSE)
  ))
  dat <- data.frame(y = rmixture2p(n = 60), cond = rep(c("a", "b"), 30))
  cases <- list(
    nlpar = list(model = mixture2p("y"), old = bmf(kappa ~ 1, thetat ~ 1),
                 new = bmf(kappa ~ 1 + cond, thetat ~ 1)),
    dpar = list(model = sdm("y"), old = bmf(c ~ 1, kappa ~ 1),
                new = bmf(c ~ 1 + cond, kappa ~ 1))
  )

  for (case_name in names(cases)) {
    case <- cases[[case_name]]
    # the fit stores brms's flat default row for the class of the old formula's
    # parameter; carried over as a user prior, it overrode bmm's prior on the slope
    up <- suppressWarnings(suppressMessages(update(
      stub_fit(case$old, case$model), formula. = case$new, newdata = dat,
      backend = "mock", mock_fit = stub, rename = FALSE, recompile = TRUE
    )))
    expect_equal(brms::stancode(up), brms::stancode(stub_fit(case$new, case$model)),
                 label = case_name)
  }
})

# read the emitted program rather than bmm's own stanvars, and assert on the
# call: the sliced function is declared in the functions block either way
sdm_likelihood_is_sliced <- function(fit) {
  grepl("target += sdm_simple_run_ldenom_slice", brms::stancode(fit), fixed = TRUE)
}

test_that("update.bmmfit configures the likelihood for the effective threading spec", {
  skip_on_cran()
  fit1 <- sdm_fixture()

  # brms::update.brmsfit falls back to the fit's own threading spec when threads
  # is not passed, so the threaded chunk must be emitted without it too
  threaded_fit <- fit1
  threaded_fit$threads <- brms::threading(2)
  expect_true(sdm_likelihood_is_sliced(update_mock(threaded_fit)))

  # brms reads an explicit threads = NULL as "turn threading off", which it
  # distinguishes from an absent argument by name presence, not by NULL-ness
  expect_false(sdm_likelihood_is_sliced(update_mock(threaded_fit, threads = NULL)))

  expect_false(sdm_likelihood_is_sliced(update_mock(fit1)))
  expect_true(sdm_likelihood_is_sliced(update_mock(fit1, threads = brms::threading(2))))

  # the spec of the fit also has to win over a global brms.threads, which brms
  # itself ignores once it has fallen back to the fit's own spec
  withr::with_options(list(brms.threads = brms::threading(2)), {
    expect_false(sdm_likelihood_is_sliced(update_mock(fit1)))

    # a fit saved before brmsfit carried a $threads field
    no_threads <- fit1
    no_threads$threads <- NULL
    expect_false(sdm_likelihood_is_sliced(update_mock(no_threads)))
  })
})

test_that("update.bmmfit keeps track of the file the fit is saved in", {
  skip_on_cran()
  fit1 <- sdm_fixture()
  fit1$file <- "some/cached/fit.rds"

  expect_equal(update_mock(fit1)$file, "some/cached/fit.rds")

  # an explicit file must hold a bmmfit -- brms writes it from inside
  # update.brmsfit(), before any of the bmm postprocessing has run
  file <- tempfile()
  up <- update_mock(fit1, file = file)
  expect_equal(up$file, paste0(file, ".rds"))
  expect_true(file.exists(paste0(file, ".rds")))
  expect_true(is_bmmfit(readRDS(paste0(file, ".rds"))))
})

test_that("update.bmmfit updates and writes even when `file` already exists", {
  skip_on_cran()
  fit1 <- sdm_fixture()

  # `file` must not reach brms::brm(), which would read an existing file and
  # return its contents instead of running the update
  file <- tempfile()
  saveRDS(fit1, paste0(file, ".rds"))
  new_data <- fit1$data
  new_data$dev_rad <- 0.1
  up <- update_mock(fit1, newdata = new_data, file = file)
  expect_equal(up$fit, 1)
  expect_equal(readRDS(paste0(file, ".rds"))$fit, 1)
  expect_true(is_bmmfit(readRDS(paste0(file, ".rds"))))
})

test_that("update() builds the initial values for the data and formula it fits", {
  skip_on_cran()
  fit1 <- sdm_fixture()
  # on the same backend brms carries every stored stan_args entry the call does
  # not name, the init among them, and it leaves stan_args untouched unless the
  # model is recompiled
  fit1$backend <- "mock"
  new_data <- fit1$data[fit1$data$set_size %in% 1:2, ]
  new_data$set_size <- factor(as.character(new_data$set_size))
  new_data$ID <- factor(rep(1:3, length.out = nrow(new_data)))

  up <- update_mock(fit1, newdata = new_data, recompile = TRUE)
  expect_length(up$stan_args$init()$b_c, brms::standata(up)$K_c)

  up <- update_mock(fit1,
    formula. = bmf(c ~ 0 + set_size + (1 | ID), kappa ~ 1), newdata = new_data, recompile = TRUE
  )
  expect_equal(ncol(up$stan_args$init()$z_1), brms::standata(up)$N_1)

  expect_equal(update_mock(fit1, newdata = new_data, init = 0, recompile = TRUE)$stan_args$init, 0)
})

test_that("update() without recompiling stores the init it sampled with (#520)", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  stub <- methods::new("stanfit", sim = list(
    iter = 10L, warmup = 5L, chains = 1L, thin = 1L,
    samples = list(structure(list(), args = list(control = list())))
  ))
  withr::local_seed(520)
  dat <- data.frame(y = rsdm(120, c = 4, kappa = 5), set_size = factor(rep(1:3, each = 40)))
  fit <- suppressMessages(bmm(bmf(c ~ 0 + set_size, kappa ~ 1), dat, sdm("y"),
                              backend = "mock", mock_fit = stub, rename = FALSE))
  new_data <- droplevels(dat[dat$set_size %in% 1:2, ])
  # rstan and cmdstanr fits reuse their compiled model, and brms then keeps the
  # old fit's stan_args; a mock fit has no model to reuse until compiled_model()
  # returns one, and a recompile, which rebuilds stan_args, fails the test
  local_mocked_bindings(
    compiled_model = function(x) NULL,
    compile_model = function(...) stop("update() recompiled the model"),
    .package = "brms"
  )
  update_reuse <- function(...) suppressWarnings(suppressMessages(
    update(fit, newdata = new_data, ..., mock_fit = stub, rename = FALSE)
  ))

  expect_length(fit$stan_args$init()$b_c, 3)
  expect_length(update_reuse()$stan_args$init()$b_c, 2)
  expect_identical(update_reuse(init = 0)$stan_args$init, 0)
  expect_identical(update_reuse(inits = 0)$stan_args$init, 0)
})

test_that("update() configures the prior and the inits with the fit's data2", {
  skip_on_cran()
  fit1 <- sdm_fixture()
  new_data <- fit1$data
  new_data$ID <- factor(rep(1:3, length.out = nrow(new_data)))
  A <- diag(3)
  dimnames(A) <- list(levels(new_data$ID), levels(new_data$ID))
  formula <- bmf(c ~ 0 + set_size + (1 | gr(ID, cov = A)), kappa ~ 1)

  # the call passes data2
  up <- update_mock(fit1, formula. = formula, newdata = new_data, data2 = list(A = A))
  expect_equal(up$data2, list(A = A))
  expect_equal(dim(up$stan_args$init()$z_1), c(1, 3))
  # the fit carries it, as brms::update.brmsfit() reads it off object$data2
  fit1$data2 <- list(A = A)
  up <- update_mock(fit1, formula. = formula, newdata = new_data)
  expect_equal(dim(up$stan_args$init()$z_1), c(1, 3))
})

test_that("update() applies the package step-size default to a fit that had none", {
  skip_on_cran()
  withr::local_options(bmm.step_size = 0.02)
  fit1 <- sdm_fixture()
  # as a fit from before the package had a step-size default
  fit1$stan_args$control <- NULL

  up <- update_mock(fit1)
  expect_equal(up$stan_args$control, list(step_size = 0.02))
  up <- update_mock(fit1, control = list(adapt_delta = 0.99))
  expect_equal(up$stan_args$control, list(adapt_delta = 0.99, step_size = 0.02))
  up <- update_mock(fit1, control = list(step_size = 0.5))
  expect_equal(up$stan_args$control, list(step_size = 0.5))

  withr::local_options(bmm.step_size = FALSE)
  expect_null(update_mock(fit1)$stan_args$control)
})

test_that("update() keeps the fit's control on the same backend and algorithm only", {
  object <- list(backend = "cmdstanr", algorithm = "sampling", stan_args = list(
    control = list(adapt_delta = 0.95)
  ))
  expect_equal(carried_control(object, list()), list(adapt_delta = 0.95))
  expect_equal(carried_control(object, list(backend = "cmdstanr")), list(adapt_delta = 0.95))
  # brms's rule: the call's keys win, the fit's remaining keys are kept
  expect_equal(
    carried_control(object, list(control = list(max_treedepth = 12))),
    list(adapt_delta = 0.95, max_treedepth = 12)
  )
  expect_equal(
    carried_control(object, list(control = list(adapt_delta = 0.8))),
    list(adapt_delta = 0.8)
  )
  expect_null(carried_control(object, list(backend = "rstan")))
  expect_null(carried_control(object, list(algorithm = "meanfield")))
  # brms resolves a missing field to its first choice and then finds it changed
  no_backend <- object
  no_backend$backend <- NULL
  expect_null(carried_control(no_backend, list()))
  no_algorithm <- object
  no_algorithm$algorithm <- NULL
  expect_null(carried_control(no_algorithm, list()))
})

# brms reuses the compiled model, and with it the old fit's stan_args, only on a
# real backend, so the mocks stand in for one: compile_model() fails the test if
# update() recompiles, and fit_model() records the control the sampler gets
local_reused_model_fit <- function(..., env = parent.frame()) {
  stub <- methods::new("stanfit", sim = list(
    iter = 10L, warmup = 5L, chains = 1L, thin = 1L,
    samples = list(structure(list(), args = list(control = list())))
  ))
  dat <- withr::with_seed(520, data.frame(y = rsdm(60, c = 4, kappa = 5)))
  fit <- suppressMessages(bmm(bmf(c ~ 1, kappa ~ 1), dat, sdm("y"),
    backend = "mock", mock_fit = stub, rename = FALSE, ...))
  sampled <- new.env()
  local_mocked_bindings(
    compiled_model = function(x) NULL,
    compile_model = function(...) stop("update() recompiled the model"),
    fit_model = function(...) {
      sampled$control <- list(...)$control
      stub
    },
    .package = "brms", .env = env
  )
  list(fit = fit, sampled = sampled)
}

update_reused <- function(object, ...) {
  suppressMessages(update(object, ..., rename = FALSE))
}

test_that("a chained update() without recompiling keeps the previous update's adapt_delta", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  run <- local_reused_model_fit(control = list(adapt_delta = 0.9))
  up1 <- update_reused(run$fit, control = list(adapt_delta = 0.95))
  expect_equal(run$sampled$control$adapt_delta, 0.95)
  expect_equal(up1$stan_args$control$adapt_delta, 0.95)

  up2 <- update_reused(up1)
  expect_equal(run$sampled$control$adapt_delta, 0.95)
  expect_equal(up2$stan_args$control, list(adapt_delta = 0.95, step_size = 0.01))
})

test_that("a chained update() without recompiling keeps the previous update's step size", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  run <- local_reused_model_fit()
  up1 <- update_reused(run$fit, control = list(step_size = 0.05))
  expect_equal(up1$stan_args$control, list(step_size = 0.05))

  up2 <- update_reused(up1)
  expect_equal(run$sampled$control, list(step_size = 0.05))
  expect_equal(up2$stan_args$control, list(step_size = 0.05))
})

test_that("a chained update() after a switch of algorithm hands the sampler no stale control", {
  # cmdstanr passes every control key to the method it runs, and $variational()
  # rejects the sampler's step size
  skip_on_cran()
  skip_if_not_installed("rstan")
  run <- local_reused_model_fit(control = list(adapt_delta = 0.9))
  up1 <- update_reused(run$fit, algorithm = "meanfield")
  expect_true("control" %in% names(up1$stan_args))
  expect_null(up1$stan_args$control)

  update_reused(up1)
  expect_equal(run$sampled$control, list())
})

test_that("an abbreviated algorithm or backend is the fit's own run, as brms reads it", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  run <- local_reused_model_fit(control = list(step_size = 0.05))
  up1 <- update_reused(run$fit, algorithm = "samp")
  expect_equal(run$sampled$control, list(step_size = 0.05))
  expect_equal(up1$stan_args$control, list(step_size = 0.05))

  update_reused(up1)
  expect_equal(run$sampled$control, list(step_size = 0.05))

  up2 <- update_reused(run$fit, backend = "mo")
  expect_equal(up2$stan_args$control, list(step_size = 0.05))
})

test_that("update() fills the control from the old rstan sampler as brms does", {
  skip_if_not_installed("rstan")
  stub <- methods::new("stanfit", sim = list(samples = list(structure(list(), args = list(
    control = list(adapt_delta = 0.8, max_treedepth = 10, stepsize = 1)
  )))))
  object <- list(backend = "rstan", algorithm = "sampling", fit = stub)
  expect_equal(
    sampler_control(object, list(), list(adapt_delta = 0.95, stepsize = 0.01)),
    list(adapt_delta = 0.95, stepsize = 0.01, max_treedepth = 10)
  )
  expect_equal(
    sampler_control(object, list(algorithm = "meanfield"), list(adapt_delta = 0.95)),
    list(adapt_delta = 0.95)
  )
})

test_that("update(init = NULL) without recompiling keeps the init slot in stan_args", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  run <- local_reused_model_fit()
  up <- update_reused(run$fit, init = NULL)
  expect_true("init" %in% names(up$stan_args))
  expect_null(up$stan_args$init)
})

test_that("update() carries the stored prior rows that are not brms defaults", {
  stored <- brms::prior(normal(0, 1), class = b) +
    brms::prior(normal(0, 2), class = b, coef = x) +
    brms::prior(normal(0, 3), class = sd) +
    brms::prior(normal(0, 4), class = Intercept)
  # "" is what brms's restructure() writes into fits that predate the column
  stored$source <- c("user", "default", "vectorized", "")
  kept <- carried_prior(stored)
  expect_s3_class(kept, "brmsprior")
  expect_equal(kept$prior, c("normal(0, 1)", "normal(0, 3)", "normal(0, 4)"))
  expect_equal(kept$source, c("user", "vectorized", ""))
})

test_that("update() carries no stored prior row whose source is unknown (#552)", {
  stored <- brms::prior(normal(0, 1), class = b) +
    brms::prior(normal(0, 2), class = sd) +
    brms::prior(normal(0, 3), class = Intercept)
  stored$source <- c("user", NA, "default")
  kept <- carried_prior(stored)
  expect_equal(kept$prior, "normal(0, 1)")
  expect_false(anyNA(kept$class))

  stored$source <- NULL
  kept <- carried_prior(stored)
  expect_s3_class(kept, "brmsprior")
  expect_equal(nrow(kept), 0)
})

test_that("a step size in the call replaces the one the fit stored, either spelling", {
  # configure_control() renames both spellings to the backend's and drops the
  # duplicate, so a merge by name alone would hand the fit's value the win
  object <- list(backend = "rstan", algorithm = "sampling",
                 stan_args = list(control = list(stepsize = 0.01, adapt_delta = 0.95)))
  expect_equal(carried_control(object, list(control = list(step_size = 0.5))),
               list(adapt_delta = 0.95, step_size = 0.5))
  object$stan_args$control <- list(step_size = 0.01)
  expect_equal(carried_control(object, list(control = list(stepsize = 0.5))),
               list(stepsize = 0.5))
})

test_that("update() of a fit with a stored control keeps it next to the step size", {
  skip_on_cran()
  withr::local_options(bmm.step_size = 0.02)
  fit1 <- sdm_fixture()
  fit1$stan_args$control <- list(adapt_delta = 0.95)
  # a new backend starts from brms's defaults, as brms::update.brmsfit() does
  expect_equal(update_mock(fit1)$stan_args$control, list(step_size = 0.02))

  # recompiling keeps brms from reusing the fit's stored stan_args wholesale
  fit1$backend <- "mock"
  control <- update_mock(fit1, recompile = TRUE)$stan_args$control
  expect_equal(control$adapt_delta, 0.95)
  expect_equal(control$step_size, 0.02)

  # a call that names another key leaves the step size the fit was run with
  fit1$stan_args$control <- list(step_size = 0.5)
  control <- update_mock(fit1, recompile = TRUE, control = list(adapt_delta = 0.99))$stan_args$control
  expect_equal(control$adapt_delta, 0.99)
  expect_equal(control$step_size, 0.5)
})

# =============================================================================
# UPDATING WITHOUT newdata (#429)
# =============================================================================

m3_fixture <- function() {
  restructure(load_fixture_fit("bmmfit_m3_ppcheck.rds"))
}

# the stored-frame cases come from the model-case providers
# (helper-model-<model>.R): one mock-fittable (model, formula, data) per model,
# each chosen so that the stored model frame is the hard case, such as an m3
# category with zero options in some rows, or a non-target model whose formula
# leaves brms no reason to keep the set_size column
# built once per file: every provider builds all its cases on each call
stored_frame_cases <- local({
  cases <- NULL
  function() {
    if (is.null(cases)) cases <<- model_test_cases("stored_frame")
    cases
  }
})

test_that("update(formula. = ) refuses only a column the fit did not keep, and says to pass newdata (#552)", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  withr::local_seed(552)
  stub <- methods::new("stanfit", sim = list(
    iter = 10L, warmup = 5L, chains = 1L, thin = 1L,
    samples = list(structure(list(), args = list(control = list())))
  ))
  stub_fit <- function(formula, data, model) suppressWarnings(suppressMessages(
    bmm(formula, data, model, backend = "mock", mock_fit = stub, rename = FALSE)
  ))
  stub_update <- function(fit, ...) suppressWarnings(suppressMessages(
    update(fit, ..., backend = "mock", mock_fit = stub, rename = FALSE, recompile = TRUE)
  ))

  dat <- data.frame(y = rsdm(60), x = rep(c("a", "b"), 30), id = rep(1:6, each = 10))
  fit <- stub_fit(bmf(c ~ 1, kappa ~ 1), dat, sdm("y"))
  expect_error(
    stub_update(fit, formula. = bmf(c ~ 1 + x, kappa ~ 1)),
    "'x'.*`newdata`"
  )
  up <- stub_update(fit, formula. = bmf(c ~ 1 + x, kappa ~ 1), newdata = dat)
  expect_equal(brms::stancode(up),
               brms::stancode(stub_fit(bmf(c ~ 1 + x, kappa ~ 1), dat, sdm("y"))))

  # the correlation ID p and the Intercept of brms are no columns to keep
  f_cor <- bmf(c ~ 1 + x + (1 | p | id), kappa ~ 1 + (1 | p | id))
  cor_fit <- stub_fit(f_cor, dat, sdm("y"))
  expect_no_error(stub_update(cor_fit, formula. = f_cor))
  expect_no_error(stub_update(cor_fit, formula. = bmf(c ~ 0 + Intercept + x + (1 | id), kappa ~ 1)))

  # the tree column is rebuilt for check_data() and dropped again afterwards
  case <- stored_frame_cases()$mpt
  mpt_fit <- stub_fit(case$formula, case$data, case$model)
  expect_error(
    stub_update(mpt_fit, formula. = bmf(D ~ 1 + item_type + (1 | id), g ~ 1)),
    "'item_type'.*`newdata`"
  )
  up <- stub_update(mpt_fit, formula. = bmf(D ~ 1 + item_type + (1 | id), g ~ 1),
                    newdata = case$data)
  expect_true("item_type" %in% colnames(up$data))
})

test_that("a formula whose variables the stored frame or data2 hold is not refused (#552)", {
  configured <- function(formula) list(formula = formula, data = data.frame(y = 1, x = 1, id = "a"))
  # neither the correlation ID p, nor Intercept, nor the cov matrix M, nor the
  # non-linear parameters a and b are data columns
  expect_silent(stop_unstored_variables(configured(brms::bf(
    y ~ 0 + Intercept + x + (1 | p | gr(id, cov = M)), sigma ~ 1 + I(x^2) + (1 | p | id)
  )), data2 = NULL))
  expect_silent(stop_unstored_variables(configured(brms::bf(
    y ~ a + b * x, a ~ 1 + (1 | id), b ~ 1, nl = TRUE
  )), data2 = NULL))
  expect_silent(stop_unstored_variables(configured(brms::bf(y ~ 1 + z)), data2 = list(z = 1)))
  expect_error(
    stop_unstored_variables(configured(brms::bf(y ~ 1 + z + w)), data2 = NULL),
    "'z', 'w'.*`newdata`"
  )
})

for (case_name in names(stored_frame_cases())) {
  test_that(paste0("update() with a fit's own formula and no newdata is not refused (", case_name, ")"), {
    skip_on_cran() # over 0.4 s in CRAN mode (ADR 0004, R10)
    skip_if_not_installed("rstan")
    stub <- methods::new("stanfit", sim = list(
      iter = 10L, warmup = 5L, chains = 1L, thin = 1L,
      samples = list(structure(list(), args = list(control = list())))
    ))
    case <- stored_frame_cases()[[case_name]]
    fit <- suppressWarnings(suppressMessages(
      bmm(case$formula, case$data, case$model, backend = "mock", mock_fit = stub, rename = FALSE)
    ))
    expect_no_error(suppressWarnings(suppressMessages(
      update(fit, formula. = case$formula, backend = "mock", mock_fit = stub,
             rename = FALSE, recompile = TRUE)
    )))
  })
}

test_that("a custom m3 fit updates to a formula with other activation parameters (#552)", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  stub <- methods::new("stanfit", sim = list(
    iter = 10L, warmup = 5L, chains = 1L, thin = 1L,
    samples = list(structure(list(), args = list(control = list())))
  ))
  stub_fit <- function(formula) suppressWarnings(suppressMessages(
    bmm(formula, oberauer_lewandowsky_2019_e1, model, backend = "mock",
        mock_fit = stub, rename = FALSE)
  ))
  model <- m3(
    resp_cats = c("corr", "other", "dist", "npl"),
    num_options = paste0("n_", c("corr", "other", "dist", "npl")),
    choice_rule = "simple",
    links = list(c = "log", a = "log", d = "log", e = "log")
  )
  f_d <- bmf(corr ~ b + a + c, other ~ b + a, dist ~ b + d, npl ~ b, c ~ 1, a ~ 1, d ~ 1)
  f_e <- bmf(corr ~ b + a + c, other ~ b + a, dist ~ b + e, npl ~ b, c ~ 1, a ~ 1, e ~ 1)

  up <- suppressWarnings(suppressMessages(update(
    stub_fit(f_d), formula. = f_e, backend = "mock", mock_fit = stub,
    rename = FALSE, recompile = TRUE
  )))
  fresh <- stub_fit(f_e)
  expect_equal(brms::stancode(up), brms::stancode(fresh))
  expect_equal(up$bmm$model$parameters, fresh$bmm$model$parameters)
  expect_false("d" %in% names(up$bmm$model$parameters))

  # a constant of the old formula goes with it
  f_dfix <- bmf(corr ~ b + a + c, other ~ b + a, dist ~ b + d, npl ~ b, c ~ 1, a ~ 1, d = 0.5)
  up <- suppressWarnings(suppressMessages(update(
    stub_fit(f_dfix), formula. = f_e, backend = "mock", mock_fit = stub,
    rename = FALSE, recompile = TRUE
  )))
  expect_equal(brms::stancode(up), brms::stancode(fresh))
  expect_false("d" %in% names(up$bmm$model$fixed_parameters))
})

test_that("update.bmmfit needs no newdata for a model that packs its response", {
  skip_on_cran()
  fit <- m3_fixture()
  up <- update_mock(fit)
  expect_true(is(up, "bmmfit"))
  expect_equal(brms::standata(up), brms::standata(fit))

  # a new formula over the predictors brms did keep needs no newdata either.
  # Dropping the cond coefficients of `a` leaves the old fit's global prior for
  # class b unused, which brms warns about whether or not newdata is passed
  up <- suppressWarnings(update_mock(fit, formula. = bmf(
    corr ~ b + a + c, other ~ b + a, dist ~ b + d, npl ~ b,
    c ~ 1 + cond + (1 + cond || ID), a ~ 1 + (1 || ID), d ~ 1 + (1 || ID)
  )))
  expect_true(is(up, "bmmfit"))
  expect_equal(brms::standata(up)$K_a, 1)
})

for (case_name in names(stored_frame_cases())) {
  test_that(paste("the stored model frame of", case_name, "passes check_data again"), {
    skip_on_cran()
    case <- stored_frame_cases()[[case_name]]
    fit <- stored_frame_fit(case)
    warned <- capture_warnings(
      data <- check_stored_data(case$model, fit$data, fit$bmm$user_formula)
    )
    expect_false(any(grepl("reserved", warned)))
    expect_equal(brms::standata(fit, newdata = data), brms::standata(fit))

    model <- fit$bmm$model
    configured <- function(data) suppressWarnings(suppressMessages({
      config <- configure_model(model, data, check_formula(model, data, fit$bmm$user_formula))
      list(
        standata = brms::make_standata(config$formula, config$data, stanvars = config$stanvars),
        prior = configure_prior(model, data, config$formula, user_prior = NULL)
      )
    }))
    expect_equal(
      configured(data),
      configured(suppressWarnings(check_data(model, case$data, fit$bmm$user_formula)))
    )
  })
}

test_that("a set_size rebuilt for check_data() stays out of the model frame", {
  skip_on_cran()
  for (case_name in c("imm", "mixture3p")) {
    case <- stored_frame_cases()[[case_name]]
    fit <- stored_frame_fit(case)
    expect_false("set_size" %in% colnames(fit$data))
    expect_false("set_size" %in% colnames(
      check_stored_data(case$model, fit$data, fit$bmm$user_formula)
    ))
  }
})

test_that("an mpt tree column rebuilt for check_data() stays out of the model frame", {
  skip_on_cran()
  case <- stored_frame_cases()$mpt
  fit <- stored_frame_fit(case)
  expect_false("item_type" %in% colnames(fit$data))
  expect_false("item_type" %in% colnames(
    check_stored_data(case$model, fit$data, fit$bmm$user_formula)
  ))
})

test_that("an m column rebuilt for check_data() stays out of the model frame", {
  skip_on_cran()
  case <- stored_frame_cases()$sdt_mafc_m_column
  fit <- stored_frame_fit(case)
  expect_false("n_afc" %in% colnames(fit$data))
  data <- check_stored_data(case$model, fit$data, fit$bmm$user_formula)
  expect_false("n_afc" %in% colnames(data))
  expect_equal(data$m_afc, as.integer(case$data$n_afc))

  case <- stored_frame_cases()$sdt_ranking_m_column
  fit <- stored_frame_fit(case)
  expect_false("set_size" %in% colnames(fit$data))
  data <- check_stored_data(case$model, fit$data, fit$bmm$user_formula)
  expect_false("set_size" %in% colnames(data))
  expect_equal(data$max_rank, as.numeric(case$data$set_size))
  expect_equal(unname(data$Y), unname(as.matrix(case$data[paste0("rank", 1:5)])))
})

test_that("an m3 frame whose Idx_ columns came from the wrong option columns is refused (#457)", {
  skip_on_cran()
  case <- stored_frame_cases()$m3_char_num_options_by_category
  fit <- stored_frame_fit(case)
  # before #457 each Idx_ column was computed from the option column at the same
  # position of num_options, here (other, corr, npl, dist); zeroing the right
  # columns by those Idx_ columns would rebuild a fit that is wrong in a new way
  idx_cols <- paste0("Idx_", c("corr", "other", "dist", "npl"))
  fit$data[idx_cols] <- fit$data[idx_cols[c(2, 1, 4, 3)]]
  expect_error(
    check_stored_data(case$model, fit$data, fit$bmm$user_formula),
    "refit with `bmm\\(\\)`", ignore.case = TRUE
  )
})

test_that("every column a revert method rebuilds is dropped again", {
  skip_on_cran()
  case <- stored_frame_cases()$mixture3p
  fit <- stored_frame_fit(case)
  local_mocked_s3_method("revert_check_data", "circular", function(model, data) {
    data$extra <- 1
    attr(data, "rebuilt") <- c(attr(data, "rebuilt"), "extra")
    NextMethod("revert_check_data")
  })
  data <- check_stored_data(case$model, fit$data, fit$bmm$user_formula)
  expect_false(any(c("set_size", "extra") %in% colnames(data)))
})

test_that("a user column named LureIdx<n> does not shift the rebuilt set size", {
  skip_on_cran()
  case <- stored_frame_cases()$mixture3p
  set_size <- as.numeric(as.character(case$data$set_size))
  case$data$LureIdx9 <- as.numeric(set_size <= 3)
  case$formula <- bmf(kappa ~ 1 + LureIdx9, thetat ~ 1, thetant ~ 1)
  fit <- stored_frame_fit(case)
  data <- check_stored_data(case$model, fit$data, fit$bmm$user_formula)
  expect_equal(data$ss_numeric, set_size)
})

# brms drops rows whose response is NA, so when every response at the largest
# set size is missing the stored frame never shows that set size (#459)
largest_set_size_cases <- function() {
  cases <- stored_frame_cases()
  nt_features <- paste0("col_nt", 1:7)
  nt_distances <- paste0("dist_nt", 1:7)
  cases <- list(
    mixture3p = cases$mixture3p,
    mixture3p_set_size = cases$mixture3p_set_size,
    imm_full = cases$imm,
    imm_bsc = list(
      model = imm("dev_rad",
        nt_features = nt_features, nt_distances = nt_distances,
        set_size = "set_size", version = "bsc"
      ),
      formula = bmf(c ~ 1, s ~ 1, kappa ~ 1), data = cases$imm$data
    ),
    imm_abc = list(
      model = imm("dev_rad",
        nt_features = nt_features, set_size = "set_size", version = "abc"
      ),
      formula = bmf(c ~ 1, a ~ 1, kappa ~ 1), data = cases$imm$data
    )
  )
  lapply(cases, function(case) {
    case$data$dev_rad[case$data$set_size == 8] <- NA
    case
  })
}

test_that("a stored frame without the largest set size passes check_data again (#459)", {
  skip_on_cran()
  for (case_name in names(largest_set_size_cases())) {
    case <- largest_set_size_cases()[[case_name]]
    fit <- stored_frame_fit(case)
    expect_false(8 %in% fit$data$set_size, label = case_name)
    data <- check_stored_data(case$model, fit$data, fit$bmm$user_formula)
    expect_equal(attr(data, "max_set_size"), 8, label = case_name)
    expect_null(attr(data, "fit_max_set_size"), label = case_name)
    expect_equal(brms::standata(fit, newdata = data), brms::standata(fit), label = case_name)
  }
})

test_that("update() without newdata refits when the largest set size has no response (#459)", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  stub <- methods::new("stanfit", sim = list(
    iter = 10L, warmup = 5L, chains = 1L, thin = 1L,
    samples = list(structure(list(), args = list(control = list())))
  ))
  stub_fit <- function(case) suppressWarnings(suppressMessages(
    bmm(case$formula, case$data, case$model, backend = "mock", mock_fit = stub, rename = FALSE)
  ))
  # a mock fit cannot be reused, so update() takes the recompile path
  stub_update <- function(fit, ...) suppressWarnings(suppressMessages(
    update(fit, ..., backend = "mock", mock_fit = stub, rename = FALSE, recompile = TRUE)
  ))

  for (case_name in c("mixture3p", "imm_full")) {
    fit <- stub_fit(largest_set_size_cases()[[case_name]])
    expect_equal(brms::standata(stub_update(fit)), brms::standata(fit), label = case_name)
  }

  fit <- stub_fit(largest_set_size_cases()$mixture3p_set_size)
  up <- stub_update(fit, formula. = bmf(kappa ~ 1, thetat ~ 1, thetant ~ 1))
  expect_equal(brms::standata(up)$K_thetant, 1)
})

test_that("data without the largest set size still fail the nt_features check", {
  skip_on_cran()
  case <- largest_set_size_cases()$mixture3p
  fit <- stored_frame_fit(case)
  short_data <- case$data[case$data$set_size != 8, ]
  msg <- "'nt_features' should equal max\\(set_size\\)-1"
  expect_error(bmm(case$formula, short_data, case$model, backend = "mock", mock_fit = 1), msg)
  expect_error(update_mock(fit, newdata = short_data), msg)
})

test_that("an updated fit carries nothing from the frame that called update()", {
  skip_on_cran() # over 0.4 s in CRAN mode (ADR 0004, R10)
  skip_if_not_installed("rstan")
  withr::local_seed(445)
  stub <- methods::new("stanfit", sim = list(
    iter = 10L, warmup = 5L, chains = 1L, thin = 1L,
    samples = list(structure(list(), args = list(control = list())))
  ))
  dat <- data.frame(y = rsdm(60, kappa = 5))
  fit <- suppressMessages(bmm(bmf(c ~ 1, kappa ~ 1), dat, sdm("y"),
                              backend = "mock", mock_fit = stub, rename = FALSE))
  # on a recompile, the only route the mock backend runs, brms stores the init
  # function update() builds in the fit, which saveRDS() would then write with
  # the caller's frame
  update_beside <- function(n) {
    ballast <- runif(n)
    suppressMessages(update(fit, backend = "mock", mock_fit = stub, rename = FALSE, recompile = TRUE))
  }
  # each fit is sized on its own, so that the other one's ballast is not counted
  serialized_size <- function(n) length(serialize(update_beside(n), NULL))
  expect_equal(serialized_size(1e6), serialized_size(1))
})

test_that("update() switches an SDT fit between criterion and bias", {
  skip_on_cran() # over 0.4 s in CRAN mode (ADR 0004, R10)
  skip_if_not_installed("rstan")
  stub <- methods::new("stanfit", sim = list(
    iter = 10L, warmup = 5L, chains = 1L, thin = 1L,
    samples = list(structure(list(), args = list(control = list())))
  ))
  dat <- data.frame(n_old = c(10, 40, 12, 38), stimulus = c(0, 1, 0, 1),
                    n_trials = 50)
  fit <- suppressMessages(bmm(bmf(d ~ 1, criterion ~ 1), dat,
                              sdt_yn("n_old", "stimulus", "n_trials"),
                              backend = "mock", mock_fit = stub, rename = FALSE))
  switched <- suppressMessages(update(
    fit, formula. = bmf(d ~ 1, bias ~ 1), backend = "mock", mock_fit = stub,
    rename = FALSE, recompile = TRUE
  ))
  expect_equal(switched$family$name, "sdt_yn_bias")
  expect_true("bias" %in% switched$prior$dpar)
  expect_false("criterion" %in% switched$prior$dpar)

  back <- suppressMessages(update(
    switched, formula. = bmf(d ~ 1, criterion ~ 1), backend = "mock",
    mock_fit = stub, rename = FALSE, recompile = TRUE
  ))
  expect_equal(back$family$name, "sdt_yn")
  expect_false("bias" %in% back$prior$dpar)
})

# =============================================================================
# UPDATING WITH THE STORED MODEL FRAME AS newdata (kfold(), reloo())
# =============================================================================

# models whose stored frame check_data() accepts as it is, so that update() has
# nothing to route
frame_is_user_layout <- c("cswald", "ddm", "ezdm", "mixture2p", "mixture3p_set_size", "sdm")

for (case_name in names(stored_frame_cases())) {
  test_that(paste0("only a slice of a fit's model frame counts as stored data (", case_name, ")"), {
    skip_on_cran()
    case <- stored_frame_cases()[[case_name]]
    fit <- stored_frame_fit(case)
    model <- fit$bmm$model
    slice <- fit$data[-1, , drop = FALSE]
    expect_false(is_stored_frame(model, fit$data, case$data), label = case_name)
    expect_equal(
      is_stored_frame(model, fit$data, slice),
      !case_name %in% frame_is_user_layout,
      label = case_name
    )
    if (case_name %in% frame_is_user_layout) {
      expect_no_error(suppressWarnings(check_data(model, slice, fit$bmm$user_formula)))
    }
  })
}

test_that("user data that hold a column of the stored layout stay user data", {
  skip_on_cran()
  for (case_name in c("m3", "sdt_rating", "sdt_ranking", "sdt_cdp", "mpt")) {
    case <- stored_frame_cases()[[case_name]]
    fit <- stored_frame_fit(case)
    model <- fit$bmm$model
    for (column in c("Y", "nTrials")) {
      user_data <- case$data
      user_data[[column]] <- 1
      expect_false(is_stored_frame(model, fit$data, user_data), label = case_name)
    }
    # all of them, still with the response columns the stored frame has given up
    user_data <- case$data
    user_data[colnames(fit$data)] <- 1
    expect_false(is_stored_frame(model, fit$data, user_data), label = case_name)
  }
})

# a mock fit cannot be reused, so update() takes the recompile path, and brms
# then reads the sampler settings of the fit it updates
stub_stanfit <- function() {
  methods::new("stanfit", sim = list(
    iter = 10L, warmup = 5L, chains = 1L, thin = 1L,
    samples = list(structure(list(), args = list(control = list())))
  ))
}

quietly <- function(x) suppressWarnings(suppressMessages(x))

stub_frame_fit <- function(case) {
  quietly(bmm(case$formula, case$data, case$model,
    backend = "mock", mock_fit = stub_stanfit(), rename = FALSE
  ))
}

stub_update <- function(fit, ...) {
  update(fit, ..., backend = "mock", mock_fit = stub_stanfit(), rename = FALSE, recompile = TRUE)
}

slice_cases <- c("m3", "m3_char_num_options_by_category", "sdt_rating", "sdt_ranking_m_column",
                 "sdt_cdp", "sdt_mafc_m_column", "mixture3p", "imm",
                 "mpt", "mpt_tree_predictor", "mpt_single_tree", "rdm")

# one stub fit per case, shared by the two tests below
stub_frame_fits <- local({
  fits <- list()
  function(case_name) {
    if (is.null(fits[[case_name]])) {
      fits[[case_name]] <<- stub_frame_fit(stored_frame_cases()[[case_name]])
    }
    fits[[case_name]]
  }
})

for (case_name in slice_cases) {
  test_that(paste0("update() takes a slice of the model frame as newdata (kfold(), reloo()) (", case_name, ")"), {
    skip_on_cran()
    skip_if_not_installed("rstan")
    fit <- stub_frame_fits(case_name)
    rows <- seq_len(nrow(fit$data))[-1]
    up <- quietly(stub_update(fit, newdata = fit$data[rows, , drop = FALSE]))
    expect_true(is(up, "bmmfit"), label = case_name)
    expect_equal(nrow(up$data), length(rows), label = case_name)
    expect_equal(brms::standata(up)$N, length(rows), label = case_name)
    expect_equal(attr(up$data, "data_name"), attr(fit$data, "data_name"), label = case_name)
  })
}

for (case_name in slice_cases) {
  test_that(paste0("update() with the whole model frame as newdata is update() without newdata (", case_name, ")"), {
    skip_on_cran()
    skip_if_not_installed("rstan")
    fit <- stub_frame_fits(case_name)
    expect_equal(
      brms::standata(quietly(stub_update(fit, newdata = fit$data))),
      brms::standata(quietly(stub_update(fit))),
      label = case_name
    )
  })
}

test_that("update() still treats data with a column Y as the user's data", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  case <- stored_frame_cases()$sdt_rating
  fit <- stub_frame_fit(case)
  user_data <- case$data
  user_data$Y <- 1
  # only the user-data route warns that Y is reserved: the stored route strips
  # the generated columns before check_data() sees them
  expect_warning(suppressMessages(stub_update(fit, newdata = user_data)), "reserved")
  expect_equal(
    brms::standata(quietly(stub_update(fit, newdata = case$data))),
    brms::standata(fit)
  )
})

# the messages message_displaced_priors() prints, whatever else update() says
displacement_messages <- function(expr) {
  msgs <- capture_messages(expr)
  msgs[grepl("constraint of the model", msgs, fixed = TRUE)]
}

# a fit made under a constraint carries the constraint in its stored prior, so a
# test-only class whose configure_prior method constrains kappa's intercept
# stands in for the shipped models, and the fixture's prior is edited to what
# such a fit stores; only the prior of the update() call is then new. The
# fixture was fitted with silent = 2, which update() restores from the stored
# options and which mutes message2(), so the calls below ask for silent = 1
constrained_fixture <- function() {
  fit <- sdm_fixture()
  class(fit$bmm$model) <- c(class(fit$bmm$model), "constraint_test")
  fit$prior$prior[fit$prior$class == "Intercept" & fit$prior$dpar == "kappa"] <- "constant(5)"
  fit
}

kappa_intercept_prior <- function(fit) {
  fit$prior$prior[fit$prior$class == "Intercept" & fit$prior$dpar == "kappa"]
}

test_that("update() says when a model constraint replaces the prior of the call (#549)", {
  skip_on_cran()
  local_mocked_s3_method(
    "configure_prior", "constraint_test",
    function(model, data, formula, user_prior, ...) {
      brms::prior_("constant(5)", class = "Intercept", dpar = "kappa")
    }
  )
  fit1 <- constrained_fixture()
  update_with <- function(...) update(fit1, ..., silent = 1, backend = "mock", mock_fit = 1, rename = FALSE)

  displaced <- brms::prior_("normal(2, 1)", class = "Intercept", dpar = "kappa")
  msgs <- displacement_messages(fit <- update_with(prior = displaced))
  expect_length(msgs, 1)
  expect_match(msgs, 'dpar = "kappa", class = "Intercept": normal(2, 1) replaced by constant(5)', fixed = TRUE)
  expect_equal(kappa_intercept_prior(fit), "constant(5)")
})

test_that("update() keeps a prior of the call that no constraint touches, silently (#549)", {
  skip_on_cran() # over 0.4 s in CRAN mode (ADR 0004, R10)
  local_mocked_s3_method(
    "configure_prior", "constraint_test",
    function(model, data, formula, user_prior, ...) {
      brms::prior_("constant(5)", class = "Intercept", dpar = "kappa")
    }
  )
  fit1 <- constrained_fixture()
  update_with <- function(...) update(fit1, ..., silent = 1, backend = "mock", mock_fit = 1, rename = FALSE)

  other <- brms::prior_("normal(1, 1)", class = "b", dpar = "c")
  expect_length(displacement_messages(fit <- update_with(prior = other)), 0)
  expect_equal(fit$prior$prior[fit$prior$class == "b" & fit$prior$dpar == "c" & fit$prior$coef == ""], "normal(1, 1)")
  expect_equal(kappa_intercept_prior(fit), "constant(5)")

  # the carried constraint restates the model's own, so nothing is displaced
  expect_length(displacement_messages(fit <- update_with()), 0)
  expect_equal(kappa_intercept_prior(fit), "constant(5)")
})

test_that("a prior of the update() call on a parameter the new formula frees still wins (#549)", {
  skip_on_cran()
  fit1 <- sdm_fixture()
  # mu is fixed at 0 in the fixture; freeing it drops the stored constant(), and
  # the prior of the call must then replace the fresh default, as before the
  # prior of the call went through configure_prior()
  up <- suppressMessages(update(
    fit1, formula. = bmf(c ~ 0 + set_size, kappa ~ 1, mu ~ 1),
    prior = brms::prior_("normal(0, 0.3)", class = "Intercept"),
    backend = "mock", mock_fit = 1, rename = FALSE
  ))
  expect_equal(
    up$prior$prior[up$prior$class == "Intercept" & up$prior$dpar == "" & up$prior$nlpar == ""],
    "normal(0, 0.3)"
  )
  expect_equal(constant_priors(up), 0)
})

# update() of a mock fit reads the sampler's arguments from a stanfit, and a
# fresh mock fit cannot be reused, so these updates recompile
stub_stanfit <- function() {
  methods::new("stanfit", sim = list(
    iter = 10L, warmup = 5L, chains = 1L, thin = 1L,
    samples = list(structure(list(), args = list(control = list())))
  ))
}

user_prior_data <- function() {
  withr::local_seed(445)
  data.frame(y = rmixture2p(n = 60), cond = rep(c("a", "b"), 30))
}

stub_bmm <- function(formula, model, prior = NULL, data = user_prior_data()) {
  suppressWarnings(suppressMessages(bmm(
    formula, data, model, prior = prior, backend = "mock", mock_fit = stub_stanfit(), rename = FALSE
  )))
}

stub_update <- function(fit, ..., data = user_prior_data()) {
  suppressWarnings(suppressMessages(update(
    fit, ..., newdata = data, backend = "mock", mock_fit = stub_stanfit(),
    rename = FALSE, recompile = TRUE
  )))
}

prior_rows <- function(fit) {
  rows <- as.data.frame(fit$prior)[c("prior", "class", "coef", "group", "dpar", "nlpar", "lb", "ub", "source")]
  rows <- rows[do.call(order, c(unname(rows), method = "radix")), ]
  rownames(rows) <- NULL
  rows
}

expect_bmm_prior <- function(up, fresh, label) {
  expect_identical(prior_rows(up), prior_rows(fresh), label = label)
  expect_identical(brms::stancode(up), brms::stancode(fresh), label = label)
}

# the main parameter switches between levels and intercept + effects, which
# changes the default prior of its class; mixture2p's kappa is an nlpar, sdm's
# c a dpar, and brms stores their intercept priors in different classes
user_prior_cases <- list(
  nlpar = list(
    model = mixture2p("y"),
    levels = bmf(kappa ~ 0 + cond, thetat ~ 1),
    effects = bmf(kappa ~ 1 + cond, thetat ~ 1),
    both_predicted = bmf(kappa ~ 1 + cond, thetat ~ 1 + cond),
    user = brms::set_prior("normal(1, 0.5)", class = "b", nlpar = "kappa"),
    other = brms::set_prior("normal(0, 0.5)", class = "b", nlpar = "thetat"),
    intercept = brms::set_prior("normal(0, 0.3)", class = "b", coef = "Intercept", nlpar = "kappa")
  ),
  dpar = list(
    model = sdm("y"),
    levels = bmf(c ~ 0 + cond, kappa ~ 1),
    effects = bmf(c ~ 1 + cond, kappa ~ 1),
    both_predicted = bmf(c ~ 1 + cond, kappa ~ 1 + cond),
    user = brms::set_prior("normal(1, 0.5)", class = "b", dpar = "c"),
    other = brms::set_prior("normal(1, 0.5)", class = "Intercept", dpar = "kappa"),
    intercept = brms::set_prior("normal(0, 0.3)", class = "Intercept", dpar = "c")
  )
)

for (case_name in names(user_prior_cases)) {
  test_that(paste0("update() to a formula gives the priors bmm() gives it (#538, ", case_name, ")"), {
    # eight mock fits, about 8 s: both directions, with and without a user prior
    skip_on_cran()
    skip_if_not_installed("rstan")
    case <- user_prior_cases[[case_name]]
    for (dir in list(c("levels", "effects"), c("effects", "levels"))) {
      for (prior in list(NULL, case$user)) {
        label <- paste(case_name, dir[1], "->", dir[2], if (is.null(prior)) "no prior" else "user prior")
        up <- stub_update(stub_bmm(case[[dir[1]]], case$model, prior), formula. = case[[dir[2]]])
        expect_bmm_prior(up, stub_bmm(case[[dir[2]]], case$model, prior), label)
      }
    }
  })

  test_that(paste0("a chained update() keeps the priors of bmm() and of each call (#538, ", case_name, ")"), {
    skip_on_cran()
    skip_if_not_installed("rstan")
    case <- user_prior_cases[[case_name]]
    both <- combine_prior(case$other, case$intercept)
    up <- stub_update(stub_bmm(case$levels, case$model, case$other),
                      formula. = case$effects, prior = case$intercept)
    expect_bmm_prior(up, stub_bmm(case$effects, case$model, both), "first update")
    expect_identical(up$bmm$user_prior, both)
    up <- stub_update(up, formula. = case$both_predicted)
    expect_bmm_prior(up, stub_bmm(case$both_predicted, case$model, both), "second update")
  })
}

test_that("update() keeps brms from adding back the defaults of the old fit (#538)", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  withr::local_seed(445)
  dat <- data.frame(y = rmixture2p(n = 60), cond = factor(rep(c("a", "b"), 30)), x = rnorm(60),
                    id = rep(1:6, each = 10))
  # brms marks bmm's defaults as user rows and brms::update.brmsfit() adds back
  # the ones the new prior does not restate: the level default of condb, set per
  # coefficient when several predictors have no intercept, and the L row of the
  # correlation next to the new cor row
  cases <- list(
    coefficient = list(model = mixture2p("y"), old = bmf(kappa ~ 0 + cond + x, thetat ~ 1),
                       new = bmf(kappa ~ 1 + cond + x, thetat ~ 1)),
    correlation = list(model = sdm("y"), old = bmf(c ~ 1 + cond + (1 | p | id), kappa ~ 1 + (1 | p | id)),
                       new = bmf(c ~ 0 + cond + (1 | p | id), kappa ~ 1 + (1 | p | id)))
  )
  for (case_name in names(cases)) {
    case <- cases[[case_name]]
    up <- stub_update(stub_bmm(case$old, case$model, data = dat), formula. = case$new, data = dat)
    expect_bmm_prior(up, stub_bmm(case$new, case$model, data = dat), case_name)
  }
})

test_that("a user prior on a parameter the new formula fixes comes back when it is freed (#538)", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  user <- brms::set_prior("normal(2, 1)", class = "Intercept", dpar = "kappa")
  free <- bmf(c ~ 1, kappa ~ 1)
  fixed <- bmf(c ~ 1, kappa = 3)
  fit <- stub_bmm(free, sdm("y"), user)

  # the constant wins and the stored row is named, as bmm() names it for the
  # same formula and prior (#598)
  msgs <- displacement_messages(up <- suppressWarnings(update(
    fit, formula. = fixed, newdata = user_prior_data(), silent = 1,
    backend = "mock", mock_fit = stub_stanfit(), rename = FALSE, recompile = TRUE
  )))
  expect_length(msgs, 1)
  expect_match(msgs, 'dpar = "kappa", class = "Intercept": normal(2, 1) replaced by constant(3) (kappa fixed at 3)',
               fixed = TRUE)
  expect_match(brms::stancode(up), "Intercept_kappa = 3;", fixed = TRUE)
  expect_bmm_prior(up, stub_bmm(fixed, sdm("y"), user), "fixed")

  up <- stub_update(up, formula. = free)
  expect_equal(kappa_intercept_prior(up), "normal(2, 1)")
  expect_bmm_prior(up, stub_bmm(free, sdm("y"), user), "freed")
})

test_that("a user prior on a coefficient an update() drops comes back with it (#538)", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  user <- brms::set_prior("normal(1, 0.5)", class = "b", dpar = "c")
  levels <- bmf(c ~ 0 + cond, kappa ~ 1)
  # bmm() refuses a prior on a coefficient its formula does not have, update()
  # drops it as brms::update.brmsfit() does
  up <- stub_update(stub_bmm(levels, sdm("y"), user), formula. = bmf(c ~ 1, kappa ~ 1))
  expect_false(any(up$prior$class == "b" & up$prior$dpar == "c"))
  up <- stub_update(up, formula. = levels)
  expect_bmm_prior(up, stub_bmm(levels, sdm("y"), user), "restored")
})

test_that("a fit saved before bmm() stored the user's prior keeps the previous update() priors (#538)", {
  skip_on_cran()
  fit1 <- sdm_fixture()
  expect_null(fit1$bmm$user_prior)
  # the bmm default for the levels of c in the old formula is carried as a user
  # prior, because nothing on the fit tells it apart from one
  up <- update_mock(fit1, formula. = bmf(c ~ 1 + set_size, kappa ~ 1))
  expect_equal(up$prior$prior[up$prior$class == "b" & up$prior$dpar == "c" & up$prior$coef == ""],
               "student_t(5, 2, 0.75)")
  # and the updated fit stays on that route when it is updated again
  expect_null(up$bmm$user_prior)
})

test_that("a new fit without a stored user prior updates as a fit saved before (#538)", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  rows <- function(fit) prior_rows(fit)[c("prior", "class", "coef", "dpar", "nlpar", "source")]
  # the expected rows are what update() gave before bmm() stored the user's prior

  fit <- stub_bmm(bmf(kappa ~ 0 + cond, thetat ~ 1), mixture2p("y"))
  fit$bmm$user_prior <- NULL
  up <- stub_update(fit, formula. = bmf(kappa ~ 1 + cond, thetat ~ 1))
  expect_identical(rows(up), data.frame(
    prior = c("", "", "", "constant(-100)", "constant(0)", "constant(0)", "logistic(0, 1)",
              "logistic(0, 1)", "normal(2, 1)"),
    class = c("b", "b", "b", "Intercept", "Intercept", "Intercept", "b", "theta2", "b"),
    coef = c("", "Intercept", "condb", "", "", "", "Intercept", "", ""),
    dpar = c("", "", "", "kappa2", "mu1", "mu2", "", "", ""),
    nlpar = c("thetat", "kappa", "kappa", "", "", "", "thetat", "", "kappa"),
    source = c("default", "default", "default", "user", "user", "user", "user", "default", "user")
  ))
  expect_null(up$bmm$user_prior)

  # a parameter the new formula fixes takes the old route's stale-row filter
  fit <- stub_bmm(bmf(c ~ 1, kappa ~ 1), sdm("y"),
                  brms::set_prior("normal(2, 1)", class = "Intercept", dpar = "kappa"))
  fit$bmm$user_prior <- NULL
  up <- stub_update(fit, formula. = bmf(c ~ 1, kappa = 3))
  expect_identical(rows(up), data.frame(
    prior = c("constant(0)", "constant(3)", "student_t(5, 2, 0.75)"),
    class = c("Intercept", "Intercept", "Intercept"),
    coef = c("", "", ""),
    dpar = c("", "kappa", "c"),
    nlpar = c("", "", ""),
    source = c("user", "user", "user")
  ))
})
