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

test_that("update() fits a 3par ezcdm fit given a driftangle formula as 4par", {
  skip_on_cran()
  fit <- load_fixture_fit("bmmfit_ezcdm_ppcheck.rds")
  expect_s3_class(fit$bmm$model, "ezcdm_3par")
  # the fixture was fitted with silent = 2, which update() would keep
  expect_message(
    up <- update(fit, formula. = bmf(drift ~ 1, driftangle ~ 1, bound ~ 1, ndt ~ 1),
                 silent = 1, backend = "mock", mock_fit = 1, rename = FALSE),
    "fitted as version '4par'"
  )
  expect_s3_class(up$bmm$model, "ezcdm_4par")
  expect_identical(up$bmm$model$version, "4par")
  expect_identical(attr(up$bmm$model, "call")$version, "4par")
  expect_equal(constant_priors(up), 1)
  expect_match(brms::stancode(up), "real Intercept_driftangle;", fixed = TRUE)
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
# changes the default prior of its class; sdt_rating's d is an nlpar, sdm's
# c a dpar, and brms stores their intercept priors in different classes
user_prior_cases <- list(
  nlpar = list(
    model = sdt_rating(paste0("r", 1:5), "stimulus"),
    data = transform(toy_rating_data(), cond = rep(c("a", "b"), each = 4)),
    levels = bmf(d ~ 0 + cond, criterion ~ 1, spacing ~ 1),
    effects = bmf(d ~ 1 + cond, criterion ~ 1, spacing ~ 1),
    both_predicted = bmf(d ~ 1 + cond, criterion ~ 1 + cond, spacing ~ 1),
    user = brms::set_prior("normal(1, 0.5)", class = "b", nlpar = "d"),
    other = brms::set_prior("normal(0, 0.5)", class = "b", nlpar = "criterion"),
    intercept = brms::set_prior("normal(0, 0.3)", class = "b", coef = "Intercept", nlpar = "d")
  ),
  dpar = list(
    model = sdm("y"),
    data = user_prior_data(),
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
        up <- stub_update(stub_bmm(case[[dir[1]]], case$model, prior, case$data),
                          formula. = case[[dir[2]]], data = case$data)
        expect_bmm_prior(up, stub_bmm(case[[dir[2]]], case$model, prior, case$data), label)
      }
    }
  })

  test_that(paste0("a chained update() keeps the priors of bmm() and of each call (#538, ", case_name, ")"), {
    skip_on_cran()
    skip_if_not_installed("rstan")
    case <- user_prior_cases[[case_name]]
    both <- combine_prior(case$other, case$intercept)
    up <- stub_update(stub_bmm(case$levels, case$model, case$other, case$data),
                      formula. = case$effects, prior = case$intercept, data = case$data)
    expect_bmm_prior(up, stub_bmm(case$effects, case$model, both, case$data), "first update")
    expect_identical(up$bmm$user_prior, both)
    up <- stub_update(up, formula. = case$both_predicted, data = case$data)
    expect_bmm_prior(up, stub_bmm(case$both_predicted, case$model, both, case$data), "second update")
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

# a prior written for bmm 1.4 is translated once, where it enters: translated
# after joining the stored prior, the old row and the stored new-form row on the
# same coefficient both reached brms, which refused them, and the untranslated
# combination was stored and warned again at every update()
test_that("a 1.4-form prior in update() replaces a stored prior on the same coefficient", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  stored <- brms::set_prior("normal(1, 0.5)", class = "b", dpar = "kappa")
  call <- brms::set_prior("normal(0, 0.1)", class = "b", nlpar = "kappa")
  fit <- stub_bmm(bmf(kappa ~ 1 + cond, thetat ~ 1), mixture2p("y"), stored)
  refit <- function(fit, ...) suppressMessages(update(
    fit, ..., newdata = user_prior_data(), backend = "mock", mock_fit = stub_stanfit(),
    rename = FALSE, recompile = TRUE
  ))
  expect_warning(up <- refit(fit, prior = call), "nlpar = \"kappa\"", fixed = TRUE)
  rows <- as.data.frame(up$prior)
  expect_identical(rows$prior[rows$class == "b" & rows$dpar == "kappa" & rows$coef == ""], "normal(0, 0.1)")
  expect_identical(up$bmm$user_prior$nlpar, c("", ""))
  expect_no_warning(refit(up))

  # bmm() stores the translation too
  fit <- stub_bmm(bmf(kappa ~ 1 + cond, thetat ~ 1), mixture2p("y"), call)
  expect_identical(fit$bmm$user_prior$nlpar, c("", ""))
  expect_no_warning(refit(fit))
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
    prior = c("", "constant(0)", "logistic(0, 1)", "normal(2, 1)", "normal(2, 1)"),
    class = c("b", "Intercept", "Intercept", "Intercept", "b"),
    coef = c("condb", "", "", "", ""),
    dpar = c("kappa", "", "thetat", "kappa", "kappa"),
    nlpar = c("", "", "", "", ""),
    source = c("default", "user", "user", "user", "user")
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
