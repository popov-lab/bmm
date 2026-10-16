save_pars <- brms::save_pars

test_that("update.bmmfit works", {
  skip_on_cran()
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
                 "mpt", "mpt_tree_predictor", "mpt_single_tree")

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
