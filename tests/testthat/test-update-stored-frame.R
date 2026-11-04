# Split from test-update.R so that no one file sets the wall time of the
# four-worker run (ADR 0004, 2026-11-03); helpers shared with it are in
# helper-update.R.

# =============================================================================
# UPDATING WITHOUT newdata (#429)
# =============================================================================

m3_fixture <- function() {
  restructure(load_fixture_fit("bmmfit_m3_ppcheck.rds"))
}

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
