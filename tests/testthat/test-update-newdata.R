# Split from test-update.R so that no one file sets the wall time of the
# four-worker run (ADR 0004, 2026-11-03); helpers shared with it are in
# helper-update.R.

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
                 "mpt", "mpt_tree_predictor", "mpt_single_tree", "tcc", "tcc_bins")

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

# =============================================================================
# DATA-SCALED DEFAULT PRIORS AND newdata
# =============================================================================

# messages of an update() whose default priors were scaled to the old data
prior_range_messages <- function(fit, newdata, ...) {
  messages <- character(0)
  withCallingHandlers(
    stub_update(fit, newdata = newdata, ...),
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  grep("priors of the original fit", messages, value = TRUE)
}

test_that("update(newdata =) says once that psychometric priors keep the old intensity range", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  dat <- data.frame(x = rep(0:4, 8), y = rep(c(0L, 0L, 1L, 1L, 1L), 8))
  fit <- stub_frame_fit(list(formula = bmf(midpoint ~ 1, width ~ 1, guess ~ 1, lapse ~ 1),
                             data = dat, model = psychometric_function("y", "x")))
  # range 0 to 4: the message needs a centre shift over 1 or a range outside 2 to 8
  expect_length(prior_range_messages(fit, dat), 0)
  expect_length(prior_range_messages(fit, transform(dat, x = x * 1.1)), 0)
  shifted <- prior_range_messages(fit, transform(dat, x = x + 2))
  expect_length(shifted, 1)
  expect_match(shifted, "formula.", fixed = TRUE)
  expect_match(shifted, "prior =", fixed = TRUE)
  expect_length(prior_range_messages(fit, transform(dat, x = (x - 2) * 3 + 2)), 1)
  expect_length(prior_range_messages(fit, transform(dat, x = x * 3 + 10)), 1)
  # kfold() and reloo() hand update() slices of the model frame, which are not
  # new data: this slice spans 3 to 4 and as user data would get the message
  expect_length(prior_range_messages(fit, fit$data[fit$data$x >= 3, ]), 0)
  expect_length(prior_range_messages(fit, dat[dat$x >= 3, ]), 1)
})

test_that("the psychometric prior-range message names only the parameters left on defaults", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  dat <- data.frame(x = rep(0:4, 8), y = rep(c(0L, 0L, 1L, 1L, 1L), 8))
  shifted <- transform(dat, x = x * 3 + 2)
  formula <- bmf(midpoint ~ 1, width ~ 1, guess ~ 1, lapse ~ 1)
  model <- psychometric_function("y", "x")
  on_midpoint <- brms::set_prior("normal(4, 1)", class = "Intercept", dpar = "midpoint")
  on_width <- brms::set_prior("normal(0, 1)", class = "Intercept", dpar = "width")
  fit <- stub_frame_fit(list(formula = formula, data = dat, model = model))

  message <- prior_range_messages(fit, shifted, prior = on_midpoint)
  expect_length(message, 1)
  expect_match(message, "stay on width.", fixed = TRUE)
  expect_no_match(message, "midpoint", fixed = TRUE)
  expect_length(prior_range_messages(fit, shifted, prior = on_midpoint + on_width), 0)

  # a user prior of the original fit is carried into the update, so it is no default either
  fit <- quietly(bmm(formula, dat, model, prior = on_width,
                     backend = "mock", mock_fit = stub_stanfit(), rename = FALSE))
  expect_match(prior_range_messages(fit, shifted), "stay on midpoint.", fixed = TRUE)
  expect_length(prior_range_messages(fit, shifted, prior = on_midpoint), 0)
})

# midpoint's slopes and SDs are scaled to the intensities too, so user priors on
# its intercept alone leave defaults that go stale
test_that("the psychometric prior-range message counts midpoint's slope and SD defaults", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  dat <- data.frame(x = rep(0:4, 8), y = rep(c(0L, 0L, 1L, 1L, 1L), 8),
                    id = factor(rep(1:4, each = 10)), cond = factor(rep(c("a", "b"), 20)))
  shifted <- transform(dat, x = x * 10 + 20)
  fit <- stub_frame_fit(list(formula = bmf(midpoint ~ 1 + cond + (1 | id), width ~ 1, guess ~ 1, lapse ~ 1),
                             data = dat, model = psychometric_function("y", "x")))
  intercepts <- brms::set_prior("normal(4, 1)", class = "Intercept", dpar = "midpoint") +
    brms::set_prior("normal(0, 1)", class = "Intercept", dpar = "width")
  message <- prior_range_messages(fit, shifted, prior = intercepts)
  expect_length(message, 1)
  expect_match(message, "stay on midpoint.", fixed = TRUE)
  midpoint_rest <- brms::set_prior("normal(0, 10)", class = "b", dpar = "midpoint") +
    brms::set_prior("exponential(0.1)", class = "sd", dpar = "midpoint")
  expect_length(prior_range_messages(fit, shifted, prior = intercepts + midpoint_rest), 0)
})

test_that("the psychometric prior-range message leaves out a parameter the formula fixes", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  dat <- data.frame(x = rep(0:4, 8), y = rep(c(0L, 0L, 1L, 1L, 1L), 8))
  fit <- stub_frame_fit(list(formula = bmf(midpoint = 2, width ~ 1, guess ~ 1, lapse ~ 1),
                             data = dat, model = psychometric_function("y", "x")))
  shifted <- transform(dat, x = x * 3 + 2)
  expect_match(prior_range_messages(fit, shifted), "stay on width.", fixed = TRUE)
  on_width <- brms::set_prior("normal(0, 1)", class = "Intercept", dpar = "width")
  expect_length(prior_range_messages(fit, shifted, prior = on_width), 0)
})

test_that("the psychometric prior range of the log-scale sigmoids is compared on the log axis", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  dat <- data.frame(x = rep(c(0.01, 0.02, 0.04, 0.08, 0.16), 8),
                    y = rep(c(0L, 0L, 1L, 1L, 1L), 8))
  model <- psychometric_function("y", "x", sigmoid = "weibull", guess = 0.5, lapse = 0)
  fit <- stub_frame_fit(list(formula = bmf(midpoint ~ 1, width ~ 1), data = dat, model = model))
  # x 1.5 moves the log centre by 0.41 of a log range of 2.77 (silent), but the
  # linear centre by 0.043 of a range of 0.15 (over a quarter)
  expect_length(prior_range_messages(fit, transform(dat, x = x * 1.5)), 0)
  # x 10 shifts the log axis without stretching it, so width's prior stays as it is
  message <- prior_range_messages(fit, transform(dat, x = x * 10))
  expect_length(message, 1)
  expect_match(message, "stay on midpoint.", fixed = TRUE)
})

# a prior of a class can leave some of its rows on the scaled defaults: the
# per-level main priors of a formula without intercept, or the SDs of a second
# grouping variable
test_that("the psychometric prior-range message names midpoint while any of its scaled rows stay", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  dat <- data.frame(x = rep(0:4, 16), y = rep(c(0L, 0L, 1L, 1L, 1L), 16),
                    id = factor(rep(1:4, each = 20)), item = factor(rep(1:5, 16)),
                    cond = factor(rep(c("a", "b"), 40)),
                    c2 = factor(rep(c("u", "v", "w", "z"), each = 2, length.out = 80)))
  shifted <- transform(dat, x = x * 10 + 20)
  model <- psychometric_function("y", "x")
  on_width <- brms::set_prior("normal(1, 1)", class = "Intercept", dpar = "width")

  fit <- stub_frame_fit(list(formula = bmf(midpoint ~ 0 + cond + c2, width ~ 1, guess ~ 1, lapse ~ 1),
                             data = dat, model = model))
  b_prior <- brms::set_prior("normal(0, 10)", class = "b", dpar = "midpoint")
  expect_match(prior_range_messages(fit, shifted, prior = on_width + b_prior),
               "stay on midpoint.", fixed = TRUE)

  fit <- stub_frame_fit(list(formula = bmf(midpoint ~ 1 + (1 | id) + (1 | item), width ~ 1,
                                           guess ~ 1, lapse ~ 1),
                             data = dat, model = model))
  on_midpoint <- brms::set_prior("normal(25, 5)", class = "Intercept", dpar = "midpoint")
  sd_id <- brms::set_prior("exponential(0.1)", class = "sd", group = "id", dpar = "midpoint")
  expect_match(prior_range_messages(fit, shifted, prior = on_width + on_midpoint + sd_id),
               "stay on midpoint.", fixed = TRUE)
  sd_all <- brms::set_prior("exponential(0.1)", class = "sd", dpar = "midpoint")
  expect_length(prior_range_messages(fit, shifted, prior = on_width + on_midpoint + sd_all), 0)
})
