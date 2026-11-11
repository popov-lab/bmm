# Split from test-update.R so that no one file sets the wall time of the
# four-worker run (ADR 0004, 2026-11-03); helpers shared with it are in
# helper-update.R.

# =============================================================================
# UPDATING WITH THE STORED MODEL FRAME AS newdata (kfold(), reloo())
# =============================================================================

# models whose stored frame check_data() accepts as it is, so that update() has
# nothing to route
frame_is_user_layout <- c(
  "cswald", "ddm", "ezdm", "mixture2p", "mixture3p_set_size", "sdm", "sdm_cd"
)

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
