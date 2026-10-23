# newdata for the brms prediction methods (#552). brms checks newdata against
# the model frame, which holds the columns check_data() built; standata() is
# where every prediction method of brms turns newdata into Stan data, so it is
# tested there. The registry's stored-frame data are the layout the user gave
# bmm(), and their rows are the rows of the frame: no case sorts or drops any.

newdata_cases <- local({
  cases <- NULL
  function() {
    if (is.null(cases)) cases <<- model_test_cases("stored_frame")
    cases
  }
})

# one mock fit per case for the whole file; the tests read the fits, so none
# sees another's state
newdata_fit <- local({
  fits <- list()
  function(case_name) {
    if (is.null(fits[[case_name]])) {
      fits[[case_name]] <<- stored_frame_fit(newdata_cases()[[case_name]])
    }
    fits[[case_name]]
  }
})

for (case_name in names(newdata_cases())) {
  test_that(paste("newdata in the layout given to bmm() reaches the Stan data of", case_name), {
    skip_on_cran()
    case <- newdata_cases()[[case_name]]
    fit <- newdata_fit(case_name)
    user <- as.data.frame(case$data)
    rows <- rev(seq_len(nrow(user)))

    expect_equal(
      suppressWarnings(brms::standata(fit, newdata = user)),
      brms::standata(fit)
    )
    # reversed rows come back reversed, not sorted into the frame's order
    expect_equal(
      suppressWarnings(brms::standata(fit, newdata = user[rows, , drop = FALSE])),
      brms::standata(fit, newdata = fit$data[rows, , drop = FALSE])
    )
  })
}

for (case_name in names(newdata_cases())) {
  test_that(paste("slices of the model frame of", case_name, "skip check_data()"), {
    skip_on_cran()
    case <- newdata_cases()[[case_name]]
    fit <- newdata_fit(case_name)
    local_mocked_bindings(check_data = function(model, data, formula) {
      stop("check_data() was called on a slice of the model frame")
    })

    expect_equal(brms::standata(fit, newdata = fit$data), brms::standata(fit))
    half <- seq_len(ceiling(nrow(fit$data) / 2))
    expect_no_error(brms::standata(fit, newdata = fit$data[half, , drop = FALSE]))
  })
}

test_that("newdata for a model that builds no columns skips check_data()", {
  skip_on_cran()
  case_names <- c("sdm", "mixture2p", "ddm", "cswald", "ezdm")
  fits <- lapply(setNames(nm = case_names), newdata_fit)
  local_mocked_bindings(check_data = function(model, data, formula) {
    stop("check_data() was called for a model that builds no columns")
  })
  for (case_name in case_names) {
    case <- newdata_cases()[[case_name]]
    fit <- fits[[case_name]]
    resp <- case$model$resp_vars[[1]]
    expect_equal(brms::standata(fit, newdata = case$data), brms::standata(fit),
                 info = case_name)
    # posterior_predict() needs no response column
    expect_no_error(brms::standata(
      fit, newdata = case$data[setdiff(names(case$data), resp)], check_response = FALSE
    ))
  }
})

# newdata for a prediction can stop below the fit's largest set size, which
# check_data() would otherwise take from the data
test_that("newdata at one set size matches the same rows of the model frame", {
  skip_on_cran()
  cases <- Filter(
    function(case) is_data_var(case$model$other_vars$set_size, case$data),
    newdata_cases()
  )
  expect_gt(length(cases), 0)
  for (case_name in names(cases)) {
    fit <- newdata_fit(case_name)
    user <- as.data.frame(cases[[case_name]]$data)
    set_size <- fit$bmm$model$other_vars$set_size
    rows <- user[[set_size]] == min(as.numeric(as.character(user[[set_size]])))
    expect_no_warning(
      from_user <- brms::standata(fit, newdata = user[rows, , drop = FALSE])
    )
    expect_equal(
      from_user, brms::standata(fit, newdata = fit$data[rows, , drop = FALSE]),
      info = case_name
    )
  }
})

# only the non-target models have a largest set size: one non-target column
# per item beyond the first. mixture2p's capacity rule takes any set size
test_that("newdata with a set size above the fit's maximum is refused", {
  skip_on_cran()
  cases <- Filter(
    function(case) inherits(case$model, "non_targets") &&
      is_data_var(case$model$other_vars$set_size, case$data),
    newdata_cases()
  )
  expect_gt(length(cases), 0)
  for (case_name in names(cases)) {
    fit <- newdata_fit(case_name)
    user <- as.data.frame(cases[[case_name]]$data)
    set_size <- fit$bmm$model$other_vars$set_size
    user[[set_size]] <- as.numeric(as.character(user[[set_size]]))
    above <- length(fit$bmm$model$other_vars$nt_features) + 2
    user[[set_size]][1] <- above
    expect_error(
      brms::standata(fit, newdata = user),
      paste0("set size of ", above, ", above the ", above - 1), fixed = TRUE
    )
  }
})

test_that("check_newdata() leaves no mark on the data it hands to brms", {
  skip_on_cran()
  # any attribute beyond those of the same rows checked on their own is a mark
  # that outlived the check
  for (case_name in names(newdata_cases())) {
    case <- newdata_cases()[[case_name]]
    fit <- newdata_fit(case_name)
    user <- as.data.frame(case$data)
    own <- suppressWarnings(suppressMessages(check_data(fit$bmm$model, user, fit$bmm$user_formula)))
    checked <- suppressWarnings(check_newdata(fit, user))
    expect_identical(
      setdiff(names(attributes(checked)), names(attributes(own))), character(0),
      label = case_name
    )
  }
})

test_that("newdata for mixture2p's capacity versions may hold a larger set size", {
  skip_on_cran()
  fit <- newdata_fit("mixture2p_slot")
  user <- as.data.frame(newdata_cases()$mixture2p_slot$data)
  user$set_size <- as.numeric(as.character(user$set_size))
  user$set_size[1] <- max(user$set_size) + 3
  sdata <- suppressWarnings(brms::standata(fit, newdata = user))
  expect_equal(sdata$vint1[1], user$set_size[1])
})

test_that("newdata without a column check_data() needs stops with its message", {
  skip_on_cran()
  case <- newdata_cases()$sdt_yn
  fit <- newdata_fit("sdt_yn")
  resp <- case$model$resp_vars$response
  expected <- tryCatch(
    check_data(case$model, case$data[setdiff(names(case$data), resp)], fit$bmm$user_formula),
    error = conditionMessage
  )
  expect_error(
    brms::standata(fit, newdata = case$data[setdiff(names(case$data), resp)],
                   check_response = FALSE),
    expected, fixed = TRUE
  )
})

test_that("log_lik() and posterior_predict() of an m3 fit take the data given to bmm()", {
  skip_unless_fixtures()
  fit <- load_fixture_fit("bmmfit_m3_ppcheck.rds")
  user <- oberauer_lewandowsky_2019_e1[oberauer_lewandowsky_2019_e1$ID %in% 1:10, ]

  expect_equal(
    brms::log_lik(fit, newdata = user, draw_ids = 1:5),
    brms::log_lik(fit, draw_ids = 1:5)
  )
  withr::local_seed(1)
  from_user <- brms::posterior_predict(fit, newdata = user, draw_ids = 1:5)
  withr::local_seed(1)
  expect_equal(from_user, brms::posterior_predict(fit, draw_ids = 1:5))
})
