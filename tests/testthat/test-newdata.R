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

for (case_name in names(newdata_cases())) {
  test_that(paste("newdata in the layout given to bmm() reaches the Stan data of", case_name), {
    skip_on_cran()
    case <- newdata_cases()[[case_name]]
    fit <- stored_frame_fit(case)
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
    fit <- stored_frame_fit(case)
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
  fits <- lapply(newdata_cases()[case_names], stored_frame_fit)
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

test_that("newdata without a column check_data() needs stops with its message", {
  skip_on_cran()
  case <- newdata_cases()$sdt_yn
  fit <- stored_frame_fit(case)
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
