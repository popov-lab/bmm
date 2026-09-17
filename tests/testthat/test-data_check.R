finding_messages <- function(res) {
  vapply(res$findings, `[[`, character(1), "message")
}

finding_severities <- function(res) {
  vapply(res$findings, `[[`, character(1), "severity")
}

test_that("bmm_data_check validates its arguments", {
  expect_error(
    bmm_data_check(kappa ~ 1, data.frame(y = 1), mixture2p(resp_error = "y")),
    "must be a bmmformula"
  )
  expect_error(
    bmm_data_check(bmf(kappa ~ 1), data.frame(y = 1), mixture2p(resp_error = "y"),
      min_trials = 0
    ),
    "must be a single positive number"
  )
  expect_error(
    bmm_data_check(bmf(kappa ~ 1), y ~ 1, mixture2p(resp_error = "y")),
    "coercible to a data.frame"
  )
})

test_that("bmm_data_check produces a clean report for valid data", {
  res <- bmm_data_check(
    bmf(kappa ~ 0 + set_size + (0 + set_size | ID), thetat ~ 0 + set_size),
    oberauer_lin_2017,
    mixture3p(
      resp_error = "dev_rad", nt_features = paste0("col_nt", 1:7),
      set_size = "set_size"
    )
  )
  expect_s3_class(res, "bmm_data_check")
  expect_null(res$pipeline$error)
  expect_length(res$findings, 0)
  expect_equal(res$response$variable, "dev_rad")
  expect_equal(res$response$expected, "radians in [-pi, pi]")
  expect_setequal(res$predictors$pred_map$kappa, c("set_size", "ID"))
  expect_equal(res$predictors$group_vars, "ID")
  expect_equal(nrow(res$cells$counts), 152)
  expect_output(print(res), "No issues detected")
  expect_output(print(res), "Hard checks")
})

test_that("bmm_data_check flags responses in degrees or coded on [0, 2*pi)", {
  model <- mixture3p(
    resp_error = "dev_rad", nt_features = paste0("col_nt", 1:7),
    set_size = "set_size"
  )
  dat_deg <- oberauer_lin_2017
  dat_deg$dev_rad <- dat_deg$dev_rad * 180 / pi
  res_deg <- bmm_data_check(bmf(kappa ~ 1, thetat ~ 1), dat_deg, model)
  expect_true(any(
    finding_severities(res_deg) == "warning" &
      grepl("degrees", finding_messages(res_deg))
  ))
  expect_true(any(grepl("degrees", res_deg$pipeline$warnings)))

  dat_wrap <- oberauer_lin_2017
  dat_wrap$dev_rad <- dat_wrap$dev_rad %% (2 * pi)
  res_wrap <- bmm_data_check(bmf(kappa ~ 1, thetat ~ 1), dat_wrap, model)
  expect_true(any(grepl("wrap\\(\\)", finding_messages(res_wrap))))
  expect_null(res_wrap$pipeline$error)
  expect_length(res_wrap$pipeline$warnings, 0)
})

test_that("bmm_data_check flags misplaced NAs and values in nt_features", {
  dat <- data.frame(
    y = 0, ss = rep(c(1, 2, 3), each = 4),
    nt1 = c(rep(NA, 4), rep(0.5, 8)),
    nt2 = c(rep(NA, 8), rep(0.5, 4))
  )
  model <- mixture3p(resp_error = "y", nt_features = c("nt1", "nt2"), set_size = "ss")
  res_clean <- bmm_data_check(bmf(kappa ~ 1, thetat ~ 1), dat, model)
  expect_length(res_clean$findings, 0)

  dat$nt1[5] <- NA
  dat$nt2[1] <- 0.2
  res <- bmm_data_check(bmf(kappa ~ 1, thetat ~ 1), dat, model)
  msgs <- finding_messages(res)
  sev <- finding_severities(res)
  expect_true(any(sev == "warning" & grepl("NA value\\(s\\) in 'nt_features'", msgs)))
  expect_true(any(sev == "note" & grepl("beyond set_size", msgs)))
})

test_that("bmm_data_check flags misplaced NAs in nt_distances for imm", {
  dat <- data.frame(
    y = 0, nt1 = 0.3, nt2 = -0.3,
    d1 = 1, d2 = c(NA, rep(1, 5))
  )
  model <- imm(
    resp_error = "y", nt_features = c("nt1", "nt2"),
    nt_distances = c("d1", "d2"), set_size = 3, version = "full"
  )
  res <- bmm_data_check(bmf(kappa ~ 1, c ~ 1, a ~ 1, s ~ 1), dat, model)
  expect_true(any(grepl("'nt_distances'", finding_messages(res))))
})

test_that("bmm_data_check captures hard check failures instead of throwing", {
  res <- bmm_data_check(
    bmf(kappa ~ 1, thetat ~ 1),
    data.frame(x = 1:10),
    mixture3p(resp_error = "y", nt_features = "x", set_size = 2)
  )
  expect_match(res$pipeline$error, "not present in the data")
  expect_true(any(grepl("not present", finding_messages(res))))
  expect_output(print(res), "FAILED")
})

test_that("bmm_data_check reports predictor coding and design cells", {
  dat <- data.frame(
    y = runif(40, -3, 3),
    cond = rep(c("a", "b"), each = 20),
    ss_num = rep(c(2, 4), 20),
    id = rep(1:4, each = 10)
  )
  res <- bmm_data_check(
    bmf(kappa ~ cond * ss_num + (1 | id), thetat ~ 1),
    dat, mixture2p(resp_error = "y"),
    min_trials = 8
  )
  coding <- res$predictors$coding
  expect_setequal(coding$variable, c("cond", "ss_num", "id"))
  expect_equal(coding$role[coding$variable == "id"], "grouping")
  msgs <- finding_messages(res)
  expect_true(any(grepl("Character predictor", msgs)))
  expect_true(any(grepl("Numeric predictor", msgs)))
  expect_true(any(grepl("fewer than 8", msgs)))
  expect_true(any(grepl("no observations", msgs)))
})

test_that("bmm_data_check flags formula predictors missing from the data", {
  res <- bmm_data_check(
    bmf(kappa ~ nonexistent, thetat ~ 1),
    data.frame(y = runif(10, -3, 3)),
    mixture2p(resp_error = "y")
  )
  expect_true(any(grepl("neither columns", finding_messages(res))))
})

test_that("data_check_findings returns an empty list for models without methods", {
  expect_identical(
    data_check_findings(
      ezdm(mean_rt = "mrt", var_rt = "vrt", n_upper = "nup", n_trials = "ntr"),
      data.frame(mrt = 0.5, vrt = 0.02, nup = 40, ntr = 50),
      bmf(drift ~ 1)
    ),
    list()
  )
})

test_that("bmm_data_check skips the min_trials note for aggregate-data models", {
  res <- bmm_data_check(
    bmf(c ~ 1 + (1 | ID), a ~ 1 + (1 | ID)),
    oberauer_lewandowsky_2019_e1,
    m3(
      resp_cats = c("corr", "other", "npl"),
      num_options = c("n_corr", "n_other", "n_npl"),
      choice_rule = "simple", version = "ss"
    )
  )
  expect_null(res$pipeline$error)
  expect_false(any(grepl("fewer than", finding_messages(res))))
  expect_equal(min(res$cells$counts$n), 3)
})

test_that("bmm_data_check needs nothing beyond check_data and check_formula", {
  dat <- data.frame(y = runif(40, -3, 3), cond = rep(c("a", "b"), 20))
  model <- sdm(resp_error = "y")
  class(model) <- setdiff(class(model), "circular")

  expect_length(response_annotations(model), 0)
  expect_identical(data_check_findings(model, dat, bmf(c ~ cond, kappa ~ 1)), list())

  res <- bmm_data_check(bmf(c ~ cond, kappa ~ 1), dat, model)
  expect_null(res$pipeline$error)
  expect_true(is.na(res$response$expected))
  expect_setequal(res$predictors$coding$variable, "cond")
  expect_output(print(res), "Response variables")
  expect_output(print(res), "Hard checks")
})

test_that("bmm_data_check captures check_model failures instead of throwing", {
  dat <- data.frame(
    corr = c(5, 6, 7, 4), other = c(2, 1, 0, 3), npl = 1,
    n_corr = 1, n_other = 2, n_npl = 3, ID = 1:4
  )
  model <- m3(
    resp_cats = c("corr", "other", "npl"),
    num_options = c("n_corr", "n_other", "n_npl"), choice_rule = "simple"
  )
  res <- bmm_data_check(bmf(b ~ 1, c ~ 1, a ~ 1), dat, model)
  expect_s3_class(res, "bmm_data_check")
  expect_match(res$pipeline$error, "link function")
  expect_output(print(res), "FAILED")

  res_class <- bmm_data_check(bmf(kappa ~ 1), data.frame(y = 1:5), list(foo = "bar"))
  expect_s3_class(res_class, "bmm_data_check")
  expect_match(res_class$pipeline$error, "model argument")
  expect_output(print(res_class), "FAILED")
})

test_that("bmm_data_check flags misplaced NAs in nt_distances for imm_bsc", {
  dat <- data.frame(
    y = 0, nt1 = 0.3, nt2 = -0.3,
    d1 = 1, d2 = c(NA, rep(1, 5))
  )
  model <- imm(
    resp_error = "y", nt_features = c("nt1", "nt2"),
    nt_distances = c("d1", "d2"), set_size = 3, version = "bsc"
  )
  res <- bmm_data_check(bmf(kappa ~ 1, c ~ 1, s ~ 1), dat, model)
  expect_true(any(grepl("'nt_distances'", finding_messages(res))))
})

# the Model: header prints the model call on one line, as print.bmmodel does,
# so this guards the sections that describe the data
test_that("bmm_data_check wraps every printed line at 80 characters", {
  dat <- data.frame(
    resp_error = seq(-3, 3, length.out = 400),
    cond = factor(rep(paste0("condition_level_", 1:20), each = 20)),
    ID = rep(1:20, each = 20)
  )
  res <- bmm_data_check(
    bmf(c ~ 0 + cond, kappa ~ 1 + (1 | ID)), dat, sdm(resp_error = "resp_error")
  )
  expect_true(all(nchar(capture.output(print(res, color = FALSE))) <= 80))
})

test_that("bmm_data_check treats design-cell crossing and value display alike", {
  dat <- data.frame(y = seq(-3, 3, length.out = 80), cond = rep(1:8, 10), ID = rep(1:4, 20))
  res <- bmm_data_check(
    bmf(kappa ~ 0 + cond, thetat ~ 1 + (1 | ID)), dat, mixture2p(resp_error = "y")
  )
  expect_true("cond" %in% res$cells$cell_vars)
  coding <- res$predictors$coding
  expect_match(coding$summary[coding$variable == "cond"], "8 unique values")
  expect_true(any(grepl("Numeric predictor", finding_messages(res))))
})

test_that("bmm_data_check reports a skipped design-cell table instead of omitting it", {
  dat <- data.frame(
    y = seq(-3, 3, length.out = 400),
    cond = factor(rep(paste0("c", 1:300), length.out = 400)),
    ID = 1:400
  )
  res <- bmm_data_check(
    bmf(kappa ~ 0 + cond, thetat ~ 1 + (1 | ID)), dat, mixture2p(resp_error = "y")
  )
  expect_null(res$cells$counts)
  expect_equal(res$cells$n_combos, 120000)
  expect_output(print(res), "not tabulated")
})

test_that("data_check_findings dispatches to the default method for RT models", {
  dat <- data.frame(rt = 0.5, resp = 1)
  expect_identical(
    data_check_findings(ddm(rt = "rt", response = "resp"), dat, bmf(drift ~ 1)),
    list()
  )
  expect_identical(
    data_check_findings(cswald(rt = "rt", response = "resp"), dat, bmf(drift ~ 1)),
    list()
  )
})

test_that("data_check_finding rejects an unrecognized severity", {
  expect_error(data_check_finding("warnings", "text"), "should be one of")
  expect_equal(data_check_finding("note", "text")$severity, "note")
})
