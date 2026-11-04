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
    small_data(),
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
  # 2 participants x 8 set sizes
  expect_equal(nrow(res$cells$counts), 16)
  expect_output(print(res), "No issues detected")
  expect_output(print(res), "Hard checks")
})

test_that("bmm_data_check flags responses in degrees or coded on [0, 2*pi)", {
  model <- mixture3p(
    resp_error = "dev_rad", nt_features = paste0("col_nt", 1:7),
    set_size = "set_size"
  )
  dat_deg <- small_data()
  dat_deg$dev_rad <- dat_deg$dev_rad * 180 / pi
  res_deg <- bmm_data_check(bmf(kappa ~ 1, thetat ~ 1), dat_deg, model)
  expect_true(any(
    finding_severities(res_deg) == "warning" &
      grepl("degrees", finding_messages(res_deg))
  ))
  expect_true(any(grepl("degrees", res_deg$pipeline$warnings)))

  dat_wrap <- small_data()
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
  withr::local_seed(445)
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
  withr::local_seed(445)
  res <- bmm_data_check(
    bmf(kappa ~ nonexistent, thetat ~ 1),
    data.frame(y = runif(10, -3, 3)),
    mixture2p(resp_error = "y")
  )
  expect_true(any(grepl("neither columns", finding_messages(res))))
})

test_that("bmm_data_check does not read a correlation ID as a data column", {
  withr::local_seed(445)
  dat <- data.frame(
    y = runif(40, -3, 3),
    cond = rep(c("a", "b"), each = 20),
    id = rep(1:4, each = 10)
  )
  model <- mixture2p(resp_error = "y")

  res <- bmm_data_check(bmf(kappa ~ cond + (1 |p| id), thetat ~ 1 + (1 |p| id)), dat, model)
  expect_false(any(grepl("neither columns", finding_messages(res))))
  expect_setequal(res$predictors$coding$variable, c("cond", "id"))
  expect_false("p" %in% unlist(res$predictors$pred_map))

  # a correlation ID named like a column is a label, not a predictor
  res_label <- bmm_data_check(
    bmf(kappa ~ 1 + (1 |cond| id), thetat ~ 1 + (1 |cond| id)), dat, model
  )
  expect_setequal(res_label$predictors$coding$variable, "id")
  expect_equal(res_label$cells$cell_vars, character(0))

  # the same name used as a predictor stays one
  res_both <- bmm_data_check(bmf(kappa ~ cond + (1 |cond| id), thetat ~ 1), dat, model)
  expect_setequal(res_both$predictors$coding$variable, c("cond", "id"))
})

test_that("re_group_vars keeps every column a grouping term names", {
  f <- bmf(kappa ~ 1 + (1 |p| gr(id, by = grp)), thetat ~ 1 + (1 | mm(g1, g2)))
  expect_setequal(re_group_vars(f), c("id", "grp", "g1", "g2"))
})

test_that("bmm_data_check does not flag the columns m3 adds to the data (#495)", {
  res <- bmm_data_check(
    bmf(
      corr ~ b + a + c, other ~ b + a, npl ~ b, c ~ 1, a ~ 1,
      dist ~ b + n_opt_dist + nTrials + Idx_dist
    ),
    oberauer_lewandowsky_2019_e1,
    m3(c("corr", "other", "dist", "npl"), c(1, 4, 5, 5), links = list(c = "log", a = "log"))
  )
  expect_null(res$pipeline$error)
  expect_length(res$predictors$unknown_vars, 0)
  expect_false(any(grepl("neither columns", finding_messages(res))))
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

sdt_yn_cells <- function(n_trials) {
  dat <- expand.grid(id = 1:10, cond = c("a", "b"), stimulus = 0:1)
  dat$n_trials <- n_trials
  dat$n_old <- n_trials %/% 2
  dat
}

sdt_yn_check <- function(dat, n_trials = "n_trials") {
  bmm_data_check(
    bmf(d ~ cond + (cond | id), criterion ~ 1),
    dat, sdt_yn("n_old", stimulus = "stimulus", n_trials = n_trials)
  )
}

test_that("bmm_data_check counts the trials of aggregated rows, not the rows (#647)", {
  res <- sdt_yn_check(sdt_yn_cells(50L))
  expect_null(res$pipeline$error)
  expect_identical(res$cells$unit, "trials")
  expect_equal(unique(res$cells$counts$n), 100)
  expect_equal(unique(res$cells$counts$n_rows), 2)
  expect_equal(res$n_trials, 2000)
  expect_false(any(grepl("fewer than", finding_messages(res))))
  out <- capture.output(print(res, color = FALSE))
  expect_true(any(grepl("40 rows, 2000 trials", out)))
  expect_true(any(grepl("Trials per design cell", out)))
  expect_true(any(grepl("trials per cell: min 100 / median 100 / max 100", out)))
})

test_that("bmm_data_check prints large trial totals without scientific notation", {
  out <- capture.output(print(sdt_yn_check(sdt_yn_cells(50000L)), color = FALSE))
  expect_true(any(grepl("40 rows, 2000000 trials", out, fixed = TRUE)))
  expect_true(any(grepl("min 100000 / median 100000 / max 100000", out, fixed = TRUE)))
})

test_that("bmm_data_check flags aggregated cells with fewer than 10 trials in total", {
  res <- sdt_yn_check(sdt_yn_cells(4L))
  msgs <- finding_messages(res)
  expect_true(any(grepl("20 design cell(s) have fewer than 10 trials", msgs, fixed = TRUE)))
  expect_true(any(grepl("id=1, cond=a (n=8)", msgs, fixed = TRUE)))

  res5 <- sdt_yn_check(sdt_yn_cells(5L))
  expect_false(any(grepl("fewer than", finding_messages(res5))))
})

test_that("bmm_data_check flags cells whose rows hold no trials, not as empty cells", {
  dat <- sdt_yn_cells(50L)
  dat[dat$id == 1, c("n_trials", "n_old")] <- NA
  msgs <- finding_messages(sdt_yn_check(dat))
  expect_true(any(grepl("2 design cell(s) have fewer than 10 trials: id=1, cond=a (n=0)",
    msgs, fixed = TRUE)))
  expect_false(any(grepl("contain no trials", msgs)))

  dat[c("n_trials", "n_old")] <- NA_integer_
  res <- sdt_yn_check(dat)
  expect_identical(res$cells$unit, "trials")
  expect_no_warning(out <- capture.output(print(res, color = FALSE)))
  expect_true(any(grepl("id x cond: 20 non-empty cell(s)", out, fixed = TRUE)))
  expect_true(any(grepl("trials per cell: min 0 / median 0 / max 0", out, fixed = TRUE)))
})

test_that("bmm_data_check counts rows and checks no cell when the trial counts are missing", {
  res <- sdt_yn_check(sdt_yn_cells(50L), n_trials = "trials")
  expect_identical(res$cells$unit, "rows")
  expect_null(res$n_trials)
  msgs <- finding_messages(res)
  expect_true(any(grepl("trial counts could not be read", msgs)))
  expect_false(any(grepl("fewer than", msgs)))
  out <- capture.output(print(res, color = FALSE))
  expect_true(any(grepl("Rows per design cell", out)))
  expect_true(any(grepl("Data:  40 rows ('dat')", out, fixed = TRUE)))
})

test_that("bmm_data_check counts m3 trials as the row sums of the response categories", {
  dat <- as.data.frame(oberauer_lewandowsky_2019_e1)
  model <- m3(
    resp_cats = c("corr", "other", "npl"),
    num_options = c("n_corr", "n_other", "n_npl"),
    choice_rule = "simple", version = "ss"
  )
  formula <- bmf(c ~ 1 + (1 | ID), a ~ 1 + (1 | ID))
  res <- bmm_data_check(formula, dat, model)
  expect_null(res$pipeline$error)
  expect_identical(res$cells$unit, "trials")
  expect_equal(
    res$cells$counts$n[res$cells$counts$ID == "1"],
    sum(dat[dat$ID == 1, c("corr", "other", "npl")])
  )
  expect_false(any(grepl("fewer than", finding_messages(res))))

  dat[dat$ID == 1, c("corr", "other", "npl")] <- 1
  sparse <- bmm_data_check(formula, dat, model)
  expect_true(any(grepl("1 design cell(s) have fewer than 10 trials: ID=1 (n=9)",
    finding_messages(sparse), fixed = TRUE)))
})

test_that("bmm_data_check checks each mpt tree when the tree column is in the formula", {
  withr::local_seed(445)
  dat <- mpt_2htm_data(n_id = 4, n_items = 4)
  res <- bmm_data_check(
    bmf(D ~ 0 + item_type + (1 | id), g ~ 1 + (1 | id)), dat,
    mpt(mpt_2htm_trees(), tree_id = "item_type")
  )
  expect_null(res$pipeline$error)
  expect_identical(res$cells$unit, "trials")
  expect_equal(unique(res$cells$counts$n), 4)
  expect_true(any(grepl("8 design cell(s) have fewer than 10 trials",
    finding_messages(res), fixed = TRUE)))
})

test_that("bmm_data_check sums the mpt trees of a cell when no formula uses the tree column", {
  dat <- data.frame(
    id = factor(rep(1:3, each = 2)), item_type = rep(c("old", "new"), 3),
    old = c(30, 1, 28, 0, 31, 1), new = c(10, 1, 12, 2, 9, 1)
  )
  res <- bmm_data_check(
    bmf(D ~ 1 + (1 | id), g ~ 1), dat,
    mpt(mpt_2htm_trees(), tree_id = "item_type")
  )
  expect_equal(res$cells$counts$n, c(42, 42, 42))
  expect_equal(res$cells$counts$n_rows, c(2, 2, 2))
  expect_false(any(grepl("fewer than", finding_messages(res))))
})

test_that("bmm_data_check takes ezdm trials from n_trials", {
  dat <- data.frame(
    id = rep(1:3, each = 2), mrt = 0.5, vrt = 0.02, nup = 3,
    ntr = c(4, 4, 30, 30, 5, 6)
  )
  res <- bmm_data_check(
    bmf(drift ~ 1 + (1 | id), bound ~ 1, ndt ~ 1),
    dat, ezdm(mean_rt = "mrt", var_rt = "vrt", n_upper = "nup", n_trials = "ntr")
  )
  expect_equal(res$cells$counts$n, c(8, 60, 11))
  expect_true(any(grepl("1 design cell(s) have fewer than 10 trials: id=1 (n=8)",
    finding_messages(res), fixed = TRUE)))
})

test_that("data_trials reads the trial counts of each aggregated SDT format", {
  counts <- data.frame(r1 = c(1, 2), r2 = c(3, NA), r3 = c(5, 6), n = c(9, 8))
  expect_equal(data_trials(sdt_mafc("r1", n_trials = "n", m = 4), counts), c(9, 8))
  expect_equal(data_trials(sdt_rating(c("r1", "r2", "r3"), stimulus = "s"), counts), c(9, 8))
  expect_equal(data_trials(sdt_ranking(c("r1", "r2", "r3"), m = 3), counts), c(9, 8))

  cols <- c("new1", "know2", "know3", "remember2", "remember3")
  cdp <- as.data.frame(matrix(1:10, nrow = 2, ncol = 5, dimnames = list(NULL, cols)))
  model <- sdt_cdp(stimulus = "s", n_new = 1, n_old = 2)
  expect_equal(data_trials(model, cdp), c(25, 30))
  cdp$guess2 <- 100
  cdp$guess3 <- 1000
  expect_equal(data_trials(model, cdp), c(1125, 1130))

  colnames(cdp) <- paste0("x_", colnames(cdp))
  cdp$new1 <- 1e6
  prefixed <- sdt_cdp(response = "x_", stimulus = "s", n_new = 1, n_old = 2)
  expect_equal(data_trials(prefixed, cdp), c(1125, 1130))
})

test_that("data_trials reports unreadable counts as NA and trial rows as NULL", {
  counts <- data.frame(r1 = c("1", "2"), n = c(9, 8))
  expect_identical(data_trials(sdt_mafc("r1", n_trials = "absent", m = 4), counts), NA)
  expect_identical(data_trials(sdt_rating(c("r1", "n", "r3"), stimulus = "s"), counts), NA)
  expect_null(data_trials(mixture2p(resp_error = "y"), counts))
})

test_that("bmm_data_check needs nothing beyond check_data and check_formula", {
  withr::local_seed(445)
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
