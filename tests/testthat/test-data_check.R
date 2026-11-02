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

test_that("bmm_data_check does not count mpt rows as trials", {
  withr::local_seed(445)
  res <- bmm_data_check(
    bmf(D ~ 0 + item_type, g ~ 1),
    mpt_2htm_data(n_id = 4),
    mpt(mpt_2htm_trees(), tree_id = "item_type")
  )
  expect_null(res$pipeline$error)
  expect_false(any(grepl("fewer than", finding_messages(res))))
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

test_that("the trial-level report crosses only the participant variable and counts items per level", {
  # 4 subjects x 10 items, subject 4 saw only 5 items
  dat <- expand.grid(id = 1:4, item = 1:10)
  dat <- dat[!(dat$id == 4 & dat$item > 5), ]
  dat$stimulus <- dat$item %% 2
  dat$resp <- (dat$id + dat$item) %% 3 == 0
  res <- bmm_data_check(
    bmf(d ~ 1 + (1 | id) + (1 | item), criterion ~ 1 + (1 | id) + (1 | item)),
    dat, sdt_yn("resp", "stimulus")
  )
  expect_null(res$pipeline$error)
  expect_equal(res$cells$group_vars, "id")
  expect_equal(
    res$cells$other_groups,
    data.frame(variable = "item", n_levels = 10L, min = 3L, median = 3.5, max = 4L)
  )
  expect_equal(
    finding_messages(res),
    "1 design cell(s) have fewer than 10 observations: id=4 (n=5)"
  )
  printed <- capture.output(print(res, color = FALSE))
  expect_match(paste(trimws(printed), collapse = " "),
               "only 'id' (the first in the formula found in the data", fixed = TRUE)
  expect_true(any(grepl(
    "item: 10 levels, observations per level: min 3 / median 3.5 / max 4", printed, fixed = TRUE
  )))
})

# rows of aggregated counts are still counted as observations per cell; the
# trial-level report must leave this report as it was
test_that("the report on aggregated sdt_yn rows is unchanged", {
  dat <- data.frame(id = rep(1:2, each = 2), stimulus = rep(0:1, 2),
                    n_old = c(10, 35, 12, 30), n_trials = 50)
  res <- bmm_data_check(
    bmf(d ~ 1 + (1 | id), criterion ~ 1 + (1 | id)),
    dat, sdt_yn("n_old", "stimulus", "n_trials")
  )
  printed <- capture.output(print(res, color = FALSE))
  cells <- which(printed == "Observations per design cell:")
  expect_identical(printed[cells + 0:4], c(
    "Observations per design cell:",
    "  id: 2 non-empty cell(s)",
    "  observations per cell: min 2 / median 2 / max 2",
    "",
    "Findings:"
  ))
  expect_identical(
    finding_messages(res),
    "2 design cell(s) have fewer than 10 observations: id=1 (n=2); id=2 (n=2)"
  )
})

test_that("a small grouping factor written first does not hide a sparse subject", {
  dat <- expand.grid(id = 1:4, item = 1:20)
  dat <- dat[!(dat$id == 4 & dat$item > 5), ]
  dat$list <- dat$id %% 2
  dat$stimulus <- (dat$item + dat$list) %% 2
  dat$resp <- (dat$id + dat$item) %% 3 == 0
  res <- bmm_data_check(
    bmf(d ~ 1 + (1 | list) + (1 | id) + (1 | item), criterion ~ 1 + (1 | id)),
    dat, sdt_yn("resp", "stimulus")
  )
  expect_equal(res$cells$group_vars, "list")
  expect_equal(finding_messages(res), "1 level(s) of 'id' have fewer than 10 observations: id=4 (n=5)")
})

# fewer items than participants: a rule that counted levels would cross items
# with the condition, report empty item x condition cells, and miss the
# participant with 3 trials in condition b
test_that("the trial-level report crosses the first grouping variable in the formula", {
  dat <- expand.grid(id = 1:60, item = 1:40)
  dat$cond <- factor(ifelse(dat$item <= 20, "a", "b"))
  dat <- dat[!(dat$id == 1 & dat$item > 23), ]
  dat$stimulus <- dat$item %% 2
  dat$resp <- as.integer((dat$id + dat$item) %% 3 == 0)
  res <- bmm_data_check(
    bmf(d ~ cond + (cond | id) + (1 | item), criterion ~ cond + (cond | id) + (1 | item)),
    dat, sdt_yn("resp", "stimulus")
  )
  expect_equal(res$cells$group_vars, "id")
  expect_equal(res$cells$cell_vars, "cond")
  expect_equal(
    finding_messages(res),
    "1 design cell(s) have fewer than 10 observations: id=1, cond=b (n=3)"
  )
  printed <- capture.output(print(res, color = FALSE))
  expect_match(paste(trimws(printed), collapse = " "), "(the first in the formula", fixed = TRUE)
})

# every subject sees a random half of the items: item counts scatter around
# their median without any item being sparse
test_that("items drawn at random for each subject are not flagged", {
  withr::local_seed(4821)
  dat <- do.call(rbind, lapply(1:20, function(i) data.frame(id = i, item = sample(100, 50))))
  dat$stimulus <- dat$item %% 2
  dat$resp <- as.integer((dat$id + dat$item) %% 3 == 0)
  res <- bmm_data_check(
    bmf(d ~ 1 + (1 | id) + (1 | item), criterion ~ 1 + (1 | id) + (1 | item)),
    dat, sdt_yn("resp", "stimulus")
  )
  expect_equal(res$cells$other_groups$median, 10)
  expect_lt(res$cells$other_groups$min, 10)
  expect_identical(finding_messages(res), character(0))
})

test_that("unused factor levels are not counted per level", {
  dat <- expand.grid(id = 1:2, item = factor(1:4, levels = 1:6))
  dat$stimulus <- as.integer(dat$item) %% 2
  dat$resp <- 1L
  res <- bmm_data_check(
    bmf(d ~ 1 + (1 | id) + (1 | item), criterion ~ 1 + (1 | id)),
    dat, sdt_yn("resp", "stimulus")
  )
  expect_equal(
    res$cells$other_groups,
    data.frame(variable = "item", n_levels = 4L, min = 2L, median = 2, max = 2L)
  )
})

# aggregated rows with two grouping variables cross both, as before trial rows
test_that("the report on aggregated rows crosses every grouping variable", {
  dat <- data.frame(id = rep(1:2, each = 2), item = rep(1:2, 2), stimulus = rep(0:1, 2),
                    n_old = c(10, 35, 12, 30), n_trials = 50)
  res <- bmm_data_check(
    bmf(d ~ 1 + (1 | id) + (1 | item), criterion ~ 1 + (1 | id)),
    dat, sdt_yn("n_old", "stimulus", "n_trials")
  )
  expect_equal(res$cells$group_vars, c("id", "item"))
  expect_null(res$cells$other_groups)
})

# subject 5 holds 8 trials and subject 6 holds 5, against a typical 14: the
# bound is min(10, 14 / 2) = 7, so only subject 6 is flagged, and the message
# names 7
test_that("a per-level finding names the bound it applies", {
  dat <- expand.grid(id = 1:6, trial = 1:14)
  dat <- dat[!(dat$id == 5 & dat$trial > 8) & !(dat$id == 6 & dat$trial > 5), ]
  dat$list <- dat$id %% 2
  dat$stimulus <- dat$trial %% 2
  dat$resp <- as.integer((dat$id + dat$trial) %% 3 == 0)
  res <- bmm_data_check(
    bmf(d ~ 1 + (1 | list) + (1 | id), criterion ~ 1 + (1 | id)),
    dat, sdt_yn("resp", "stimulus")
  )
  expect_equal(finding_messages(res), "1 level(s) of 'id' have fewer than 7 observations: id=6 (n=5)")
})

test_that("a grouping variable that is also a fixed predictor stays crossed with the participant", {
  dat <- expand.grid(id = 1:4, item = 1:24)
  dat$grp <- factor(ifelse(dat$item <= 12, "A", "B"))
  dat <- dat[!(dat$id == 1 & dat$item > 15), ]
  dat$stimulus <- dat$item %% 2
  dat$resp <- as.integer((dat$id + dat$item) %% 3 == 0)
  res <- bmm_data_check(
    bmf(d ~ grp + (1 | id) + (1 | grp) + (1 | item), criterion ~ 1 + (1 | id)),
    dat, sdt_yn("resp", "stimulus")
  )
  expect_equal(res$cells$group_vars, "id")
  expect_equal(res$cells$cell_vars, "grp")
  expect_equal(res$cells$other_groups$variable, "item")
  expect_equal(finding_messages(res), "1 design cell(s) have fewer than 10 observations: id=1, grp=B (n=3)")
})

test_that("the per-level bound is printed in full", {
  counts <- list(id = table(rep(1:6, c(2001, 2001, 2002, 2002, 1000, 3))))
  expect_match(group_level_findings(counts, 1500)[[1]]$message,
               "fewer than 1000.5 observations", fixed = TRUE)
})

# a grouping variable crossed as a predictor is not also judged per level, so
# its sparse level C is reported once, as a design cell; only subject 1 saw C
test_that("a grouping variable crossed as a predictor is not counted per level too", {
  dat <- expand.grid(id = 1:4, item = 1:26)
  dat$grp <- factor(ifelse(dat$item <= 13, "A", "B"), levels = c("A", "B", "C"))
  dat <- rbind(dat, data.frame(id = 1, item = 27:28, grp = "C"))
  dat$stimulus <- dat$item %% 2
  dat$resp <- as.integer((dat$id + dat$item) %% 3 == 0)
  res <- bmm_data_check(
    bmf(d ~ grp + (1 | id) + (1 | grp) + (1 | item), criterion ~ 1 + (1 | id)),
    dat, sdt_yn("resp", "stimulus")
  )
  expect_equal(res$cells$other_groups$variable, "item")
  expect_equal(finding_messages(res), c(
    "1 design cell(s) have fewer than 10 observations: id=1, grp=C (n=2)",
    "3 design cell combination(s) contain no observations - fine if the design is intentionally incomplete."
  ))
})

test_that("a grouping variable used as a random slope stays crossed with the participant", {
  dat <- expand.grid(id = 1:4, item = 1:24)
  dat$grp <- factor(ifelse(dat$item <= 12, "A", "B"))
  dat <- dat[!(dat$id == 1 & dat$item > 15), ]
  dat$stimulus <- dat$item %% 2
  dat$resp <- as.integer((dat$id + dat$item) %% 3 == 0)
  res <- bmm_data_check(
    bmf(d ~ 1 + (grp | id) + (1 | grp) + (1 | item), criterion ~ 1 + (1 | id)),
    dat, sdt_yn("resp", "stimulus")
  )
  expect_equal(res$cells$cell_vars, "grp")
  expect_equal(res$cells$other_groups$variable, "item")
  expect_equal(finding_messages(res), "1 design cell(s) have fewer than 10 observations: id=1, grp=B (n=3)")
})
