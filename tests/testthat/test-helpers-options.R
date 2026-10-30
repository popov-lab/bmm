# The option-set contract of sdm() (#558, D3, D4): a response grid of
# n_options equally spaced options, and per-trial sets of offered options

# five trials with sets of different size and spacing; the target (0) is in
# every set and each response is one of its trial's options, in degrees
options_data <- function() {
  data.frame(
    y = deg2rad(c(0, 90, -45, 0, 180)),
    o1 = deg2rad(c(0, 0, 0, -90, 0)),
    o2 = deg2rad(c(180, 90, -45, 0, 180)),
    o3 = deg2rad(c(NA, 180, 45, 90, NA)),
    o4 = deg2rad(c(NA, -90, NA, 180, NA)),
    m = c(2, 4, 3, 4, 2)
  )
}

check_options <- function(model, data, formula = bmf(c ~ 1, kappa ~ 1)) {
  check_data(model, data, formula)
}

test_that(".bin_options() puts each response on the nearest option, -pi and pi on one", {
  expect_identical(
    .bin_options(c(-pi, pi, 0, 0.0087, 0.0088, -0.0088, 3.1329), 360, "radians"),
    c(180L, 180L, 0L, 0L, 1L, -1L, 180L)
  )
  expect_identical(.bin_options(c(-180, 180, 181, -179, 0, 359), 360, "bins"), c(180L, 180L, -179L, -179L, 0L, -1L))
  expect_identical(.bin_options(c(NA, 1), 4, "bins"), c(NA, 1L))
})

test_that("sdm() checks the option arguments", {
  expect_error(sdm("y", n_options = 7), "even whole number of at least 4")
  expect_error(sdm("y", n_options = 2), "even whole number of at least 4")
  expect_error(sdm("y", n_options = c(360, 180)), "even whole number of at least 4")
  expect_error(sdm("y", n_options = 360.5), "even whole number of at least 4")
  expect_error(sdm("y", options = c("a", "b")), "need `n_options`")
  expect_error(sdm("y", n_alternatives = 4), "need `n_options`")
  expect_error(sdm("y", units = "bins"), "need `n_options`")
  expect_error(sdm("y", n_options = 360, n_alternatives = "m"), "needs the `options` columns")
  expect_error(sdm("y", n_options = 360, n_alternatives = 4, options = c("a", "b")), "not both")
  expect_error(sdm("y", n_options = 360, n_alternatives = 7), "multiple of 7")
  expect_error(sdm("y", n_options = 360, n_alternatives = 1), "whole number of at least 2")
  expect_error(sdm("y", n_options = 360, n_alternatives = c(2, 4)), "whole number of at least 2")
  expect_error(sdm("y", n_options = 360, options = 1:2), "column names")
  expect_error(sdm("y", n_options = 360, units = "steps"), "should be one of")
})

test_that("only a model with n_options has the option_set class and other_vars", {
  expect_identical(class(sdm("y")), c("bmmodel", "circular", "sdm", "sdm_simple"))
  expect_identical(sdm("y")$other_vars, list(activation = "vonmises", scaling = "density"))
  grid <- sdm("y", n_options = 360)
  expect_identical(class(grid), c("bmmodel", "circular", "option_set", "sdm", "sdm_simple"))
  expect_identical(
    grid$other_vars,
    list(activation = "vonmises", scaling = "density", n_options = 360L, units = "radians")
  )
  sets <- sdm("y", n_options = 360, options = c("o1", "o2"), n_alternatives = "m", units = "bins")
  expect_identical(sets$other_vars$options, c("o1", "o2"))
  expect_identical(sets$other_vars$n_alternatives, "m")
  expect_identical(sets$other_vars$units, "bins")
  expect_identical(sdm("y", n_options = 360, n_alternatives = 4)$other_vars$n_alternatives, 4L)
})

test_that("the full wheel bins the response and builds no option columns", {
  dat <- data.frame(y = c(-pi, 0.01, 0.0087, 2))
  checked <- check_options(sdm("y", n_options = 360), dat)
  expect_identical(checked$sdm_bin, c(180L, 1L, 0L, 115L))
  expect_false(any(grepl("^OptIdx", names(checked))))
  expect_null(checked$opt_set_id)
  expect_identical(built_data_columns(sdm("y", n_options = 360)), "sdm_bin")
  expect_identical(built_data_columns(sdm("y")), character(0))
})

test_that("option columns become integer offsets and a set id from the sorted offsets", {
  model <- sdm("y", n_options = 360, options = paste0("o", 1:4), n_alternatives = "m")
  checked <- check_options(model, options_data())
  absent <- .option_absent(360)
  expect_identical(checked$sdm_bin, c(0L, 90L, -45L, 0L, 180L))
  expect_identical(checked$OptIdx1, c(0L, 0L, 0L, -90L, 0L))
  expect_identical(checked$OptIdx2, c(180L, 90L, -45L, 0L, 180L))
  expect_identical(checked$OptIdx3, c(absent, 180L, 45L, 90L, absent))
  expect_identical(checked$OptIdx4, c(absent, -90L, absent, 180L, absent))
  # trials 2 and 4 offer the same set in another order, as do 1 and 5
  expect_identical(checked$opt_set_id, c(1L, 2L, 3L, 2L, 1L))
  expect_setequal(built_data_columns(model), c("sdm_bin", paste0("OptIdx", 1:4), "opt_set_id"))
})

test_that("NA padding may sit anywhere in a row, and the n_alternatives column is optional", {
  dat <- options_data()
  dat[3, c("o2", "o3")] <- dat[3, c("o3", "o2")]
  with_m <- check_options(sdm("y", n_options = 360, options = paste0("o", 1:4), n_alternatives = "m"), dat)
  without <- check_options(sdm("y", n_options = 360, options = paste0("o", 1:4)), dat)
  expect_identical(with_m$opt_set_id, without$opt_set_id)
  expect_identical(with_m[paste0("OptIdx", 1:4)], without[paste0("OptIdx", 1:4)])
  expect_identical(with_m$OptIdx2[3], 45L)
})

test_that("the scalar n_alternatives is the equally spaced set written out as columns", {
  dat <- data.frame(y = deg2rad(c(0, 90, -90, 180, 0)), cond = c("a", "a", "b", "b", "b"))
  explicit <- dat
  explicit[paste0("o", 1:4)] <- matrix(deg2rad(c(0, 90, 180, 270)), nrow(dat), 4, byrow = TRUE)
  formula <- bmf(c ~ cond, kappa ~ 1)
  short <- sdm("y", n_options = 360, n_alternatives = 4)
  long <- sdm("y", n_options = 360, options = paste0("o", 1:4))
  checked_short <- check_options(short, dat, formula)
  checked_long <- check_options(long, explicit, formula)
  built <- c("sdm_bin", paste0("OptIdx", 1:4), "opt_set_id")
  expect_identical(checked_short[built], checked_long[built])
  expect_identical(built_data_columns(short), built)
  expect_identical(standata(formula, dat, short), standata(formula, explicit, long))
})

test_that("the option checks refuse what the likelihood cannot take", {
  dat <- options_data()
  model <- sdm("y", n_options = 360, options = paste0("o", 1:4), n_alternatives = "m")

  no_target <- dat
  no_target$o1[2] <- deg2rad(-135)
  expect_error(check_options(model, no_target), "target \\(an offset of 0\\) is not among the options in row 2")

  not_offered <- dat
  not_offered$y[3] <- deg2rad(10)
  expect_error(check_options(model, not_offered), "response is not one of the options in row 3")

  off_grid <- dat
  off_grid$o2[1] <- deg2rad(180.4)
  expect_error(check_options(model, off_grid), "row 1 .* 0.4 of a step from the nearest of the 360 positions")

  duplicate <- dat
  duplicate$o3[1] <- deg2rad(-180)
  duplicate$m[1] <- 3
  expect_error(check_options(model, duplicate), "offers the same option twice in row 1")

  wrong_m <- dat
  wrong_m$m[4] <- 3
  expect_error(check_options(model, wrong_m), "'m' says 3 options in row 4, but the option columns hold 4")

  na_m <- dat
  na_m$m[2] <- NA
  expect_error(check_options(model, na_m), "'m' is missing in row 2")

  degrees <- dat
  degrees$o2 <- rad2deg(degrees$o2)
  expect_error(check_options(model, degrees), "must be in radians")

  missing_col <- dat
  missing_col$o4 <- NULL
  expect_error(check_options(model, missing_col), "'o4' is not in the data")
})

test_that("options a little off the grid are rounded onto it with a message", {
  dat <- options_data()
  dat$o2 <- round(dat$o2, 3)
  dat$o3 <- round(dat$o3, 3)
  dat$y <- round(dat$y, 3)
  model <- sdm("y", n_options = 360, options = paste0("o", 1:4), n_alternatives = "m")
  expect_message(checked <- check_options(model, dat), "rounded onto the 360 positions")
  expect_identical(checked$OptIdx2, c(180L, 90L, -45L, 0L, 180L))
  # floating noise from deg2rad() is no reason to say so
  expect_no_message(check_options(model, options_data()))
})

test_that("rows with a missing response are not checked but keep their place", {
  dat <- options_data()
  dat$y[2] <- NA
  dat$o1[2] <- deg2rad(45)
  model <- sdm("y", n_options = 360, options = paste0("o", 1:4), n_alternatives = "m")
  checked <- check_options(model, dat)
  expect_identical(nrow(checked), nrow(dat))
  expect_true(is.na(checked$sdm_bin[2]))
})

test_that("units = 'bins' takes whole option steps for the response and the options", {
  dat <- data.frame(y = c(0, 90, 270), o1 = 0, o2 = c(180, 90, -90))
  model <- sdm("y", n_options = 360, options = c("o1", "o2"), units = "bins")
  checked <- check_options(model, dat)
  expect_identical(checked$sdm_bin, c(0L, 90L, -90L))
  expect_identical(checked$OptIdx2, c(180L, 90L, -90L))
  dat$o2[1] <- 180.5
  expect_error(check_options(model, dat), "whole numbers")
  expect_error(check_options(sdm("y", n_options = 360, units = "bins"), data.frame(y = 1.5)), "whole numbers")
  # option steps beyond 2 pi are no sign of degrees here
  expect_no_warning(check_options(sdm("y", n_options = 360, units = "bins"), data.frame(y = c(0, 100))))
})

test_that("a normaliser cell is a run of equal parameters times an option set", {
  dat <- options_data()
  dat$cond <- c("a", "a", "a", "b", "b")
  model <- sdm("y", n_options = 360, options = paste0("o", 1:4), n_alternatives = "m")
  sdata <- suppressMessages(standata(bmf(c ~ cond, kappa ~ 1), dat, model))
  # runs a = rows 1-3 (sets 1, 2, 3) and b = rows 4-5 (sets 2, 1)
  expect_identical(as.integer(sdata$G_sdm_cells), 5L)
  expect_identical(as.integer(sdata$sdm_cell_first), 1:5)
  expect_identical(as.integer(sdata$sdm_cell_count), rep(1L, 5))
  expect_identical(as.integer(sdata$sdm_cell_set), c(1L, 2L, 3L, 2L, 1L))
  expect_identical(as.integer(sdata$sdm_set_size), c(2L, 4L, 3L))
  expect_identical(as.integer(sdata$sdm_set_off[2, ]), c(-90L, 0L, 90L, 180L))
  expect_identical(as.integer(sdata$sdm_K), 360L)

  # one run: the two trials with set 1 share a cell, wherever they are
  sdata <- suppressMessages(standata(bmf(c ~ 1, kappa ~ 1), dat, model))
  expect_identical(as.integer(sdata$sdm_cell_first), 1:3)
  expect_identical(as.integer(sdata$sdm_cell_count), c(2L, 2L, 1L))

  # the full wheel: one cell per run, set 0
  sdata <- standata(bmf(c ~ cond, kappa ~ 1), dat, sdm("y", n_options = 360))
  expect_identical(as.integer(sdata$sdm_cell_first), c(1L, 4L))
  expect_identical(as.integer(sdata$sdm_cell_count), c(3L, 2L))
  expect_identical(as.integer(sdata$sdm_cell_set), c(0L, 0L))
})

test_that("cells follow the rows brms keeps after dropping missing values", {
  dat <- options_data()
  dat$y[2] <- NA
  model <- sdm("y", n_options = 360, options = paste0("o", 1:4), n_alternatives = "m")
  # brms warns that it dropped the row
  sdata <- suppressWarnings(suppressMessages(standata(bmf(c ~ 1, kappa ~ 1), dat, model)))
  expect_identical(as.integer(sdata$N), 4L)
  expect_identical(as.integer(sdata$sdm_cell_first), c(1L, 2L, 3L))
  expect_identical(as.integer(sdata$sdm_cell_count), c(2L, 1L, 1L))
  expect_identical(as.integer(sdata$sdm_set_size), c(2L, 3L, 4L))
})

test_that("more normaliser cells than half the rows give a message", {
  dat <- data.frame(y = deg2rad(rep(c(0, 90), 10)), x = seq(0, 1, length.out = 20))
  expect_message(
    standata(bmf(c ~ x, kappa ~ 1), dat, sdm("y", n_options = 360)),
    "20 of the 20 rows"
  )
  expect_no_message(standata(bmf(c ~ 1, kappa ~ 1), dat, sdm("y", n_options = 360)))
})

test_that("the stored frame gives the response and the options back at the option centres", {
  dat <- options_data()
  dat$y[1] <- deg2rad(0.3)
  model <- sdm("y", n_options = 360, options = paste0("o", 1:4), n_alternatives = "m")
  checked <- check_options(model, dat)
  stored <- checked[c("sdm_bin", paste0("OptIdx", 1:4))]
  reverted <- revert_check_data(model, stored)
  expect_equal(reverted$y, deg2rad(c(0, 90, -45, 0, 180)))
  expect_equal(reverted$o3, deg2rad(c(NA, 180, 45, 90, NA)))
  expect_identical(reverted$m, c(2L, 4L, 3L, 4L, 2L))
  expect_false(any(c("sdm_bin", "OptIdx1", "opt_set_id") %in% names(reverted)))
  expect_setequal(attr(reverted, "rebuilt"), c("y", paste0("o", 1:4), "m"))
  rechecked <- check_options(model, reverted)
  expect_identical(rechecked[names(stored)], checked[names(stored)])

  bins <- sdm("y", n_options = 360, units = "bins")
  expect_identical(revert_check_data(bins, data.frame(sdm_bin = c(3L, -2L)))$y, c(3L, -2L))
})
