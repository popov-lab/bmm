test_that("init argument is overwritten if the user supplies it", {
  config_args <- list(formula = "a", family = "b", data = "d", stanvars = "e", init = 1)
  dots <- list(init = 2)
  out <- combine_args(nlist(config_args, dots))
  expect_equal(out, list(formula = "a", family = "b", data = "d", stanvars = "e", init = 2))
})

test_that("user cannot overwrite the custom family", {
  config_args <- list(formula = "a", family = "b", data = "d", stanvars = "e", init = 1)
  dots <- list(family = "c")
  expect_error(combine_args(nlist(config_args, dots)), "Unsupported argument")
})

test_that("empty dots don't crash the function", {
  config_args <- list(formula = "a", family = "b", data = "d", stanvars = "e", init = 1)
  out <- combine_args(nlist(config_args))
  expect_equal(out, list(formula = "a", family = "b", data = "d", stanvars = "e", init = 1))
})

test_that("missing arguments in models are handled correctly", {
  expect_error(mixture2p(), "arguments are missing in mixture2p\\(\\)\\: resp_error")
  expect_error(sdm(), "arguments are missing in sdm\\(\\)\\: resp_error")
  expect_error(mixture3p("y"), "arguments are missing in mixture3p\\(\\)\\: nt_features, set_size")
  expect_error(mixture3p(set_size = "y"), "arguments are missing in mixture3p\\(\\)\\: resp_error, nt_features")
})

test_that("get_variables works", {
  expect_equal(get_variables("a", c("a", "b", "c")), "a")
  expect_equal(get_variables("a", c("a", "b", "c"), regex = TRUE), "a")
  expect_equal(get_variables("a", c("a", "b", "c"), regex = FALSE), "a")
  expect_equal(get_variables("a|b", c("a", "b", "c"), regex = TRUE), c("a", "b"))
  expect_equal(
    get_variables("abc", c("abc1", "abc2", "abc3", "other"), regex = TRUE),
    c("abc1", "abc2", "abc3")
  )
  expect_equal(
    get_variables("^abc", c("abc1", "abc2", "abc3", "other_abc4"), regex = TRUE),
    c("abc1", "abc2", "abc3")
  )
  expect_equal(
    get_variables("abc$", c("nt1_abc", "nt2_abc", "nt3_abc", "other_abc4"), regex = TRUE),
    c("nt1_abc", "nt2_abc", "nt3_abc")
  )
  expect_equal(
    get_variables("nt.*_abc", c("nt1_abc", "nt2_abc", "nt3_abc", "other_abc4"), regex = TRUE),
    c("nt1_abc", "nt2_abc", "nt3_abc")
  )
  expect_equal(get_variables("a|b", c("a", "b", "c"), regex = FALSE), "a|b")
  expect_error(get_variables("d", c("a", "b", "c"), regex = TRUE))
})

test_that("bmm_options works", {
  withr::defer(suppressMessages(bmm_options()))
  expect_message(bmm_options(), "Current bmm options")
  expect_message(bmm_options(sort_data = TRUE), "sort_data = TRUE")
  expect_equal(getOption("bmm.sort_data"), TRUE)
  op <- suppressMessages(bmm_options(sort_data = FALSE))
  expect_equal(getOption("bmm.sort_data"), FALSE)
  options(op)
  expect_equal(getOption("bmm.sort_data"), TRUE)
})


test_that("check_rds_file works", {
  good_files <- list("a.rds", "abc/a.rds", "a", "abc/a", "a.M")
  bad_files <- list(1, mean, c("a", "b"), TRUE)

  for (f in good_files) {
    expect_silent(res <- check_rds_file(f))
    expect_equal(fs::path_ext(res), "rds")
  }

  for (f in bad_files) {
    expect_error(check_rds_file(f))
  }

  expect_null(check_rds_file(NULL))
})

test_that("try_read_bmmfit works", {
  withr::local_options(bmm.sort_data = FALSE)
  mock_fit <- bmm(bmf(c ~ 1, kappa ~ 1), oberauer_lin_2017, sdm("dev_rad"),
    backend = "mock", mock_fit = 1, rename = F
  )
  file <- tempfile()
  mock_fit$file <- paste0(file, ".rds")
  saveRDS(mock_fit, paste0(file, ".rds"))
  expect_equal(try_read_bmmfit(file, FALSE), mock_fit,
    ignore_function_env = TRUE,
    ignore_formula_env = TRUE
  )

  x <- 1
  saveRDS(x, paste0(file, ".rds"))
  expect_error(try_read_bmmfit(file, FALSE), "not of class 'bmmfit'")
})

test_that("try_save_bmmfit works", {
  withr::local_options(bmm.sort_data = FALSE)
  file <- tempfile()
  mock_fit <- bmm(bmf(c ~ 1, kappa ~ 1), oberauer_lin_2017, sdm("dev_rad"),
    backend = "mock", mock_fit = 1, rename = F,
    file = file
  )
  rds_file <- paste0(file, ".rds")
  expect_true(file.exists(rds_file))
  expect_equal(readRDS(rds_file), mock_fit, ignore_function_env = TRUE, ignore_formula_env = TRUE)

  mock_fit2 <- bmm(bmf(c ~ 1, kappa ~ 1), oberauer_lin_2017, sdm("dev_rad"),
    backend = "mock", mock_fit = 2, rename = F,
    file = file
  )
  expect_equal(mock_fit, mock_fit2, ignore_attr = TRUE)

  # they should not be the same if file_refit = TRUE
  mock_fit3 <- bmm(bmf(c ~ 1, kappa ~ 1), oberauer_lin_2017, sdm("dev_rad"),
    backend = "mock", mock_fit = 3, rename = F,
    file = file, file_refit = TRUE
  )
  expect_error(expect_equal(mock_fit, mock_fit3))
})

test_that("is_namedlist works", {
  expect_true(is_namedlist(list(a = 1)))
  expect_true(is_namedlist(list(a = 1, b = 2)))
  expect_true(is_namedlist(nlist()))
  expect_true(is_namedlist(nlist(y ~ 1)))
  arg <- "hello"
  expect_true(is_namedlist(nlist(arg)))

  expect_false(is_namedlist(list(a = 1, 2)))
  expect_false(is_namedlist(list(1, 2)))
  expect_false(is_namedlist(list()))
  expect_false(is_namedlist(list(y ~ 1)))

})
