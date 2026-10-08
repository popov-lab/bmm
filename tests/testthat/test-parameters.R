library(testthat)

# ===========================================================================
# parameter_info() for bmmodel objects
# ===========================================================================

test_that("parameter_info() returns correct structure for sdm", {
  m <- sdm(resp_error = "y")
  p <- parameter_info(m)

  expect_s3_class(p, "bmm_parameters")
  expect_s3_class(p, "data.frame")
  expect_true(all(c("parameter", "description", "fixed", "value", "link") %in% names(p)))
  expect_equal(nrow(p), 3)
  expect_true(all(c("mu", "c", "kappa") %in% p$parameter))
})

test_that("parameter_info() flags fixed parameters for sdm", {
  m <- sdm(resp_error = "y")
  p <- parameter_info(m)

  mu_row <- p[p$parameter == "mu", ]
  expect_true(mu_row$fixed)
  expect_equal(mu_row$value, "0")

  c_row <- p[p$parameter == "c", ]
  expect_false(c_row$fixed)
  expect_true(is.na(c_row$value))
})

test_that("parameter_info() shows correct link functions", {
  m <- sdm(resp_error = "y")
  p <- parameter_info(m)

  expect_equal(p$link[p$parameter == "c"], "log")
  expect_equal(p$link[p$parameter == "kappa"], "log")
  expect_equal(p$link[p$parameter == "mu"], "tan_half")
})

test_that("parameter_info() works for imm model", {
  m <- imm(
    resp_error = "y", nt_features = "nt",
    nt_distances = "d", set_size = 2
  )
  p <- parameter_info(m)

  expect_s3_class(p, "bmm_parameters")
  expect_true(all(c("kappa", "a", "c", "s") %in% p$parameter))
})

test_that("parameter_info() works for mixture3p model", {
  m <- mixture3p(resp_error = "y", nt_features = "nt", set_size = 2)
  p <- parameter_info(m)

  expect_s3_class(p, "bmm_parameters")
  expect_true(all(c("thetat", "thetant", "kappa") %in% p$parameter))
})

test_that("parameter_info() works for m3 ss version", {
  m <- m3(
    resp_cats = c("corr", "other", "npl"),
    num_options = c(1, 4, 5), version = "ss"
  )
  p <- parameter_info(m)

  expect_s3_class(p, "bmm_parameters")
  expect_true(all(c("b", "c", "a") %in% p$parameter))
})

test_that("parameter_info() for m3 custom without formula shows note", {
  m <- m3(
    resp_cats = c("corr", "other", "npl"),
    num_options = c(1, 4, 5), version = "custom",
    links = list(c = "log", a = "log")
  )
  p <- parameter_info(m)

  expect_false(is.null(attr(p, "m3_note")))
  expect_match(attr(p, "m3_note"), "custom M3 model")
})

test_that("parameter_info() for m3 custom with formula discovers params", {
  m <- m3(
    resp_cats = c("corr", "other", "npl"),
    num_options = c(1, 4, 5), version = "custom",
    links = list(c = "log", a = "log")
  )
  ff <- bmf(corr ~ b + a + c, other ~ b + a, npl ~ b, c ~ 1, a ~ 1)
  p <- parameter_info(m, formula = ff)

  expect_true("a" %in% p$parameter)
  expect_true("c" %in% p$parameter)
  expect_null(attr(p, "m3_note"))
})

test_that("parameter_info() for m3 custom lists data columns and counts not as parameters (#495)", {
  m <- m3(
    resp_cats = c("corr", "other", "dist", "npl"),
    num_options = c(1, 4, 5, 5),
    links = list(c = "log", a = "log", x = "log")
  )
  ff <- bmf(
    corr ~ b + a + c * time, other ~ b + a, dist ~ b + x * n_opt_dist, npl ~ b,
    c ~ 1, a ~ 1, x ~ 1
  )
  p <- parameter_info(m, formula = ff)

  expect_true("x" %in% p$parameter)
  expect_false(any(c("time", "n_opt_dist") %in% p$parameter))
  checked <- suppressWarnings(check_model(m, data.frame(time = 1), ff))
  expect_setequal(p$parameter, names(checked$parameters))
})

test_that("parameter_info() identifies free parameters for sdm", {
  m <- sdm(resp_error = "y")
  p <- parameter_info(m)

  free_pars <- p[!p$fixed, ]
  expect_equal(nrow(free_pars), 2)
  expect_true(all(c("c", "kappa") %in% free_pars$parameter))
})

test_that("parameter_info() includes descriptions", {
  m <- sdm(resp_error = "y")
  p <- parameter_info(m)

  expect_true(all(nchar(p$description) > 0))
  expect_true(is.character(p$description))
})

test_that("parameter_info() reads the model and user formula of a fit", {
  path <- test_path("assets/bmmfit_m3_ppcheck.rds")
  skip_if_not(file.exists(path), "M3 fixture not available (excluded by .Rbuildignore)")
  fit <- readRDS(path)

  expect_identical(
    parameter_info(fit),
    parameter_info(restructure(fit)$bmm$model, formula = fit$bmm$user_formula)
  )
})

# ===========================================================================
# parameters(), deprecated in 1.4.0
# ===========================================================================

test_that("parameters() warns once and returns what parameter_info() returns", {
  m <- sdm(resp_error = "y")
  warnings <- capture_warnings(out <- parameters(m))

  expect_length(warnings, 1)
  expect_match(warnings, "`parameters()` is deprecated as of bmm 1.4.0; use `parameter_info()`", fixed = TRUE)
  expect_identical(out, parameter_info(m))
})

test_that("parameters() passes the formula on to parameter_info()", {
  m <- m3(
    resp_cats = c("corr", "other", "npl"),
    num_options = c(1, 4, 5), version = "custom",
    links = list(c = "log", a = "log")
  )
  ff <- bmf(corr ~ b + a + c, other ~ b + a, npl ~ b, c ~ 1, a ~ 1)

  expect_warning(out <- parameters(m, formula = ff), "deprecated")
  expect_identical(out, parameter_info(m, formula = ff))
})

# ===========================================================================
# print.bmmodel()
# ===========================================================================

test_that("print.bmmodel() produces output", {
  m <- sdm(resp_error = "y")
  out <- capture.output(print(m))

  expect_true(any(grepl("Parameters:", out)))
  expect_true(any(grepl("parameter_info\\(\\)", out)))
})

test_that("print.bmmodel() shows fixed parameters", {
  m <- sdm(resp_error = "y")
  out <- capture.output(print(m))

  expect_true(any(grepl("Fixed:", out)))
  expect_true(any(grepl("mu", out)))
})

test_that("print.bmmodel() shows response variables and links", {
  m <- sdm(resp_error = "y")
  out <- capture.output(print(m))

  expect_true(any(grepl("Response:.*resp_error = y", out)))
  expect_true(any(grepl("Links:.*mu = tan_half; c = log; kappa = log", out)))
})

test_that("print.bmmodel() lists response variables one per line with annotations", {
  m <- ddm(rt = "myrt", response = "myresp")
  out <- capture.output(print(m))

  expect_true(any(grepl("Response:   rt = myrt (seconds)", out, fixed = TRUE)))
  expect_true(any(grepl(
    "            response = myresp (0/1 or logical; 1 = upper boundary)",
    out,
    fixed = TRUE
  )))
})

test_that("print.bmmodel() shows response bounds for circular models", {
  out <- capture.output(print(sdm(resp_error = "y")))

  expect_true(any(grepl("resp_error = y (radians in [-pi, pi])", out, fixed = TRUE)))
})

test_that("print.bmmodel() annotates aggregated ezdm response variables", {
  m <- ezdm(
    mean_rt = "mrt", var_rt = "vrt", n_upper = "nu", n_trials = "nt"
  )
  out <- capture.output(print(m))

  expect_true(any(grepl("mean_rt = mrt (seconds)", out, fixed = TRUE)))
  expect_true(any(grepl("var_rt = vrt (seconds^2)", out, fixed = TRUE)))
  expect_true(any(grepl(
    "n_upper = nu (count of upper-boundary responses)",
    out,
    fixed = TRUE
  )))
})

test_that("print.bmmodel() annotates aggregated sdt_yn response counts", {
  m <- sdt_yn(response = "n_old", stimulus = "stim", n_trials = "nt")
  out <- capture.output(print(m))

  expect_true(any(grepl(
    "response = n_old (count of 'old'/'signal' responses per cell)",
    out,
    fixed = TRUE
  )))
  expect_true(any(grepl("Fixed:.*sdratio = 0", out)))
})

test_that("print.bmmodel() lists vector-valued response variables", {
  m <- m3(
    resp_cats = c("corr", "other", "npl"),
    num_options = c(1, 4, 5),
    choice_rule = "simple",
    version = "ss"
  )
  out <- capture.output(print(m))

  expect_true(any(grepl("Response:.*resp_cats = corr, other, npl", out)))
  expect_true(any(grepl("Links:.*c = log; a = log", out)))
})
