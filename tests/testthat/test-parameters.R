library(testthat)

# ===========================================================================
# parameters() for bmmodel objects
# ===========================================================================

test_that("parameters() returns correct structure for sdm", {
  m <- sdm(resp_error = "y")
  p <- parameters(m)

  expect_s3_class(p, "bmm_parameters")
  expect_s3_class(p, "data.frame")
  expect_true(all(c("parameter", "description", "fixed", "value", "link") %in% names(p)))
  expect_equal(nrow(p), 3)
  expect_true(all(c("mu", "c", "kappa") %in% p$parameter))
})

test_that("parameters() flags fixed parameters for sdm", {
  m <- sdm(resp_error = "y")
  p <- parameters(m)

  mu_row <- p[p$parameter == "mu", ]
  expect_true(mu_row$fixed)
  expect_equal(mu_row$value, "0")

  c_row <- p[p$parameter == "c", ]
  expect_false(c_row$fixed)
  expect_true(is.na(c_row$value))
})

test_that("parameters() shows correct link functions", {
  m <- sdm(resp_error = "y")
  p <- parameters(m)

  expect_equal(p$link[p$parameter == "c"], "log")
  expect_equal(p$link[p$parameter == "kappa"], "log")
  expect_equal(p$link[p$parameter == "mu"], "tan_half")
})

test_that("parameters() works for imm model", {
  m <- imm(
    resp_error = "y", nt_features = "nt",
    nt_distances = "d", set_size = 2
  )
  p <- parameters(m)

  expect_s3_class(p, "bmm_parameters")
  expect_true(all(c("kappa", "a", "c", "s") %in% p$parameter))
})

test_that("parameters() works for mixture3p model", {
  m <- mixture3p(resp_error = "y", nt_features = "nt", set_size = 2)
  p <- parameters(m)

  expect_s3_class(p, "bmm_parameters")
  expect_true(all(c("thetat", "thetant", "kappa") %in% p$parameter))
})

test_that("parameters() works for m3 ss version", {
  m <- m3(
    resp_cats = c("corr", "other", "npl"),
    num_options = c(1, 4, 5), version = "ss"
  )
  p <- parameters(m)

  expect_s3_class(p, "bmm_parameters")
  expect_true(all(c("b", "c", "a") %in% p$parameter))
})

test_that("parameters() for m3 custom without formula shows note", {
  m <- m3(
    resp_cats = c("corr", "other", "npl"),
    num_options = c(1, 4, 5), version = "custom",
    links = list(c = "log", a = "log")
  )
  p <- parameters(m)

  expect_false(is.null(attr(p, "m3_note")))
  expect_match(attr(p, "m3_note"), "custom M3 model")
})

test_that("parameters() for m3 custom with formula discovers params", {
  m <- m3(
    resp_cats = c("corr", "other", "npl"),
    num_options = c(1, 4, 5), version = "custom",
    links = list(c = "log", a = "log")
  )
  ff <- bmf(corr ~ b + a + c, other ~ b + a, npl ~ b, c ~ 1, a ~ 1)
  p <- parameters(m, formula = ff)

  expect_true("a" %in% p$parameter)
  expect_true("c" %in% p$parameter)
  expect_null(attr(p, "m3_note"))
})

test_that("parameters() identifies free parameters for sdm", {
  m <- sdm(resp_error = "y")
  p <- parameters(m)

  free_pars <- p[!p$fixed, ]
  expect_equal(nrow(free_pars), 2)
  expect_true(all(c("c", "kappa") %in% free_pars$parameter))
})

test_that("parameters() includes descriptions", {
  m <- sdm(resp_error = "y")
  p <- parameters(m)

  expect_true(all(nchar(p$description) > 0))
  expect_true(is.character(p$description))
})

# ===========================================================================
# print.bmmodel()
# ===========================================================================

test_that("print.bmmodel() produces output", {
  m <- sdm(resp_error = "y")
  out <- capture.output(print(m))

  expect_true(any(grepl("Parameters:", out)))
  expect_true(any(grepl("parameters\\(\\)", out)))
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
