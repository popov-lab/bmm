# =============================================================================
# Tests for rdm model (model-specific tests)
# Distribution function tests are in test-distributions.R
# =============================================================================

# -----------------------------------------------------------------------------
# Model construction tests
# -----------------------------------------------------------------------------

test_that("rdm() creates simple model with correct structure", {
  model <- rdm(rt = "rt", response = "response", n_alternatives = 2)

  expect_s3_class(model, "bmmodel")
  expect_s3_class(model, "rdm")
  expect_s3_class(model, "rdm_simple")
  expect_equal(model$resp_vars$rt, "rt")
  expect_equal(model$resp_vars$response, "response")
  expect_equal(model$other_vars$n_alternatives, 2L)
  expect_equal(model$version, "simple")
})

test_that("rdm simple version has correct parameters", {
  model <- rdm(rt = "rt", response = "response", n_alternatives = 4)

  expect_true(all(c("driftc", "drifte", "bound", "ndt", "s", "A") %in%
                    names(model$parameters)))
  expect_equal(model$other_vars$n_alternatives, 4L)
})

test_that("rdm simple version has all log links", {
  model <- rdm(rt = "rt", response = "response", n_alternatives = 2)

  expect_equal(model$links$driftc, "log")
  expect_equal(model$links$drifte, "log")
  expect_equal(model$links$bound, "log")
  expect_equal(model$links$ndt, "log")
  expect_equal(model$links$s, "log")
  expect_equal(model$links$A, "log")
})

test_that("rdm has correct fixed parameters", {
  model <- rdm(rt = "rt", response = "response", n_alternatives = 2)

  expect_equal(model$fixed_parameters$mu, 0)
  expect_equal(model$fixed_parameters$s, 0)
  expect_equal(model$fixed_parameters$A, -100)
})

test_that("rdm accepts custom links", {
  model <- rdm(rt = "rt", response = "response", n_alternatives = 2,
               links = list(driftc = "identity"))
  expect_equal(model$links$driftc, "identity")
  expect_equal(model$links$drifte, "log")
})

test_that("rdm errors on invalid n_alternatives", {
  expect_error(rdm(rt = "rt", response = "response", n_alternatives = 1))
  expect_error(rdm(rt = "rt", response = "response", n_alternatives = 2.5))
})

test_that("rdm errors on missing required arguments", {
  expect_error(rdm(response = "response", n_alternatives = 2))
  expect_error(rdm(rt = "rt", n_alternatives = 2))
})

test_that("rdm() creates custom model with correct structure", {
  model <- rdm(rt = "rt", response = "resp", version = "custom")

  expect_s3_class(model, "rdm_custom")
  expect_equal(model$version, "custom")
  expect_null(model$other_vars$n_alternatives)
})

test_that("rdm custom accepts num_alternatives", {
  model <- rdm(rt = "rt", response = "resp", version = "custom",
               num_alternatives = c(corr = 1, err = 3))
  expect_equal(model$other_vars$num_alternatives, c(corr = 1, err = 3))
})

# -----------------------------------------------------------------------------
# Data validation tests (check_data — simple version)
# -----------------------------------------------------------------------------

test_that("check_data.rdm errors when required variables missing", {
  model <- rdm(rt = "rt", response = "response", n_alternatives = 2)

  expect_error(
    check_data(model, data.frame(x = 1), bmf(driftc ~ 1)),
    "RT variable 'rt' is not present"
  )

  expect_error(
    check_data(model, data.frame(rt = 1), bmf(driftc ~ 1)),
    "response variable 'response' is not present"
  )
})

test_that("check_data.rdm errors when RT contains NA", {
  model <- rdm(rt = "rt", response = "response", n_alternatives = 2)
  dat <- data.frame(rt = c(0.5, NA), response = c(1, 2))
  expect_error(check_data(model, dat, bmf(driftc ~ 1)), "NA values")
})

test_that("check_data.rdm errors when response contains NA", {
  model <- rdm(rt = "rt", response = "response", n_alternatives = 2)
  dat <- data.frame(rt = c(0.5, 0.6), response = c(1, NA))
  expect_error(check_data(model, dat, bmf(driftc ~ 1)), "NA values")
})

test_that("check_data.rdm errors when RT contains negative values", {
  model <- rdm(rt = "rt", response = "response", n_alternatives = 2)
  dat <- data.frame(rt = c(-0.1, 0.5), response = c(1, 2))
  expect_error(check_data(model, dat, bmf(driftc ~ 1)), "lower than zero")
})

test_that("check_data.rdm warns when RT > 10 seconds", {
  model <- rdm(rt = "rt", response = "response", n_alternatives = 2)
  dat <- data.frame(rt = c(0.5, 15), response = c(1, 2))
  expect_warning(check_data(model, dat, bmf(driftc ~ 1)), "larger than 10")
})

test_that("check_data.rdm warns when RT < 0.1 seconds", {
  model <- rdm(rt = "rt", response = "response", n_alternatives = 2)
  dat <- data.frame(rt = c(0.05, 0.5), response = c(1, 2))
  expect_warning(check_data(model, dat, bmf(driftc ~ 1)), "smaller than 0.100")
})

test_that("check_data.rdm_simple errors on response out of range", {
  model <- rdm(rt = "rt", response = "response", n_alternatives = 2)
  dat <- data.frame(rt = c(0.5, 0.6), response = c(1, 3))
  expect_error(check_data(model, dat, bmf(driftc ~ 1)), "integers in 1:2")
})

test_that("check_data.rdm_simple creates category mapping columns", {
  model <- rdm(rt = "rt", response = "response", n_alternatives = 3)
  dat <- data.frame(rt = c(0.5, 0.6, 0.7), response = c(1, 2, 3))
  result <- check_data(model, dat, bmf(driftc ~ 1))

  expect_equal(result$.rdm_cat, c(1L, 2L, 2L))
  expect_equal(result$.rdm_n1, c(1L, 1L, 1L))
  expect_equal(result$.rdm_n2, c(2L, 2L, 2L))
})

test_that("check_data.rdm_simple handles factor responses", {
  model <- rdm(rt = "rt", response = "response", n_alternatives = 2)
  dat <- data.frame(rt = c(0.5, 0.6), response = factor(c("1", "2")))
  result <- check_data(model, dat, bmf(driftc ~ 1))

  expect_equal(result$.rdm_cat, c(1L, 2L))
})

test_that("check_data.rdm returns a data.frame", {
  model <- rdm(rt = "rt", response = "response", n_alternatives = 2)
  dat <- data.frame(rt = c(0.5, 0.6), response = c(1, 2))
  result <- check_data(model, dat, bmf(driftc ~ 1))
  expect_s3_class(result, "data.frame")
})

# -----------------------------------------------------------------------------
# Data validation tests (check_data — custom version)
# -----------------------------------------------------------------------------

test_that("check_data.rdm_custom maps character responses to integers", {
  model <- rdm(rt = "rt", response = "resp", version = "custom")
  f <- bmf(corr ~ 1, err ~ 1, ndt ~ 1)
  model <- check_model(model, formula = f)
  dat <- data.frame(rt = c(0.5, 0.6), resp = c("corr", "err"))
  result <- check_data(model, dat, f)

  expect_equal(result$.rdm_cat, c(1L, 2L))
  expect_equal(result$.rdm_n1, c(1L, 1L))
  expect_equal(result$.rdm_n2, c(1L, 1L))
})

test_that("check_data.rdm_custom handles num_alternatives (integer)", {
  model <- rdm(rt = "rt", response = "resp", version = "custom",
               num_alternatives = c(corr = 1, err = 3))
  f <- bmf(corr ~ 1, err ~ 1, ndt ~ 1)
  model <- check_model(model, formula = f)
  dat <- data.frame(rt = c(0.5, 0.6), resp = c("corr", "err"))
  result <- check_data(model, dat, f)

  expect_equal(result$.rdm_n1, c(1L, 1L))
  expect_equal(result$.rdm_n2, c(3L, 3L))
})

test_that("check_data.rdm_custom handles num_alternatives (column names)", {
  model <- rdm(rt = "rt", response = "resp", version = "custom",
               num_alternatives = c(corr = "n_corr", err = "n_err"))
  f <- bmf(corr ~ 1, err ~ 1, ndt ~ 1)
  model <- check_model(model, formula = f)
  dat <- data.frame(rt = c(0.5, 0.6), resp = c("corr", "err"),
                    n_corr = c(1L, 1L), n_err = c(3L, 5L))
  result <- check_data(model, dat, f)

  expect_equal(result$.rdm_n1, c(1L, 1L))
  expect_equal(result$.rdm_n2, c(3L, 5L))
})

test_that("check_data.rdm_custom errors on mismatched response levels", {
  model <- rdm(rt = "rt", response = "resp", version = "custom")
  f <- bmf(corr ~ 1, err ~ 1, ndt ~ 1)
  model <- check_model(model, formula = f)
  dat <- data.frame(rt = c(0.5, 0.6), resp = c("corr", "unknown"))
  expect_error(check_data(model, dat, f), "not specified in the formula")
})

test_that("check_data.rdm_custom errors on non-character responses", {
  model <- rdm(rt = "rt", response = "resp", version = "custom")
  f <- bmf(corr ~ 1, err ~ 1, ndt ~ 1)
  model <- check_model(model, formula = f)
  dat <- data.frame(rt = c(0.5, 0.6), resp = c(1, 2))
  expect_error(check_data(model, dat, f), "character labels")
})

# -----------------------------------------------------------------------------
# check_model tests (custom version)
# -----------------------------------------------------------------------------

test_that("check_model.rdm_custom discovers category params from formula", {
  model <- rdm(rt = "rt", response = "resp", version = "custom")
  f <- bmf(corr ~ 1, err ~ 1, ndt ~ 1)
  model <- check_model(model, formula = f)

  expect_true("corr" %in% names(model$parameters))
  expect_true("err" %in% names(model$parameters))
  expect_equal(model$links$corr, "log")
  expect_equal(model$links$err, "log")
  expect_equal(model$other_vars$resp_cats, c("corr", "err"))
})

test_that("check_model.rdm_custom errors on Stan reserved words", {
  model <- rdm(rt = "rt", response = "resp", version = "custom")
  f <- bmf(corr ~ 1, real ~ 1, ndt ~ 1)
  expect_error(check_model(model, formula = f), "Stan reserved words")
})

# -----------------------------------------------------------------------------
# Formula conversion tests
# -----------------------------------------------------------------------------

test_that("bmf2bf.rdm_simple creates correct brms formula", {
  model <- rdm(rt = "rt", response = "response", n_alternatives = 2)
  bf <- bmf2bf(model, bmf(driftc ~ 1))
  resp_term <- deparse(bf$formula[[2]])
  expect_true(grepl("vint", resp_term))
  expect_true(grepl("\\.rdm_cat", resp_term))
  expect_true(grepl("\\.rdm_n1", resp_term))
  expect_true(grepl("\\.rdm_n2", resp_term))
})

test_that("bmf2bf.rdm_custom creates correct brms formula", {
  model <- rdm(rt = "rt", response = "resp", version = "custom")
  f <- bmf(corr ~ 1, err ~ 1, ndt ~ 1)
  model <- check_model(model, formula = f)
  bf <- bmf2bf(model, f)
  resp_term <- deparse(bf$formula[[2]])
  expect_true(grepl("vint", resp_term))
  expect_true(grepl("\\.rdm_cat", resp_term))
})

# -----------------------------------------------------------------------------
# Stan code generation tests
# -----------------------------------------------------------------------------

test_that(".rdm_stan_code generates valid Stan for A=0 (2 categories)", {
  code <- .rdm_stan_code("rdm_simple", c("driftc", "drifte"), has_A = FALSE)
  expect_true(grepl("rdm_simple_lpdf", code))
  expect_true(grepl("swald_lpdf", code))
  expect_true(grepl("swald_lccdf", code))
  expect_false(grepl("Phi\\(", code))
})

test_that(".rdm_stan_code generates valid Stan for A>0 (2 categories)", {
  code <- .rdm_stan_code("rdm_simple", c("driftc", "drifte"), has_A = TRUE)
  expect_true(grepl("rdm_simple_lpdf", code))
  expect_true(grepl("Phi\\(", code))
  expect_true(grepl("std_normal_lpdf", code))
  expect_true(grepl("b = A \\+ bound", code))
  expect_false(grepl("swald_lpdf", code))
})

test_that(".rdm_stan_code generates valid Stan for 4 categories", {
  cats <- c("cat1", "cat2", "cat3", "cat4")
  code <- .rdm_stan_code("rdm_custom", cats, has_A = FALSE)
  expect_true(grepl("rdm_custom_lpdf", code))
  expect_true(grepl("array\\[4\\]", code))
  expect_true(grepl("real cat1", code))
  expect_true(grepl("real cat4", code))
})

test_that(".rdm_stan_code generates A>0 Stan for custom version", {
  cats <- c("target", "lure", "npl")
  code <- .rdm_stan_code("rdm_custom", cats, has_A = TRUE)
  expect_true(grepl("rdm_custom_lpdf", code))
  expect_true(grepl("array\\[3\\]", code))
  expect_true(grepl("Phi\\(", code))
})

# -----------------------------------------------------------------------------
# Model configuration tests
# -----------------------------------------------------------------------------

test_that("configure_model.rdm_simple returns correct components", {
  model <- rdm(rt = "rt", response = "response", n_alternatives = 2)
  dat <- data.frame(rt = c(0.5, 0.6), response = c(1, 2))
  f <- bmf(driftc ~ 1, drifte ~ 1, bound ~ 1, ndt ~ 1)
  model <- check_model(model, data = dat, formula = f)
  dat <- check_data(model, dat, f)
  config <- configure_model(model, dat, f)

  expect_true("formula" %in% names(config))
  expect_true("data" %in% names(config))
  expect_true("stanvars" %in% names(config))
  expect_equal(config$formula$family$name, "rdm_simple")
  expect_true(all(c("mu", "driftc", "drifte", "bound", "ndt", "s", "A") %in%
                    config$formula$family$dpars))
})

test_that("configure_model.rdm_simple loads cswald_helper_functions", {
  model <- rdm(rt = "rt", response = "response", n_alternatives = 2)
  dat <- data.frame(rt = c(0.5, 0.6), response = c(1, 2))
  f <- bmf(driftc ~ 1, drifte ~ 1, bound ~ 1, ndt ~ 1)
  model <- check_model(model, data = dat, formula = f)
  dat <- check_data(model, dat, f)
  config <- configure_model(model, dat, f)

  stanvar_code <- paste(
    vapply(config$stanvars, function(x) x$scode, character(1)),
    collapse = "\n"
  )
  expect_true(grepl("swald_lpdf", stanvar_code))
})

# -----------------------------------------------------------------------------
# Integration tests with mock backend
# -----------------------------------------------------------------------------

test_that("rdm simple version runs with mock backend (2-choice)", {
  skip_on_cran()
  dat <- rrdm(n = 200, drift = c(3, 1.5), bound = 1, ndt = 0.2)
  model <- rdm(rt = "rt", response = "response", n_alternatives = 2)
  f <- bmf(driftc ~ 1, drifte ~ 1, bound ~ 1, ndt ~ 1)
  expect_no_error(
    bmm(f, dat, model, backend = "mock", mock_fit = 1, rename = FALSE)
  )
})

test_that("rdm simple version runs with mock backend (4-choice)", {
  skip_on_cran()
  dat <- rrdm(n = 200, drift = c(3, 1.5, 1.5, 1.5), bound = 1, ndt = 0.2)
  model <- rdm(rt = "rt", response = "response", n_alternatives = 4)
  f <- bmf(driftc ~ 1, drifte ~ 1, bound ~ 1, ndt ~ 1)
  expect_no_error(
    bmm(f, dat, model, backend = "mock", mock_fit = 1, rename = FALSE)
  )
})

test_that("rdm simple version runs with predictor", {
  skip_on_cran()
  dat <- rrdm(n = 200, drift = c(3, 1.5), bound = 1, ndt = 0.2)
  dat$cond <- rep(c("A", "B"), 100)
  model <- rdm(rt = "rt", response = "response", n_alternatives = 2)
  f <- bmf(driftc ~ 1 + cond, drifte ~ 1, bound ~ 1, ndt ~ 1)
  expect_no_error(
    bmm(f, dat, model, backend = "mock", mock_fit = 1, rename = FALSE)
  )
})

test_that("rdm custom version runs with mock backend", {
  skip_on_cran()
  dat <- rrdm(n = 200, drift = c(3, 1.5, 1), bound = 1, ndt = 0.2)
  cats <- c("corr", "err", "npl")
  dat$resp <- cats[dat$response]
  model <- rdm(rt = "rt", response = "resp", version = "custom",
               num_alternatives = c(corr = 1, err = 1, npl = 1))
  f <- bmf(corr ~ 1, err ~ 1, npl ~ 1, bound ~ 1, ndt ~ 1)
  expect_no_error(
    bmm(f, dat, model, backend = "mock", mock_fit = 1, rename = FALSE)
  )
})

test_that("rdm custom version runs with predictor", {
  skip_on_cran()
  dat <- rrdm(n = 200, drift = c(3, 1.5), bound = 1, ndt = 0.2)
  dat$resp <- ifelse(dat$response == 1, "corr", "err")
  dat$cond <- rep(c("A", "B"), 100)
  model <- rdm(rt = "rt", response = "resp", version = "custom")
  f <- bmf(corr ~ 1 + cond, err ~ 1, bound ~ 1, ndt ~ 1)
  expect_no_error(
    bmm(f, dat, model, backend = "mock", mock_fit = 1, rename = FALSE)
  )
})

test_that("rdm simple with A estimated runs with mock backend", {
  skip_on_cran()
  dat <- rrdm(n = 200, drift = c(3, 1.5), bound = 1, ndt = 0.2, A = 0.3)
  model <- rdm(rt = "rt", response = "response", n_alternatives = 2)
  f <- bmf(driftc ~ 1, drifte ~ 1, bound ~ 1, ndt ~ 1, A ~ 1)
  expect_no_error(
    bmm(f, dat, model, backend = "mock", mock_fit = 1, rename = FALSE)
  )
})
