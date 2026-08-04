# =============================================================================
# Tests for lnr model (model-specific tests)
# Distribution function tests are in test-distributions.R
# =============================================================================

# -----------------------------------------------------------------------------
# Model construction tests
# -----------------------------------------------------------------------------

test_that("lnr() creates simple model with correct structure", {
  model <- lnr(rt = "rt", response = "response", n_choices = 2)

  expect_s3_class(model, "bmmodel")
  expect_s3_class(model, "lnr")
  expect_s3_class(model, "lnr_simple")
  expect_equal(model$resp_vars$rt, "rt")
  expect_equal(model$resp_vars$response, "response")
  expect_equal(model$other_vars$n_choices, 2L)
  expect_equal(model$version, "simple")
})

test_that("lnr simple version has correct parameters", {
  model <- lnr(rt = "rt", response = "response", n_choices = 4)

  expect_true(all(c("correct", "error", "ndt", "s") %in%
                    names(model$parameters)))
  expect_equal(model$other_vars$n_choices, 4L)
})

test_that("lnr accepts custom links", {
  model <- lnr(rt = "rt", response = "response", n_choices = 2,
               links = list(correct = "log"))
  expect_equal(model$links$correct, "log")
  expect_equal(model$links$error, "identity")
})

test_that("lnr errors on invalid n_choices", {
  expect_error(lnr(rt = "rt", response = "response", n_choices = 1))
  expect_error(lnr(rt = "rt", response = "response", n_choices = 2.5))
})

test_that("lnr errors on missing required arguments", {
  expect_error(lnr(response = "response", n_choices = 2))
  expect_error(lnr(rt = "rt", n_choices = 2))
})

test_that("lnr() creates custom model with correct structure", {
  model <- lnr(rt = "rt", response = "resp", version = "custom")

  expect_s3_class(model, "lnr_custom")
  expect_equal(model$version, "custom")
  expect_null(model$other_vars$n_choices)
})

test_that("lnr custom accepts accumulators", {
  model <- lnr(rt = "rt", response = "resp", version = "custom",
               accumulators = c(target = 1, lure = 3))
  expect_equal(model$other_vars$accumulators, c(target = 1, lure = 3))
})

test_that("lnr deprecates n_alternatives and num_alternatives", {
  expect_warning(
    model <- lnr(rt = "rt", response = "response", n_alternatives = 2),
    "n_alternatives.*deprecated.*n_choices"
  )
  expect_equal(model$other_vars$n_choices, 2L)

  expect_warning(
    model <- lnr(rt = "rt", response = "resp", version = "custom",
                 num_alternatives = c(target = 1, lure = 3)),
    "num_alternatives.*deprecated.*accumulators"
  )
  expect_equal(model$other_vars$accumulators, c(target = 1, lure = 3))
})

test_that("lnr errors on version-specific alternative arguments", {
  expect_error(
    lnr(rt = "rt", response = "response", n_choices = 2,
        accumulators = c(correct = 1, error = 1)),
    "only supported for version 'custom'"
  )
  expect_error(
    lnr(rt = "rt", response = "resp", version = "custom", n_choices = 2),
    "only supported for version 'simple'"
  )
})

# -----------------------------------------------------------------------------
# Data validation tests (check_data — simple version)
# -----------------------------------------------------------------------------

test_that("check_data.lnr errors when required variables missing", {
  model <- lnr(rt = "rt", response = "response", n_choices = 2)

  expect_error(
    check_data(model, data.frame(x = 1), bmf(correct ~ 1)),
    "RT variable 'rt' is not present"
  )

  expect_error(
    check_data(model, data.frame(rt = 1), bmf(correct ~ 1)),
    "response variable 'response' is not present"
  )
})

test_that("check_data.lnr errors when RT contains NA", {
  model <- lnr(rt = "rt", response = "response", n_choices = 2)
  dat <- data.frame(rt = c(0.5, NA, 0.8), response = c(1, 1, 2))

  expect_error(
    check_data(model, dat, bmf(correct ~ 1)),
    "RT variable 'rt' contains.*NA"
  )
})

test_that("check_data.lnr errors when response contains NA", {
  model <- lnr(rt = "rt", response = "response", n_choices = 2)
  dat <- data.frame(rt = c(0.5, 0.6, 0.8), response = c(1, NA, 2))

  expect_error(
    check_data(model, dat, bmf(correct ~ 1)),
    "response variable 'response' contains.*NA"
  )
})

test_that("check_data.lnr errors when RT contains negative values", {
  model <- lnr(rt = "rt", response = "response", n_choices = 2)
  dat <- data.frame(rt = c(-0.5, 0.6, 0.8), response = c(1, 1, 2))

  expect_error(
    check_data(model, dat, bmf(correct ~ 1)),
    "reaction times are lower than zero"
  )
})

test_that("check_data.lnr warns when RT > 10 seconds", {
  model <- lnr(rt = "rt", response = "response", n_choices = 2)
  dat <- data.frame(
    rt = c(runif(99, 0.4, 1.5), 15),
    response = rep(c(1L, 2L), 50)
  )

  expect_warning(
    check_data(model, dat, bmf(correct ~ 1)),
    "larger than 10 seconds"
  )
})

test_that("check_data.lnr warns when RT < 0.1 seconds", {
  model <- lnr(rt = "rt", response = "response", n_choices = 2)
  dat <- data.frame(
    rt = c(0.05, runif(99, 0.4, 1.5)),
    response = rep(c(1L, 2L), 50)
  )

  expect_warning(
    check_data(model, dat, bmf(correct ~ 1)),
    "smaller than 0.100 seconds"
  )
})

test_that("check_data.lnr_simple errors on response out of range", {
  model <- lnr(rt = "rt", response = "response", n_choices = 2)

  dat <- data.frame(
    rt = runif(10, 0.4, 1.5), response = c(rep(1, 5), rep(3, 5))
  )
  expect_error(check_data(model, dat, bmf(correct ~ 1)), "integers in 1:2")

  dat <- data.frame(
    rt = runif(10, 0.4, 1.5), response = c(rep(0, 5), rep(1, 5))
  )
  expect_error(check_data(model, dat, bmf(correct ~ 1)), "integers in 1:2")
})

test_that("check_data.lnr_simple creates category mapping columns", {
  model <- lnr(rt = "rt", response = "response", n_choices = 4)
  dat <- data.frame(
    rt = runif(100, 0.4, 1.5),
    response = sample(1:4, 100, replace = TRUE)
  )

  result <- check_data(model, dat, bmf(correct ~ 1))
  expect_true(".lnr_cat" %in% names(result))
  expect_true(".lnr_n1" %in% names(result))
  expect_true(".lnr_n2" %in% names(result))
  expect_true(all(result$.lnr_n1 == 1L))
  expect_true(all(result$.lnr_n2 == 3L))
  expect_true(all(result$.lnr_cat[result$response == 1] == 1L))
  expect_true(all(result$.lnr_cat[result$response > 1] == 2L))
})

test_that("check_data.lnr_simple handles factor responses", {
  model <- lnr(rt = "rt", response = "response", n_choices = 3)
  dat <- data.frame(
    rt = runif(90, 0.4, 1.5),
    response = factor(rep(c("1", "2", "3"), 30))
  )

  result <- check_data(model, dat, bmf(correct ~ 1))
  expect_true(is.integer(result$response))
})

test_that("check_data.lnr returns a data.frame", {
  model <- lnr(rt = "rt", response = "response", n_choices = 2)
  dat <- data.frame(
    rt = runif(100, 0.4, 1.5), response = rep(c(1, 2), 50)
  )

  result <- check_data(model, dat, bmf(correct ~ 1))
  expect_s3_class(result, "data.frame")
})

# -----------------------------------------------------------------------------
# Data validation tests (check_data — custom version)
# -----------------------------------------------------------------------------

test_that("check_data.lnr_custom maps character responses to integers", {
  model <- lnr(rt = "rt", response = "resp", version = "custom")
  model$other_vars$resp_cats <- c("target", "lure")
  dat <- data.frame(
    rt = runif(100, 0.4, 1.5),
    resp = rep(c("target", "lure"), 50)
  )

  result <- check_data(model, dat, bmf(target ~ 1, lure ~ 1, ndt ~ 1))
  expect_true(all(result$.lnr_cat[result$resp == "target"] == 1L))
  expect_true(all(result$.lnr_cat[result$resp == "lure"] == 2L))
  expect_true(all(result$.lnr_n1 == 1L))
  expect_true(all(result$.lnr_n2 == 1L))
})

test_that("check_data.lnr_custom handles accumulators (integer)", {
  model <- lnr(rt = "rt", response = "resp", version = "custom",
               accumulators = c(correct = 1, other = 3, npl = 5))
  model$other_vars$resp_cats <- c("correct", "other", "npl")
  dat <- data.frame(
    rt = runif(90, 0.4, 1.5),
    resp = rep(c("correct", "other", "npl"), 30)
  )

  result <- check_data(model, dat, bmf(
    correct ~ 1, other ~ 1, npl ~ 1, ndt ~ 1
  ))
  expect_true(all(result$.lnr_n1 == 1L))
  expect_true(all(result$.lnr_n2 == 3L))
  expect_true(all(result$.lnr_n3 == 5L))
})

test_that("check_data.lnr_custom handles accumulators (column names)", {
  model <- lnr(rt = "rt", response = "resp", version = "custom",
               accumulators = c(target = "n_tgt", lure = "n_lure"))
  model$other_vars$resp_cats <- c("target", "lure")
  dat <- data.frame(
    rt = runif(100, 0.4, 1.5),
    resp = rep(c("target", "lure"), 50),
    n_tgt = 1L,
    n_lure = rep(c(3L, 5L), 50)
  )

  result <- check_data(model, dat, bmf(target ~ 1, lure ~ 1, ndt ~ 1))
  expect_true(all(result$.lnr_n1 == 1L))
  expect_equal(result$.lnr_n2, dat$n_lure)
})

test_that("check_data.lnr_custom requires exact accumulators names", {
  model <- lnr(rt = "rt", response = "resp", version = "custom",
               accumulators = c(target = 1, extra = 3))
  model$other_vars$resp_cats <- c("target", "lure")
  dat <- data.frame(
    rt = runif(100, 0.4, 1.5),
    resp = rep(c("target", "lure"), 50)
  )

  expect_error(
    check_data(model, dat, bmf(target ~ 1, lure ~ 1, ndt ~ 1)),
    "exactly the formula categories"
  )
})

test_that("check_data.lnr_custom errors on invalid accumulators counts", {
  model <- lnr(rt = "rt", response = "resp", version = "custom",
               accumulators = c(target = 1, lure = 0))
  model$other_vars$resp_cats <- c("target", "lure")
  dat <- data.frame(
    rt = runif(100, 0.4, 1.5),
    resp = rep(c("target", "lure"), 50)
  )

  expect_error(
    check_data(model, dat, bmf(target ~ 1, lure ~ 1, ndt ~ 1)),
    "positive integers"
  )

  model$other_vars$accumulators <- c(target = 1, lure = Inf)
  expect_error(
    check_data(model, dat, bmf(target ~ 1, lure ~ 1, ndt ~ 1)),
    "positive integers"
  )
})

test_that("check_data.lnr_custom validates accumulators columns", {
  model <- lnr(rt = "rt", response = "resp", version = "custom",
               accumulators = c(target = "n_tgt", lure = "n_lure"))
  model$other_vars$resp_cats <- c("target", "lure")
  dat <- data.frame(
    rt = runif(100, 0.4, 1.5),
    resp = rep(c("target", "lure"), 50),
    n_tgt = 1L,
    n_lure = c(rep(3L, 99), NA_integer_)
  )

  expect_error(
    check_data(model, dat, bmf(target ~ 1, lure ~ 1, ndt ~ 1)),
    "contains NA values"
  )

  dat$n_lure <- c(rep(3, 99), 2.5)
  expect_error(
    check_data(model, dat, bmf(target ~ 1, lure ~ 1, ndt ~ 1)),
    "integers >= 1"
  )
})

test_that("check_data.lnr_custom errors on unsupported accumulators types", {
  model <- lnr(rt = "rt", response = "resp", version = "custom",
               accumulators = list(target = 1, lure = 3))
  model$other_vars$resp_cats <- c("target", "lure")
  dat <- data.frame(
    rt = runif(100, 0.4, 1.5),
    resp = rep(c("target", "lure"), 50)
  )

  expect_error(
    check_data(model, dat, bmf(target ~ 1, lure ~ 1, ndt ~ 1)),
    "must be NULL, a named numeric vector"
  )
})

test_that("check_data.lnr_custom errors on mismatched response levels", {
  model <- lnr(rt = "rt", response = "resp", version = "custom")
  model$other_vars$resp_cats <- c("target", "lure")
  dat <- data.frame(
    rt = runif(100, 0.4, 1.5),
    resp = rep(c("target", "unknown"), 50)
  )

  expect_error(
    check_data(model, dat, bmf(target ~ 1, lure ~ 1, ndt ~ 1)),
    "not specified in the formula"
  )
})

test_that("check_data.lnr_custom errors on non-character responses", {
  model <- lnr(rt = "rt", response = "resp", version = "custom")
  model$other_vars$resp_cats <- c("target", "lure")
  dat <- data.frame(rt = runif(10, 0.4, 1.5), resp = rep(c(1, 2), 5))

  expect_error(
    check_data(model, dat, bmf(target ~ 1, lure ~ 1, ndt ~ 1)),
    "character labels"
  )
})

test_that("check_data.lnr_custom errors on reserved internal response levels", {
  model <- lnr(rt = "rt", response = "resp", version = "custom")
  model$other_vars$resp_cats <- c("target", "lure")
  dat <- data.frame(
    rt = runif(10, 0.4, 1.5),
    resp = rep(c("target", "ndt"), 5)
  )

  expect_error(
    check_data(model, dat, bmf(target ~ 1, lure ~ 1, ndt ~ 1)),
    "reserved internal parameter names"
  )
})

# -----------------------------------------------------------------------------
# check_model tests (custom version)
# -----------------------------------------------------------------------------

test_that("check_model.lnr_custom discovers category params from formula", {
  model <- lnr(rt = "rt", response = "resp", version = "custom")
  formula <- bmf(old ~ 1, new ~ 1, ndt ~ 1)

  updated <- check_model(model, formula = formula)
  expect_true("old" %in% names(updated$parameters))
  expect_true("new" %in% names(updated$parameters))
  expect_equal(updated$other_vars$resp_cats, c("old", "new"))
  expect_equal(updated$links$old, "identity")
  expect_equal(updated$links$new, "identity")
})

test_that("check_model.lnr_custom errors on Stan reserved words", {
  model <- lnr(rt = "rt", response = "resp", version = "custom")
  formula <- bmf(target ~ 1, ndt ~ 1)

  # "target" is a Stan reserved word
  expect_error(
    check_model(model, formula = formula),
    "Stan reserved words"
  )
})

test_that("check_model.lnr_custom errors on reserved internal parameter names", {
  model <- lnr(rt = "rt", response = "resp", version = "custom")
  formula <- bmf(mu ~ 1, ndt ~ 1)

  expect_error(
    check_model(model, formula = formula),
    "reserved internal parameter names"
  )
})

test_that("check_model.lnr_custom errors on category names ending in numbers", {
  model <- lnr(rt = "rt", response = "resp", version = "custom")
  formula <- bmf(correct ~ 1, error1 ~ 1, ndt ~ 1)

  expect_error(
    check_model(model, formula = formula),
    "cannot end in a number"
  )
})

test_that("check_model.lnr_custom errors on category names containing underscores", {
  model <- lnr(rt = "rt", response = "resp", version = "custom")
  formula <- bmf(correct ~ 1, error_a ~ 1, ndt ~ 1)

  expect_error(
    check_model(model, formula = formula),
    "cannot contain underscores"
  )
})

# -----------------------------------------------------------------------------
# Formula conversion tests
# -----------------------------------------------------------------------------

test_that("bmf2bf.lnr_simple creates correct brms formula", {
  model <- lnr(rt = "rt", response = "response", n_choices = 2)
  formula <- bmf(correct ~ 1, error ~ 1, ndt ~ 1)

  bf <- bmf2bf(model, formula)

  expect_s3_class(bf, "brmsformula")
  expect_true(grepl("vint", deparse(bf$formula)))
})

test_that("bmf2bf.lnr_custom creates correct brms formula", {
  model <- lnr(rt = "rt", response = "resp", version = "custom")
  model$other_vars$resp_cats <- c("fast", "slow", "guess")
  formula <- bmf(fast ~ 1, slow ~ 1, guess ~ 1, ndt ~ 1)

  bf <- bmf2bf(model, formula)
  bf_str <- deparse(bf$formula)

  expect_s3_class(bf, "brmsformula")
  expect_true(grepl("vint", bf_str))
  expect_true(grepl(".lnr_cat", bf_str))
  expect_true(grepl(".lnr_n1", bf_str))
  expect_true(grepl(".lnr_n3", bf_str))
})

# -----------------------------------------------------------------------------
# Stan code generation tests
# -----------------------------------------------------------------------------

test_that(".lnr_stan_code generates valid Stan for 2 categories", {
  code <- bmm:::.lnr_stan_code("lnr_simple", c("correct", "error"))

  expect_true(grepl("lnr_simple_lpdf", code))
  expect_true(grepl("vector correct", code))
  expect_true(grepl("vector error", code))
  expect_true(grepl("array\\[\\] int n1", code))
  expect_true(grepl("array\\[\\] int n2", code))
  # single-pass race combination shared with lba/rdm
  expect_true(grepl("log\\(n\\[win\\]\\)", code))
  expect_true(grepl("reps", code))
})

test_that(".lnr_stan_code generates valid Stan for 4 categories", {
  cats <- c("target", "similar", "other", "guess")
  code <- bmm:::.lnr_stan_code("lnr_custom", cats)

  expect_true(grepl("lnr_custom_lpdf", code))
  for (cat in cats) expect_true(grepl(paste0("vector ", cat), code))
  expect_true(grepl("array\\[\\] int n4", code))
  expect_true(grepl("array\\[4\\] real m", code))
  expect_true(grepl("array\\[4\\] int n", code))
})

test_that("LNR generated Stan code uses vectorized custom likelihoods", {
  dat <- rlnr(n = 20, m = c(-1, 0), s = c(1, 1), ndt = 0.2)
  simple_model <- lnr(rt = "rt", response = "response", n_choices = 2)
  simple_formula <- bmf(correct ~ 1, error ~ 1, ndt ~ 1)
  simple_code <- stancode(simple_formula, data = dat, model = simple_model,
                          backend = "cmdstanr")

  expect_true(grepl(
    "target \\+= lnr_simple_lpdf\\(Y \\| mu, correct, error, ndt, s, vint1, vint2, vint3\\);",
    simple_code
  ))
  expect_false(grepl(
    "for \\(n in 1:N\\) \\{\\s+target \\+= lnr_simple_lpdf",
    simple_code,
    perl = TRUE
  ))

  dat$response <- ifelse(dat$response == 1, "corr", "wrong")
  custom_model <- lnr(rt = "rt", response = "response", version = "custom")
  custom_formula <- bmf(corr ~ 1, wrong ~ 1, ndt ~ 1)
  custom_code <- stancode(custom_formula, data = dat, model = custom_model,
                          backend = "cmdstanr")

  expect_true(grepl("target \\+= lnr_custom_lpdf\\(Y \\|", custom_code))
  expect_false(grepl(
    "for \\(n in 1:N\\) \\{\\s+target \\+= lnr_custom_lpdf",
    custom_code,
    perl = TRUE
  ))
})

test_that("a threading request switches the vectorized family to sliced vars", {
  skip_on_cran()

  dat <- data.frame(
    rt = runif(100, 0.4, 1.5),
    response = sample(c(1L, 2L), 100, replace = TRUE)
  )

  # reduce_sum slices Y but passes vars through whole, so the vectorized
  # family must slice the vint columns itself when threading is requested
  model <- lnr(rt = "rt", response = "response", n_choices = 2)
  formula <- bmf(correct ~ 1, error ~ 1, ndt ~ 1)
  dat <- check_data(model, dat, formula)

  config <- configure_model(model, dat, formula)
  expect_equal(config$formula$family$vars, c("vint1", "vint2", "vint3"))

  attr(model, "threads") <- TRUE
  config <- configure_model(model, dat, formula)
  expect_equal(
    config$formula$family$vars,
    c("vint1[start:end]", "vint2[start:end]", "vint3[start:end]")
  )
})

test_that("the custom version slices every vint column when threading", {
  skip_on_cran()

  model <- lnr(rt = "rt", response = "resp", version = "custom",
               accumulators = c(fast = 1, mid = 2, slow = 3))
  formula <- bmf(fast ~ 1, mid ~ 1, slow ~ 1, ndt ~ 1)
  dat <- data.frame(
    rt = runif(99, 0.4, 1.5),
    resp = rep(c("fast", "mid", "slow"), 33)
  )
  model <- check_model(model, dat, formula)
  dat <- check_data(model, dat, formula)

  attr(model, "threads") <- TRUE
  config <- configure_model(model, dat, formula)
  expect_equal(
    config$formula$family$vars,
    paste0("vint", 1:4, "[start:end]")
  )
})

test_that("stancode emits thread-safe slicing when threads is passed", {
  skip_on_cran()

  dat <- rlnr(n = 100, m = c(-1, 0), s = c(1, 1), ndt = 0.2)
  model <- lnr(rt = "rt", response = "response", n_choices = 2)
  formula <- bmf(correct ~ 1, error ~ 1, ndt ~ 1)

  sc <- stancode(formula, dat, model = model, backend = "cmdstanr")
  expect_true(grepl(
    "lnr_simple_lpdf(Y | mu, correct, error, ndt, s, vint1, vint2, vint3);",
    sc, fixed = TRUE
  ))

  sc_threaded <- stancode(formula, dat, model = model, backend = "cmdstanr",
                          threads = brms::threading(2))
  expect_true(grepl(
    paste0("lnr_simple_lpdf(Y[start:end] | mu, correct, error, ndt, s, ",
           "vint1[start:end], vint2[start:end], vint3[start:end]);"),
    sc_threaded, fixed = TRUE
  ))
  expect_true(grepl("reduce_sum", sc_threaded))
})

# -----------------------------------------------------------------------------
# Model configuration tests
# -----------------------------------------------------------------------------

test_that("configure_model.lnr_simple returns correct components", {
  skip_on_cran()

  model <- lnr(rt = "rt", response = "response", n_choices = 2)
  dat <- data.frame(
    rt = runif(100, 0.4, 1.5),
    response = sample(c(1L, 2L), 100, replace = TRUE)
  )
  dat <- check_data(model, dat, bmf(correct ~ 1))
  formula <- bmf(correct ~ 1, error ~ 1, ndt ~ 1)

  config <- configure_model(model, dat, formula)

  expect_true(all(c("formula", "data", "stanvars") %in% names(config)))
  expect_s3_class(config$formula, "brmsformula")
  expect_s3_class(config$formula$family, "customfamily")
  expect_equal(config$formula$family$name, "lnr_simple")
  expect_true("correct" %in% config$formula$family$dpars)
  expect_true("error" %in% config$formula$family$dpars)
})

test_that("configure_model.lnr_custom returns correct components", {
  skip_on_cran()

  model <- lnr(rt = "rt", response = "resp", version = "custom",
               accumulators = c(fast = 1, slow = 3))
  dat <- data.frame(
    rt = runif(100, 0.4, 1.5),
    resp = rep(c("fast", "slow"), 50)
  )
  formula <- bmf(fast ~ 1, slow ~ 1, ndt ~ 1)

  model <- check_model(model, formula = formula)
  dat <- check_data(model, dat, formula)
  config <- configure_model(model, dat, formula)

  expect_true(all(c("formula", "data", "stanvars") %in% names(config)))
  expect_s3_class(config$formula, "brmsformula")
  expect_s3_class(config$formula$family, "customfamily")
  expect_equal(config$formula$family$name, "lnr_custom")
  expect_equal(config$formula$family$dpars, c("mu", "fast", "slow", "ndt", "s"))
})

test_that("log_lik_lnr_simple matches dlnr for a 2-choice observation", {
  prep <- structure(
    list(
      data = list(
        Y = c(0.8, 0.9),
        vint1 = c(1L, 2L),
        vint2 = c(1L, 1L),
        vint3 = c(1L, 1L)
      ),
      dpars = list(
        correct = matrix(c(-1, -1), nrow = 1),
        error = matrix(c(0, 0), nrow = 1),
        ndt = matrix(c(0.2, 0.2), nrow = 1),
        s = matrix(c(1, 1), nrow = 1)
      ),
      family = list(dpars = c("mu", "correct", "error", "ndt", "s"))
    ),
    class = "brmsprep"
  )

  expect_equal(
    log_lik_lnr_simple(1, prep),
    dlnr(0.8, 1, m = c(-1, 0), s = 1, ndt = 0.2, log = TRUE)
  )
  expect_equal(
    log_lik_lnr_simple(2, prep),
    dlnr(0.9, 2, m = c(-1, 0), s = 1, ndt = 0.2, log = TRUE)
  )
})

# -----------------------------------------------------------------------------
# Integration tests with mock backend
# -----------------------------------------------------------------------------

test_that("lnr simple version runs with mock backend (2-choice)", {
  skip_on_cran()

  dat <- rlnr(n = 100, m = c(-1, 0), s = c(1, 1), ndt = 0.2)
  model <- lnr(rt = "rt", response = "response", n_choices = 2)
  formula <- bmf(correct ~ 1, error ~ 1, ndt ~ 1)

  expect_silent(
    bmm(formula, dat, model, backend = "mock", mock = 1, rename = FALSE)
  )
})

test_that("lnr simple version runs with mock backend (4-choice)", {
  skip_on_cran()

  dat <- rlnr(n = 200, m = c(-1, 0, 0, 0), s = c(1, 1, 1, 1), ndt = 0.2)
  model <- lnr(rt = "rt", response = "response", n_choices = 4)
  formula <- bmf(correct ~ 1, error ~ 1, ndt ~ 1)

  expect_silent(
    bmm(formula, dat, model, backend = "mock", mock = 1, rename = FALSE)
  )
})

test_that("lnr simple with predictors runs with mock backend", {
  skip_on_cran()

  dat <- rlnr(n = 200, m = c(-1, 0), s = c(1, 1), ndt = 0.2)
  dat$condition <- rep(c("A", "B"), each = 100)
  model <- lnr(rt = "rt", response = "response", n_choices = 2)
  formula <- bmf(correct ~ condition, error ~ 1, ndt ~ 1)

  expect_silent(
    bmm(formula, dat, model, backend = "mock", mock = 1, rename = FALSE)
  )
})

test_that("lnr custom version runs with mock backend", {
  skip_on_cran()

  dat <- data.frame(
    rt = runif(150, 0.3, 2.0),
    resp = rep(c("fast", "medium", "slow"), 50)
  )
  model <- lnr(rt = "rt", response = "resp", version = "custom")
  formula <- bmf(fast ~ 1, medium ~ 1, slow ~ 1, ndt ~ 1)

  expect_silent(
    bmm(formula, dat, model, backend = "mock", mock = 1, rename = FALSE)
  )
})

test_that("lnr custom with accumulators runs with mock backend", {
  skip_on_cran()

  dat <- data.frame(
    rt = runif(100, 0.3, 2.0),
    resp = rep(c("correct", "other", "npl"), length.out = 100)
  )
  model <- lnr(rt = "rt", response = "resp", version = "custom",
               accumulators = c(correct = 1, other = 3, npl = 5))
  formula <- bmf(correct ~ 1, other ~ 1, npl ~ 1, ndt ~ 1)

  expect_silent(
    bmm(formula, dat, model, backend = "mock", mock = 1, rename = FALSE)
  )
})

test_that("lnr custom with predictors runs with mock backend", {
  skip_on_cran()

  dat <- data.frame(
    rt = runif(200, 0.3, 2.0),
    resp = rep(c("old", "new"), 100),
    condition = rep(c("easy", "hard"), each = 100)
  )
  model <- lnr(rt = "rt", response = "resp", version = "custom")
  formula <- bmf(old ~ condition, new ~ 1, ndt ~ 1)

  expect_silent(
    bmm(formula, dat, model, backend = "mock", mock = 1, rename = FALSE)
  )
})
