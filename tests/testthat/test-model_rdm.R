# =============================================================================
# Tests for rdm model (model-specific tests)
# The rdm distribution functions are tested at the bottom of this file
# =============================================================================

# -----------------------------------------------------------------------------
# Model construction tests
# -----------------------------------------------------------------------------

test_that("rdm() creates simple model with correct structure", {
  model <- rdm(rt = "rt", response = "response", n_choices = 2)

  expect_s3_class(model, "bmmodel")
  expect_s3_class(model, "rdm")
  expect_s3_class(model, "rdm_simple")
  expect_equal(model$resp_vars$rt, "rt")
  expect_equal(model$resp_vars$response, "response")
  expect_equal(model$other_vars$n_choices, 2L)
  expect_equal(model$version, "simple")
})

test_that("rdm simple version has correct parameters", {
  model <- rdm(rt = "rt", response = "response", n_choices = 4)

  expect_true(all(c("driftc", "drifte", "gap", "ndt", "s", "sp") %in%
                    names(model$parameters)))
  expect_equal(model$other_vars$n_choices, 4L)
})

test_that("rdm simple version has correct links", {
  model <- rdm(rt = "rt", response = "response", n_choices = 2)

  expect_equal(model$links$driftc, "log")
  expect_equal(model$links$drifte, "log")
  expect_equal(model$links$gap, "log")
  expect_equal(model$links$ndt, "log")
  expect_equal(model$links$s, "log")
  expect_equal(model$links$sp, "log")
})

test_that("rdm has correct fixed parameters", {
  model <- rdm(rt = "rt", response = "response", n_choices = 2)

  expect_equal(model$fixed_parameters$mu, 0)
  expect_equal(model$fixed_parameters$s, 0)
  expect_equal(model$fixed_parameters$sp, -100)
})

test_that("rdm refuses every link but log, at construction and afterwards", {
  expect_error(
    rdm(rt = "rt", response = "response", n_choices = 2,
        links = list(driftc = "identity")),
    "Unknown link function"
  )
  expect_error(
    rdm(rt = "rt", response = "resp", version = "custom",
        links = list(ndt = "softplus")),
    "Unknown link function"
  )

  # the documented idiom for changing a link after construction goes through
  # check_model(), which re-validates whatever differs from the constructor
  model <- rdm(rt = "rt", response = "response", n_choices = 2)
  model$links$ndt <- "identity"
  expect_error(
    check_model(model, formula = bmf(driftc ~ 1)),
    "Unknown link function"
  )
})

test_that("rdm refuses a link target it does not have", {
  expect_error(
    rdm(rt = "rt", response = "response", n_choices = 2,
        links = list(typo = "log")),
    "Unrecognized link target"
  )
  # the custom version's accumulators are named by the formula, so no link
  # target can be refused at construction
  model <- rdm(rt = "rt", response = "resp", version = "custom",
               links = list(target = "log"))
  expect_equal(model$links$target, "log")
})

test_that("report_priors() does not report the technical mu of the rdm family", {
  skip_on_cran()

  dat <- rrdm(n = 100, drift = c(3, 1.5), gap = 1, ndt = 0.2)
  fit <- bmm(
    bmf(driftc ~ 1, drifte ~ 1, gap ~ 1, ndt ~ 1), dat,
    rdm(rt = "rt", response = "response", n_choices = 2),
    backend = "mock", mock_fit = 1, rename = FALSE
  )

  out <- report_priors(fit)
  expect_false("mu" %in% out$parameter)
  expect_true(all(c("driftc", "drifte", "gap", "ndt") %in% out$parameter))
})

test_that("rdm errors on invalid n_choices", {
  expect_error(rdm(rt = "rt", response = "response", n_choices = 1))
  expect_error(rdm(rt = "rt", response = "response", n_choices = 2.5))
})

test_that("rdm errors on missing required arguments", {
  expect_error(rdm(response = "response", n_choices = 2))
  expect_error(rdm(rt = "rt", n_choices = 2))
})

test_that("rdm() creates custom model with correct structure", {
  model <- rdm(rt = "rt", response = "resp", version = "custom")

  expect_s3_class(model, "rdm_custom")
  expect_equal(model$version, "custom")
  expect_null(model$other_vars$n_choices)
})

test_that("rdm custom accepts accumulators", {
  model <- rdm(rt = "rt", response = "resp", version = "custom",
               accumulators = c(corr = 1, err = 3))
  expect_equal(model$other_vars$accumulators, c(corr = 1, err = 3))
})

test_that("rdm deprecates n_alternatives and num_alternatives", {
  expect_warning(
    model <- rdm(rt = "rt", response = "response", n_alternatives = 2),
    "n_alternatives.*deprecated.*n_choices"
  )
  expect_equal(model$other_vars$n_choices, 2L)

  expect_warning(
    model <- rdm(rt = "rt", response = "resp", version = "custom",
                 num_alternatives = c(corr = 1, err = 3)),
    "num_alternatives.*deprecated.*accumulators"
  )
  expect_equal(model$other_vars$accumulators, c(corr = 1, err = 3))
})

# -----------------------------------------------------------------------------
# Data validation tests (check_data — simple version)
# -----------------------------------------------------------------------------

test_that("check_data.rdm errors when required variables missing", {
  model <- rdm(rt = "rt", response = "response", n_choices = 2)

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
  model <- rdm(rt = "rt", response = "response", n_choices = 2)
  dat <- data.frame(rt = c(0.5, NA), response = c(1, 2))
  expect_error(check_data(model, dat, bmf(driftc ~ 1)), "NA values")
})

test_that("check_data.rdm errors when response contains NA", {
  model <- rdm(rt = "rt", response = "response", n_choices = 2)
  dat <- data.frame(rt = c(0.5, 0.6), response = c(1, NA))
  expect_error(check_data(model, dat, bmf(driftc ~ 1)), "NA values")
})

test_that("check_data.rdm errors when RT contains negative values", {
  model <- rdm(rt = "rt", response = "response", n_choices = 2)
  dat <- data.frame(rt = c(-0.1, 0.5), response = c(1, 2))
  expect_error(check_data(model, dat, bmf(driftc ~ 1)), "lower than zero")
})

test_that("check_data.rdm warns when RT > 10 seconds", {
  model <- rdm(rt = "rt", response = "response", n_choices = 2)
  dat <- data.frame(rt = c(0.5, 15), response = c(1, 2))
  expect_warning(check_data(model, dat, bmf(driftc ~ 1)), "larger than 10")
})

test_that("check_data.rdm warns when RT < 0.1 seconds", {
  model <- rdm(rt = "rt", response = "response", n_choices = 2)
  dat <- data.frame(rt = c(0.05, 0.5), response = c(1, 2))
  expect_warning(check_data(model, dat, bmf(driftc ~ 1)), "smaller than 0.100")
})

test_that("check_data.rdm_simple errors on response out of range", {
  model <- rdm(rt = "rt", response = "response", n_choices = 2)
  dat <- data.frame(rt = c(0.5, 0.6), response = c(1, 3))
  expect_error(check_data(model, dat, bmf(driftc ~ 1)), "integers in 1:2")
})

test_that("check_data.rdm_simple creates category mapping columns", {
  model <- rdm(rt = "rt", response = "response", n_choices = 3)
  dat <- data.frame(rt = c(0.5, 0.6, 0.7), response = c(1, 2, 3))
  result <- check_data(model, dat, bmf(driftc ~ 1))

  expect_equal(result$.rdm_cat, c(1L, 2L, 2L))
  expect_equal(result$.rdm_n1, c(1L, 1L, 1L))
  expect_equal(result$.rdm_n2, c(2L, 2L, 2L))
})

test_that("check_data.rdm_simple handles factor responses", {
  model <- rdm(rt = "rt", response = "response", n_choices = 2)
  dat <- data.frame(rt = c(0.5, 0.6), response = factor(c("1", "2")))
  result <- check_data(model, dat, bmf(driftc ~ 1))

  expect_equal(result$.rdm_cat, c(1L, 2L))
})

test_that("check_data.rdm_simple names the coding for a labelled response", {
  model <- rdm(rt = "rt", response = "response", n_choices = 2)
  dat <- data.frame(rt = c(0.5, 0.6), response = factor(c("yes", "no")))
  expect_error(check_data(model, dat, bmf(driftc ~ 1)),
               "1 = correct")
  expect_error(check_data(model, dat, bmf(driftc ~ 1)),
               "'yes', 'no'")
  expect_error(
    check_data(model, transform(dat, response = as.character(response)),
               bmf(driftc ~ 1)),
    "non-numeric label"
  )
})

test_that("check_data.rdm returns a data.frame", {
  model <- rdm(rt = "rt", response = "response", n_choices = 2)
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

test_that("check_data.rdm_custom handles accumulators (integer)", {
  model <- rdm(rt = "rt", response = "resp", version = "custom",
               accumulators = c(corr = 1, err = 3))
  f <- bmf(corr ~ 1, err ~ 1, ndt ~ 1)
  model <- check_model(model, formula = f)
  dat <- data.frame(rt = c(0.5, 0.6), resp = c("corr", "err"))
  result <- check_data(model, dat, f)

  expect_equal(result$.rdm_n1, c(1L, 1L))
  expect_equal(result$.rdm_n2, c(3L, 3L))
})

test_that("check_data.rdm_custom handles accumulators (column names)", {
  model <- rdm(rt = "rt", response = "resp", version = "custom",
               accumulators = c(corr = "n_corr", err = "n_err"))
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

# A count of zero, a negative count or a fraction reaches the Stan likelihood as
# log(n[response]) = -Inf / NaN, which ends the chain with no usable message
custom_counts_model <- function(accumulators) {
  model <- rdm(rt = "rt", response = "resp", version = "custom",
               accumulators = accumulators)
  check_model(model, formula = bmf(corr ~ 1, err ~ 1, ndt ~ 1))
}

test_that("check_data.rdm_custom refuses constant counts that are not positive integers", {
  f <- bmf(corr ~ 1, err ~ 1, ndt ~ 1)
  dat <- data.frame(rt = c(0.5, 0.6), resp = c("corr", "err"))

  expect_error(
    check_data(custom_counts_model(c(corr = 0, err = 2)), dat, f),
    "corr = 0"
  )
  expect_error(
    check_data(custom_counts_model(c(corr = -1, err = 2.5)), dat, f),
    "err = 2.5"
  )
  expect_error(
    check_data(custom_counts_model(c(corr = 1, err = Inf)), dat, f),
    "positive integer"
  )
  expect_error(
    check_data(custom_counts_model(list(corr = 1, err = 2)), dat, f),
    "must be NULL, a named numeric vector"
  )
  ok <- check_data(custom_counts_model(c(corr = 1, err = 2)), dat, f)
  expect_equal(ok$.rdm_n2, c(2L, 2L))
})

test_that("check_data.rdm_custom validates the values of an accumulators column", {
  f <- bmf(corr ~ 1, err ~ 1, ndt ~ 1)
  model <- custom_counts_model(c(corr = "nc", err = "ne"))
  dat <- data.frame(rt = c(0.5, 0.6, 0.7), resp = c("corr", "err", "corr"),
                    nc = c(1L, 1L, 1L), ne = c(2L, 2L, 2L))

  expect_error(check_data(model, transform(dat, ne = c(2, 2.5, 2)), f),
               "'ne' must contain integers >= 0")
  expect_error(check_data(model, transform(dat, nc = c(1, -1, 1)), f),
               "'nc' must contain integers >= 0")
  expect_error(check_data(model, transform(dat, nc = c(1, NA, 1)), f),
               "'nc' contains NA or non-finite")
  expect_error(check_data(model, transform(dat, nc = c("a", "b", "c")), f),
               "'nc' must be numeric")

  # a category may sit a trial out; only the winner needs an accumulator
  loser_out <- check_data(model, transform(dat, ne = c(0L, 2L, 0L)), f)
  expect_equal(loser_out$.rdm_n2, c(0L, 2L, 0L))
  expect_error(check_data(model, transform(dat, nc = c(1L, 1L, 0L)), f),
               "'corr' wins on trials where accumulators")
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

# .stan_reserved is the measured list: every candidate was handed to stanc as a
# category of a generated program. `N`, `Y` and `lprior` are names the
# surrounding brms/bmm code declares, and `Intercept` collides with the
# `real Intercept` brms writes for the response's own intercept ("Identifier
# "Intercept" is already in use", stanc 2.40, measured on this family). `log` is
# only a function name and stanc accepts a vector argument that shadows it.
test_that("check_model.rdm_custom refuses the names the program already uses", {
  model <- rdm(rt = "rt", response = "resp", version = "custom")
  refused <- function(name) {
    f <- bmf(corr ~ 1, gap ~ 1, ndt ~ 1)
    f[[name]] <- stats::as.formula(paste(name, "~ 1"))
    check_model(model, formula = f)
  }
  for (name in c("N", "Y", "lprior", "Intercept")) {
    expect_error(refused(name), "reserved", info = name)
  }
  expect_error(refused("mu"), "reserved internal parameter names")
  expect_silent(refused("log"))

  # Stan is case-sensitive and so is the comparison: stanc 2.40 accepts every
  # one of these as a category of the generated program (measured), and a
  # case-folded comparison against .stan_reserved would refuse them all
  for (name in c("Real", "Data", "Vector", "y", "Lprior", "Target")) {
    expect_silent(refused(name))
  }
  expect_error(refused("real"), "Stan reserved words")
})

# gap, ndt, s and sp are consumed as the model's own parameters before the
# category names are read, so a response level of that name can only ever be
# reported as missing from the formula
test_that("check_data.rdm_custom refuses a response level named after a parameter", {
  f <- bmf(corr ~ 1, err ~ 1, ndt ~ 1)
  model <- check_model(rdm(rt = "rt", response = "resp", version = "custom"),
                       formula = f)
  for (level in c("ndt", "s", "sp", "mu", "Intercept")) {
    dat <- data.frame(rt = c(0.5, 0.6), resp = c("corr", level))
    expect_error(check_data(model, dat, f),
                 "Response levels cannot use reserved", info = level)
  }
})

test_that("check_model.rdm_custom errors on category names ending in numbers", {
  model <- rdm(rt = "rt", response = "resp", version = "custom")
  f <- bmf(correct ~ 1, error1 ~ 1, gap ~ 1, ndt ~ 1)
  expect_error(check_model(model, data = NULL, formula = f), "cannot end in a number")
})

test_that("check_model.rdm_custom errors on category names containing underscores", {
  model <- rdm(rt = "rt", response = "resp", version = "custom")
  f <- bmf(correct ~ 1, error_a ~ 1, gap ~ 1, ndt ~ 1)
  expect_error(check_model(model, data = NULL, formula = f), "cannot contain underscores")
})

# -----------------------------------------------------------------------------
# check_formula tests
# -----------------------------------------------------------------------------

test_that("check_formula.rdm warns only where s has an intercept", {
  model <- rdm(rt = "rt", response = "response", n_choices = 2)
  dat <- data.frame(rt = rep(c(0.5, 0.6), 5), response = rep(c(1, 2), 5),
                    cond = rep(c("a", "b"), 5))
  warns <- function(f) {
    checked <- suppressMessages(check_model(model, data = dat, formula = f))
    suppressMessages(check_formula(checked, dat, f))
  }

  expect_warning(warns(bmf(driftc ~ 1, drifte ~ 1, gap ~ 1, ndt ~ 1, s ~ 1)),
                 "only their ratios are identified")
  expect_warning(warns(bmf(driftc ~ 1, drifte ~ 1, gap ~ 1, ndt ~ 1, s ~ cond)),
                 "scale")
  expect_no_warning(
    warns(bmf(driftc ~ 1, drifte ~ 1, gap ~ 1, ndt ~ 1, s ~ 0 + cond))
  )
  expect_no_warning(warns(bmf(driftc ~ 1, drifte ~ 1, gap ~ 1, ndt ~ 1)))
  expect_no_warning(
    warns(bmf(driftc ~ 1, drifte ~ 1, gap ~ 1, ndt ~ 1, s = log(0.9)))
  )
})

# -----------------------------------------------------------------------------
# Formula conversion tests
# -----------------------------------------------------------------------------

test_that("bmf2bf.rdm_simple creates correct brms formula", {
  model <- rdm(rt = "rt", response = "response", n_choices = 2)
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

# The family is unrolled over the categories (no per-row array temporaries,
# which cost 8-22% per gradient); the plain-Wald path is taken only when the
# generator is told the start point is off
test_that(".rdm_stan_code unrolls the plain-Wald race for sp off (2 categories)", {
  code <- .rdm_stan_code("rdm_simple", c("driftc", "drifte"), start_var = FALSE)
  expect_match(code, "real rdm_simple_lpdf(vector rt, vector mu, vector driftc, vector drifte, vector gap, vector ndt, vector s, vector sp, array[] int response, array[] int n1, array[] int n2)", fixed = TRUE)
  expect_match(code, "if (response[i] == 1) lp = log(n1[i]) + swald_lpdf(rt[i] | driftc[i], gap[i] + sp[i], ndt[i], s[i]);", fixed = TRUE)
  expect_match(code, "else lp = log(n2[i]) + swald_lpdf(rt[i] | drifte[i], gap[i] + sp[i], ndt[i], s[i]);", fixed = TRUE)
  expect_match(code, "if (reps2 > 0) lp += reps2 * swald_lccdf(rt[i] | drifte[i], gap[i] + sp[i], ndt[i], s[i]);", fixed = TRUE)
  expect_false(grepl("rdm_log_pdf", code, fixed = TRUE))
  expect_false(grepl("array[", code, fixed = TRUE) && grepl("drift_i", code, fixed = TRUE))
})

test_that(".rdm_stan_code unrolls the start-point race for sp on (2 categories)", {
  code <- .rdm_stan_code("rdm_simple", c("driftc", "drifte"), start_var = TRUE)
  expect_match(code, "int reps1 = (response[i] == 1) ? n1[i] - 1 : n1[i];", fixed = TRUE)
  expect_match(code, "if (response[i] == 1) lp = log(n1[i]) + rdm_log_pdf(t, driftc[i], gap[i], sp[i], s[i]);", fixed = TRUE)
  expect_match(code, "if (reps1 > 0) lp += reps1 * rdm_log_surv(t, driftc[i], gap[i], sp[i], s[i]);", fixed = TRUE)
  expect_false(grepl("swald_", code, fixed = TRUE))
})

test_that(".rdm_stan_code chains the winner branches for 4 categories", {
  cats <- c("cat1", "cat2", "cat3", "cat4")
  code <- .rdm_stan_code("rdm_custom", cats, start_var = FALSE)
  expect_match(code, "real rdm_custom_lpdf(", fixed = TRUE)
  expect_match(code, "array[] int n4)", fixed = TRUE)
  expect_match(code, "else if (response[i] == 3) lp = log(n3[i]) + swald_lpdf(rt[i] | cat3[i]", fixed = TRUE)
  expect_match(code, "else lp = log(n4[i]) + swald_lpdf(rt[i] | cat4[i]", fixed = TRUE)
  expect_equal(lengths(regmatches(code, gregexpr("if (reps", code, fixed = TRUE))), 4)
})

test_that(".rdm_stan_code generates the start-point race for the custom version", {
  cats <- c("target", "lure", "npl")
  code <- .rdm_stan_code("rdm_custom", cats, start_var = TRUE)
  expect_match(code, "rdm_log_pdf(t, lure[i], gap[i], sp[i], s[i])", fixed = TRUE)
  expect_match(code, "reps3 * rdm_log_surv(t, npl[i], gap[i], sp[i], s[i])", fixed = TRUE)
})

# -----------------------------------------------------------------------------
# Model configuration tests
# -----------------------------------------------------------------------------

# std_normal_lccdf(z) is -Inf from z = 8.26 with an infinite partial on Stan
# Math 5.3 (rstan / StanHeaders 2.39), which made a loser's survival NaN from
# t = 4 s at drift 5; the chunk spells every upper tail as std_normal_lcdf(-z)
test_that("the rdm Stan chunk never calls std_normal_lccdf", {
  chunk <- read_lines2(system.file("stan_chunks", "rdm_functions.stan", package = "bmm"))
  expect_false(grepl("std_normal_lccdf", chunk, fixed = TRUE))
  expect_true(grepl("std_normal_lcdf(", chunk, fixed = TRUE))
})

test_that("configure_model.rdm_simple returns correct components", {
  model <- rdm(rt = "rt", response = "response", n_choices = 2)
  dat <- data.frame(rt = c(0.5, 0.6), response = c(1, 2))
  f <- bmf(driftc ~ 1, drifte ~ 1, gap ~ 1, ndt ~ 1)
  model <- check_model(model, data = dat, formula = f)
  dat <- check_data(model, dat, f)
  config <- configure_model(model, dat, f)

  expect_true("formula" %in% names(config))
  expect_true("data" %in% names(config))
  expect_true("stanvars" %in% names(config))
  expect_equal(config$formula$family$name, "rdm_simple")
  expect_true(all(c("mu", "driftc", "drifte", "gap", "ndt", "s", "sp") %in%
                    config$formula$family$dpars))
  expect_false(config$formula$family$loop)
  expect_true(grepl("exp\\(ndt\\)", brms::stancode(
    config$formula, data = config$data, family = config$formula$family
  )))
})

test_that("configure_model.rdm_simple loads RDM helper functions", {
  model <- rdm(rt = "rt", response = "response", n_choices = 2)
  dat <- data.frame(rt = c(0.5, 0.6), response = c(1, 2))
  f <- bmf(driftc ~ 1, drifte ~ 1, gap ~ 1, ndt ~ 1)
  model <- check_model(model, data = dat, formula = f)
  dat <- check_data(model, dat, f)
  config <- configure_model(model, dat, f)

  stanvar_code <- paste(
    vapply(config$stanvars, function(x) x$scode, character(1)),
    collapse = "\n"
  )
  expect_true(grepl("swald_lpdf", stanvar_code))
  expect_true(grepl("rdm_log_surv", stanvar_code))
  expect_false(grepl("rdm_log_lik_one", stanvar_code))
})

test_that("configure_model.rdm_simple keeps ndt as a regular log-linked parameter", {
  model <- rdm(rt = "rt", response = "response", n_choices = 2)
  dat <- data.frame(rt = c(0.5, 0.6), response = c(1, 2), cond = c(0, 1))
  f <- bmf(driftc ~ cond, drifte ~ 1, gap ~ 1, ndt ~ cond)
  model <- check_model(model, data = dat, formula = f)
  dat <- check_data(model, dat, f)
  config <- configure_model(model, dat, f)
  stan_code <- brms::stancode(config$formula, data = config$data, family = config$formula$family)

  expect_false(grepl("ndtraw", stan_code))
  expect_false(grepl("inv_logit\\(nlp_ndtraw", stan_code))
  expect_true(grepl("X_ndt", stan_code))
  expect_true(grepl("X_driftc", stan_code))
})

test_that("create_initfun for rdm keeps intercept ndt draws in 10 to 50 ms", {
  model <- rdm(rt = "rt", response = "response", n_choices = 2)
  dat <- data.frame(rt = c(0.45, 0.62), response = c(1, 2))
  f <- bmf(driftc ~ 1, drifte ~ 1, gap ~ 1, ndt ~ 1)
  model <- check_model(model, data = dat, formula = f)
  dat <- check_data(model, dat, f)
  config <- configure_model(model, dat, f)

  init_fun <- create_initfun(model, dat, config$formula)
  vals <- replicate(100, exp(init_fun()[["Intercept_ndt"]]))

  expect_true(all(vals >= 0.01 & vals <= 0.05))
})

test_that("create_initfun for rdm keeps no-intercept ndt draws in 10 to 50 ms", {
  model <- rdm(rt = "rt", response = "response", n_choices = 2)
  dat <- data.frame(
    rt = c(0.45, 0.62, 0.58, 0.71),
    response = c(1, 2, 1, 2),
    cond = factor(c("A", "B", "A", "B"))
  )
  f <- bmf(driftc ~ 1, drifte ~ 1, gap ~ 1, ndt ~ 0 + cond)
  model <- check_model(model, data = dat, formula = f)
  dat <- check_data(model, dat, f)
  config <- configure_model(model, dat, f)

  init_fun <- create_initfun(model, dat, config$formula)
  vals <- exp(init_fun()[["b_ndt"]])

  expect_true(all(vals >= 0.01 & vals <= 0.05))
})

test_that("create_initfun for rdm ndt draws do not depend on observed RT range", {
  model <- rdm(rt = "rt", response = "response", n_choices = 2)
  dat_small_rt <- data.frame(rt = c(0.18, 0.23), response = c(1, 2))
  dat_large_rt <- data.frame(rt = c(0.75, 0.92), response = c(1, 2))
  f <- bmf(driftc ~ 1, drifte ~ 1, gap ~ 1, ndt ~ 1)

  model_small <- check_model(model, data = dat_small_rt, formula = f)
  data_small <- check_data(model_small, dat_small_rt, f)
  config_small <- configure_model(model_small, data_small, f)
  init_small <- create_initfun(model_small, data_small, config_small$formula)

  model_large <- check_model(model, data = dat_large_rt, formula = f)
  data_large <- check_data(model_large, dat_large_rt, f)
  config_large <- configure_model(model_large, data_large, f)
  init_large <- create_initfun(model_large, data_large, config_large$formula)

  ndt_small <- withr::with_seed(123, exp(init_small()[["Intercept_ndt"]]))
  ndt_large <- withr::with_seed(123, exp(init_large()[["Intercept_ndt"]]))

  expect_equal(ndt_small, ndt_large)
})

test_that("stancode for rdm includes user predictors for drift parameters", {
  dat <- data.frame(
    rt = c(0.5, 0.6),
    response = c(1, 2),
    cond = factor(c("A", "B"))
  )
  model <- rdm(rt = "rt", response = "response", n_choices = 2)
  f <- bmf(driftc ~ cond, drifte ~ 1, gap ~ 1, ndt ~ 1)
  code <- suppressWarnings(stancode(f, dat, model))

  expect_true(grepl("X_driftc", code))
  expect_true(grepl("Xc_driftc \\* b_driftc", code))
  expect_true(grepl("rdm_simple_lpdf\\(Y", code))
})

test_that("a threading request switches the vectorized family to sliced vars", {
  dat <- data.frame(rt = c(0.5, 0.6), response = c(1, 2))
  f <- bmf(driftc ~ 1, drifte ~ 1, gap ~ 1, ndt ~ 1)
  model <- rdm(rt = "rt", response = "response", n_choices = 2)
  model <- check_model(model, data = dat, formula = f)
  dat <- check_data(model, dat, f)

  # reduce_sum slices Y but passes vars through whole, so the vectorized
  # family must slice the vint columns itself when threading is requested
  config <- configure_model(model, dat, f)
  expect_equal(config$formula$family$vars, c("vint1", "vint2", "vint3"))

  withr::with_options(list(brms.threads = brms::threading(2)), {
    expect_equal(
      configure_model(model, dat, f)$formula$family$vars,
      c("vint1[start:end]", "vint2[start:end]", "vint3[start:end]")
    )
  })

  # force = TRUE tells brms to compile with threading but leave the generated
  # code alone, so start/end are never defined and slicing would not compile
  withr::with_options(list(brms.threads = brms::threading(2, force = TRUE)), {
    expect_equal(
      configure_model(model, dat, f)$formula$family$vars,
      c("vint1", "vint2", "vint3")
    )
  })

  # brms accepts a bare number for the option
  withr::with_options(list(brms.threads = 2), {
    expect_equal(
      configure_model(model, dat, f)$formula$family$vars,
      c("vint1[start:end]", "vint2[start:end]", "vint3[start:end]")
    )
  })
})

test_that("the custom version slices every vint column when threading", {
  dat <- data.frame(
    rt = runif(99, 0.4, 1.5),
    resp = rep(c("fast", "mid", "slow"), 33)
  )
  f <- bmf(fast ~ 1, mid ~ 1, slow ~ 1, gap ~ 1, ndt ~ 1)
  model <- rdm(rt = "rt", response = "resp", version = "custom",
               accumulators = c(fast = 1, mid = 2, slow = 3))
  model <- check_model(model, data = dat, formula = f)
  dat <- check_data(model, dat, f)

  expect_equal(
    configure_model(model, dat, f)$formula$family$vars,
    paste0("vint", 1:4)
  )
  withr::with_options(list(brms.threads = brms::threading(2)), {
    expect_equal(
      configure_model(model, dat, f)$formula$family$vars,
      paste0("vint", 1:4, "[start:end]")
    )
  })
})

test_that("an explicit threads = NULL beats a global threading option", {
  skip_on_cran()

  dat <- rrdm(n = 100, drift = c(3, 1.5), gap = 1, ndt = 0.2)
  model <- rdm(rt = "rt", response = "response", n_choices = 2)
  f <- bmf(driftc ~ 1, drifte ~ 1, gap ~ 1, ndt ~ 1)

  # brms reads threads = NULL as "threading off" and generates serial code, so
  # the family must not slice its vint columns even though the option is set
  withr::local_options(brms.threads = brms::threading(2))
  code <- suppressWarnings(stancode(f, dat, model, backend = "cmdstanr",
                                    threads = NULL))
  expect_true(grepl(
    "rdm_simple_lpdf(Y | mu, driftc, drifte, gap, ndt, s, sp, vint1, vint2, vint3);",
    code, fixed = TRUE
  ))
  expect_false(grepl("start:end", code, fixed = TRUE))
})

test_that("stancode emits thread-safe slicing when threads is passed", {
  skip_on_cran()

  dat <- rrdm(n = 100, drift = c(3, 1.5), gap = 1, ndt = 0.2)
  model <- rdm(rt = "rt", response = "response", n_choices = 2)
  f <- bmf(driftc ~ 1, drifte ~ 1, gap ~ 1, ndt ~ 1)

  code <- suppressWarnings(stancode(f, dat, model, backend = "cmdstanr"))
  expect_true(grepl(
    "rdm_simple_lpdf(Y | mu, driftc, drifte, gap, ndt, s, sp, vint1, vint2, vint3);",
    code, fixed = TRUE
  ))

  threaded <- suppressWarnings(stancode(f, dat, model, backend = "cmdstanr",
                                        threads = brms::threading(2)))
  expect_true(grepl(
    paste0("rdm_simple_lpdf(Y[start:end] | mu, driftc, drifte, gap, ndt, ",
           "s, sp, vint1[start:end], vint2[start:end], vint3[start:end]);"),
    threaded, fixed = TRUE
  ))
  expect_true(grepl("reduce_sum", threaded))
})

# -----------------------------------------------------------------------------
# Integration tests with mock backend
# -----------------------------------------------------------------------------

test_that("rdm simple version runs with mock backend (2-choice)", {
  skip_on_cran()
  dat <- rrdm(n = 200, drift = c(3, 1.5), gap = 1, ndt = 0.2)
  model <- rdm(rt = "rt", response = "response", n_choices = 2)
  f <- bmf(driftc ~ 1, drifte ~ 1, gap ~ 1, ndt ~ 1)
  expect_no_error(
    bmm(f, dat, model, backend = "mock", mock_fit = 1, rename = FALSE)
  )
})

test_that("rdm simple version runs with mock backend (4-choice)", {
  skip_on_cran()
  dat <- rrdm(n = 200, drift = c(3, 1.5, 1.5, 1.5), gap = 1, ndt = 0.2)
  model <- rdm(rt = "rt", response = "response", n_choices = 4)
  f <- bmf(driftc ~ 1, drifte ~ 1, gap ~ 1, ndt ~ 1)
  expect_no_error(
    bmm(f, dat, model, backend = "mock", mock_fit = 1, rename = FALSE)
  )
})

test_that("rdm simple version runs with predictor", {
  skip_on_cran()
  dat <- rrdm(n = 200, drift = c(3, 1.5), gap = 1, ndt = 0.2)
  dat$cond <- rep(c("A", "B"), 100)
  model <- rdm(rt = "rt", response = "response", n_choices = 2)
  f <- bmf(driftc ~ 1 + cond, drifte ~ 1, gap ~ 1, ndt ~ 1)
  expect_no_error(
    bmm(f, dat, model, backend = "mock", mock_fit = 1, rename = FALSE)
  )
})

test_that("rdm custom version runs with mock backend", {
  skip_on_cran()
  dat <- rrdm(n = 200, drift = c(3, 1.5, 1), gap = 1, ndt = 0.2)
  cats <- c("corr", "err", "npl")
  dat$resp <- cats[dat$response]
  model <- rdm(rt = "rt", response = "resp", version = "custom",
               accumulators = c(corr = 1, err = 1, npl = 1))
  f <- bmf(corr ~ 1, err ~ 1, npl ~ 1, gap ~ 1, ndt ~ 1)
  expect_no_error(
    bmm(f, dat, model, backend = "mock", mock_fit = 1, rename = FALSE)
  )
})

test_that("rdm custom version runs with predictor", {
  skip_on_cran()
  dat <- rrdm(n = 200, drift = c(3, 1.5), gap = 1, ndt = 0.2)
  dat$resp <- ifelse(dat$response == 1, "corr", "err")
  dat$cond <- rep(c("A", "B"), 100)
  model <- rdm(rt = "rt", response = "resp", version = "custom")
  f <- bmf(corr ~ 1 + cond, err ~ 1, gap ~ 1, ndt ~ 1)
  expect_no_error(
    bmm(f, dat, model, backend = "mock", mock_fit = 1, rename = FALSE)
  )
})

# A constant sp other than the default is a start-point range the user asked
# for, not a request for the plain-Wald path; the R posterior methods read
# A = sp either way, so they must agree with the density at that range
test_that("a fixed sp is honoured as the start-point range it names", {
  dat <- data.frame(rt = c(0.5, 0.6), response = c(1, 2))
  model <- rdm(rt = "rt", response = "response", n_choices = 2)

  configure <- function(f) {
    checked <- check_model(model, data = dat, formula = f)
    configure_model(checked, check_data(checked, dat, f), f)
  }
  default <- configure(bmf(driftc ~ 1, drifte ~ 1, gap ~ 1, ndt ~ 1))
  fixed <- configure(bmf(driftc ~ 1, drifte ~ 1, gap ~ 1, ndt ~ 1, sp = log(0.3)))
  free <- configure(bmf(driftc ~ 1, drifte ~ 1, gap ~ 1, ndt ~ 1, sp ~ 1))
  code_of <- function(config) config$stanvars[[3]]$scode
  expect_match(code_of(default), "swald_lpdf(", fixed = TRUE)
  expect_false(grepl("rdm_log_pdf(", code_of(default), fixed = TRUE))
  expect_match(code_of(fixed), "rdm_log_pdf(", fixed = TRUE)
  expect_match(code_of(free), "rdm_log_pdf(", fixed = TRUE)

  prep <- structure(
    list(
      ndraws = 2L,
      data = list(Y = c(0.5, 0.6), vint1 = c(1L, 2L), vint2 = c(1L, 1L), vint3 = c(1L, 1L)),
      dpars = list(driftc = c(3, 3.2), drifte = c(1.5, 1.4), gap = c(0.7, 0.8),
                   ndt = c(0.2, 0.22), s = c(1, 1), sp = c(0.3, 0.3)),
      family = list(dpars = c("mu", "driftc", "drifte", "gap", "ndt", "s", "sp"))
    ),
    class = "brmsprep"
  )
  for (i in 1:2) {
    per_draw <- vapply(1:2, function(k) {
      drdm(prep$data$Y[i], prep$data$vint1[i],
           drift = c(prep$dpars$driftc[k], prep$dpars$drifte[k]),
           gap = prep$dpars$gap[k], sp = 0.3, ndt = prep$dpars$ndt[k], s = 1,
           log = TRUE)
    }, numeric(1))
    expect_equal(
      .rdm_log_lik(i, prep, cat_names = c("driftc", "drifte"), n_cats = 2),
      per_draw, tolerance = 1e-12
    )
  }
})

test_that("rdm simple with sp estimated runs with mock backend", {
  skip_on_cran()
  dat <- rrdm(n = 200, drift = c(3, 1.5), gap = 0.7, sp = 0.3, ndt = 0.2)
  model <- rdm(rt = "rt", response = "response", n_choices = 2)
  f <- bmf(driftc ~ 1, drifte ~ 1, gap ~ 1, ndt ~ 1, sp ~ 1)
  expect_no_error(
    bmm(f, dat, model, backend = "mock", mock_fit = 1, rename = FALSE)
  )
})

# -----------------------------------------------------------------------------
# Post-processing tests
# -----------------------------------------------------------------------------

# three draws that differ in every parameter, and three observations that differ
# in how many error accumulators race
rdm_epred_prep <- function() {
  structure(
    list(
      ndraws = 3L, nobs = 3L,
      data = list(Y = c(0.6, 0.7, 0.8), vint1 = c(1L, 2L, 1L),
                  vint2 = c(1L, 1L, 1L), vint3 = c(1L, 2L, 3L)),
      dpars = list(
        driftc = c(3, 2.5, 4), drifte = c(1.5, 1.2, 2),
        gap = c(0.8, 1.2, 0.6), ndt = c(0.2, 0.15, 0.25),
        s = c(1, 0.9, 1.1), sp = c(1e-10, 0.3, 0.2)
      ),
      family = list(dpars = c("mu", "driftc", "drifte", "gap", "ndt", "s", "sp"))
    ),
    class = "brmsprep"
  )
}

test_that("posterior_epred for rdm is deterministic and integrates the race", {
  withr::local_seed(7)
  prep <- rdm_epred_prep()
  epred <- posterior_epred_rdm_simple(prep)
  expect_identical(epred, posterior_epred_rdm_simple(prep))

  n_mc <- 1e5
  for (i in 1:3) {
    for (k in 1:3) {
      race <- .rdm_race(
        drift = matrix(c(prep$dpars$driftc[k], prep$dpars$drifte[k]), n_mc, 2,
                       byrow = TRUE),
        gap = rep(prep$dpars$gap[k], n_mc), A = rep(prep$dpars$sp[k], n_mc),
        s = rep(prep$dpars$s[k], n_mc),
        counts = matrix(c(prep$data$vint2[i], prep$data$vint3[i]), n_mc, 2,
                        byrow = TRUE)
      )
      mc <- mean(race$rt) + prep$dpars$ndt[k]
      se <- stats::sd(race$rt) / sqrt(n_mc)
      expect_lt(abs(epred[k, i] - mc), 4 * se)
    }
  }
})

# The fake prep above carries every parameter as a draw-length vector; on a real
# fit brms stores the fixed s and sp as scalars, which only this path grows to
# ndraws. posterior_epred() picks a random subset of draws unless draw_ids is
# given, so the determinism is pinned at fixed draw_ids.
test_that("posterior_epred for rdm recovers the mean RT of a fitted model", {
  path <- test_path("assets", "bmmfit_rdm_ppcheck.rds")
  skip_if_not(file.exists(path), "fixture not available")
  fit <- readRDS(path)

  epred <- brms::posterior_epred(fit, draw_ids = 1:20)
  expect_identical(epred, brms::posterior_epred(fit, draw_ids = 1:20))
  expect_identical(dim(epred), c(20L, nrow(fit$data)))
  expect_lt(abs(mean(epred) - mean(fit$data$rt)), 0.02)
})

# The winner of a K = 3 simple race is "some error accumulator", which is
# log(K - 1) above one named accumulator winning. Constant in the parameters, so
# posteriors are unaffected, but loo() and waic() are shifted by it; pinned here
# so that a change of convention is deliberate (see ?rdm)
test_that("log_lik for the simple version carries the category constant", {
  prep <- structure(
    list(
      ndraws = 1L,
      data = list(Y = c(0.6, 0.7), vint1 = c(1L, 2L),
                  vint2 = c(1L, 1L), vint3 = c(2L, 2L)),
      dpars = list(driftc = 3, drifte = 1.5, gap = 0.9, ndt = 0.2, s = 1,
                   sp = 0.25),
      family = list(dpars = c("mu", "driftc", "drifte", "gap", "ndt", "s", "sp"))
    ),
    class = "brmsprep"
  )
  per_accumulator <- function(i, response) {
    drdm(prep$data$Y[i], response,
         drift = c(prep$dpars$driftc, rep(prep$dpars$drifte, 2)),
         gap = prep$dpars$gap, ndt = prep$dpars$ndt, s = prep$dpars$s,
         sp = prep$dpars$sp, log = TRUE)
  }
  ll <- vapply(1:2, .rdm_log_lik, numeric(1), prep = prep,
               cat_names = c("driftc", "drifte"), n_cats = 2)
  expect_equal(ll[1] - per_accumulator(1, 1), 0, tolerance = 1e-12)
  expect_equal(ll[2] - per_accumulator(2, 2), log(2), tolerance = 1e-12)
  expect_equal(per_accumulator(2, 2), per_accumulator(2, 3), tolerance = 1e-12)
})

# .rdm_race() is what pp_simulate.rdm() and posterior_predict() both draw from,
# so its joint law over (response, rt) is the one the pp_check panels show. A
# simulator that drew the response from the marginal choice probabilities would
# match both margins and fail here.
test_that(".rdm_race() reproduces the joint law of drdm()", {
  withr::local_seed(11)
  drift <- c(3, 1.4)
  counts <- c(1L, 2L)
  gap <- 1
  sp <- 0.25
  ndt <- 0.2
  n <- 2e5
  race <- .rdm_race(
    drift = matrix(drift, n, 2, byrow = TRUE), gap = rep(gap, n),
    A = rep(sp, n), s = rep(1, n),
    counts = matrix(counts, n, 2, byrow = TRUE)
  )
  rt <- race$rt + ndt
  cut <- stats::median(rt)

  for (j in 1:2) {
    density_j <- function(t) {
      exp(.rdm_race_lpdf(
        t = t - ndt, response = rep(j, length(t)),
        drift = matrix(drift, length(t), 2, byrow = TRUE),
        counts = matrix(counts, length(t), 2, byrow = TRUE),
        gap = rep(gap, length(t)), A = rep(sp, length(t)), s = rep(1, length(t))
      ))
    }
    p_sim <- mean(race$response == j & rt < cut)
    p_exact <- stats::integrate(density_j, ndt, cut, rel.tol = 1e-10)$value
    expect_lt(abs(p_sim - p_exact) / sqrt(p_sim * (1 - p_sim) / n), 4)
  }
})

# -----------------------------------------------------------------------------
# Distribution function tests
# -----------------------------------------------------------------------------

test_that("drdm returns positive densities for valid inputs", {
  d <- drdm(c(0.5, 0.6), c(1, 2), drift = c(3, 1.5),
            gap = 1, ndt = 0.2)
  expect_true(all(d > 0))
})

test_that("drdm remains finite close to ndt when sp > 0", {
  ll <- drdm(
    rt = c(0.201, 0.202, 0.205, 0.210),
    response = rep(1, 4),
    drift = c(3, 1.5, 1),
    gap = 1,
    ndt = 0.2,
    sp = 0.05,
    log = TRUE
  )
  expect_true(all(is.finite(ll)))
})

test_that("drdm returns -Inf when rt is at or below ndt", {
  ll <- drdm(
    rt = c(0.2, 0.199),
    response = c(1, 1),
    drift = c(3, 1.5),
    gap = 1,
    ndt = 0.2,
    log = TRUE
  )

  expect_equal(ll, c(-Inf, -Inf))
})

test_that("drdm with positive sp approaches the sp=0 likelihood as sp -> 0", {
  ll_sp0 <- drdm(
    rt = 0.6,
    response = 1,
    drift = c(3, 1.5, 1),
    gap = 1,
    ndt = 0.2,
    sp = 0,
    log = TRUE
  )
  ll_sp_eps <- drdm(
    rt = 0.6,
    response = 1,
    drift = c(3, 1.5, 1),
    gap = 1,
    ndt = 0.2,
    sp = 1e-8,
    log = TRUE
  )
  expect_equal(ll_sp_eps, ll_sp0, tolerance = 1e-5)
})

# The likelihood of a Wald race depends on (drift, gap, sp) only relative to the
# diffusion scale s, so scaling all four by one constant leaves it unchanged.
# The identity holds by construction at sp = 0 and is broken at sp > 0 by any
# term of the start-point CDF whose coefficient carries the wrong power of s.
test_that("the race likelihood is invariant to a common scale of drift, gap, sp and s", {
  rt <- c(0.35, 0.5, 0.8, 1.4)
  response <- c(1L, 2L, 1L, 2L)
  eval_one <- function(i, scale, sp) {
    drdm(rt[i], response[i], drift = scale * c(3, 1.5), gap = scale * 0.8,
         sp = scale * sp, ndt = 0.2, s = scale, log = TRUE)
  }
  for (sp in c(0, 0.25)) {
    base <- vapply(seq_along(rt), eval_one, numeric(1), scale = 1, sp = sp)
    for (scale in c(0.5, 2, 7)) {
      scaled <- vapply(seq_along(rt), eval_one, numeric(1), scale = scale, sp = sp)
      expect_equal(scaled, base, tolerance = 1e-10)
    }
  }
})

test_that("drdm() integrates to one over responses and time when s != 1", {
  density <- function(rt, response, s) {
    vapply(rt, function(x) {
      drdm(x, response, drift = c(3, 1.5), gap = 0.7, sp = 0.3, ndt = 0.2, s = s)
    }, numeric(1))
  }
  for (s in c(0.5, 0.8, 1.5)) {
    total <- sum(vapply(1:2, function(response) {
      stats::integrate(density, 0.2, Inf, response = response, s = s,
                       rel.tol = 1e-9)$value
    }, numeric(1)))
    expect_equal(total, 1, tolerance = 1e-6)
  }
})

# log_lik() hands the kernels one observation with a vector of draws, so the
# parameters vary along the vector and different draws need different branches
# of the same kernel (start-point range below and above the midpoint switch,
# decision time near and far from the non-decision time, either sign of the
# exponential difference in the survival). Every branch must subset every
# argument by the same mask; a stray full-length vector recycles silently.
test_that("a mixed vector of draws takes each branch it needs", {
  n_draws <- 6L
  prep <- structure(
    list(
      ndraws = n_draws,
      data = list(Y = 0.55, vint1 = 1L, vint2 = 1L, vint3 = 2L),
      dpars = list(
        driftc = c(3, 4, 2, 3.5, 5, 2.5),
        drifte = c(1.5, 1.2, 1, 1.8, 2, 1.1),
        gap = c(1, 1.1, 1, 0.8, 0.4, 0.6),
        ndt = c(0.25, 0.548, 0.2, 0.549, 0.1, 0.3),
        s = c(1, 0.9, 1.1, 1, 0.5, 1.4),
        sp = c(0.05, 0.1, 1e-7, 0.06, 0.2, 0.3)
      ),
      family = list(dpars = c("mu", "driftc", "drifte", "gap", "ndt", "s", "sp"))
    ),
    class = "brmsprep"
  )
  vectorised <- .rdm_log_lik(1, prep, cat_names = c("driftc", "drifte"), n_cats = 2)
  per_draw <- vapply(seq_len(n_draws), function(k) {
    one <- prep
    one$ndraws <- 1L
    one$dpars <- lapply(prep$dpars, `[`, k)
    .rdm_log_lik(1, one, cat_names = c("driftc", "drifte"), n_cats = 2)
  }, numeric(1))
  expect_true(all(is.finite(vectorised)))
  expect_equal(vectorised, per_draw, tolerance = 1e-12)

  # the same holds for the exported density with trial-varying parameters
  rt <- c(0.3, 0.552, 0.9, 2.5, 0.35, 0.7)
  response <- c(1L, 2L, 1L, 2L, 1L, 2L)
  vec <- drdm(rt, response, drift = c(3, 1.5),
              gap = prep$dpars$gap, sp = prep$dpars$sp, ndt = prep$dpars$ndt,
              s = prep$dpars$s, log = TRUE)
  each <- vapply(seq_along(rt), function(k) {
    drdm(rt[k], response[k], drift = c(3, 1.5), gap = prep$dpars$gap[k],
         sp = prep$dpars$sp[k], ndt = prep$dpars$ndt[k], s = prep$dpars$s[k], log = TRUE)
  }, numeric(1))
  expect_equal(vec, each, tolerance = 1e-12)
})

test_that("drdm() rejects a response outside the accumulators", {
  expect_error(drdm(0.5, 3, drift = c(3, 1.5), gap = 1, ndt = 0.2), "1:2")
  expect_error(drdm(0.5, 0, drift = c(3, 1.5), gap = 1, ndt = 0.2), "1:2")
})

test_that("rrdm returns valid data.frame", {
  dat <- rrdm(100, drift = c(3, 1.5), gap = 1, ndt = 0.2)
  expect_s3_class(dat, "data.frame")
  expect_true(all(c("rt", "response") %in% names(dat)))
  expect_true(all(dat$rt > 0.2))
  expect_true(all(dat$response %in% 1:2))
})

test_that("rrdm with sp > 0 returns valid data", {
  dat <- rrdm(100, drift = c(3, 1.5), gap = 0.7, sp = 0.3, ndt = 0.2)
  expect_true(all(dat$rt > 0.2))
})

test_that("validate_rdm_parameters catches invalid inputs", {
  expect_error(drdm(0.5, 1, drift = c(3, 1.5), gap = -1, ndt = 0.2))
  expect_error(drdm(0.5, 1, drift = c(-1, 1.5), gap = 1, ndt = 0.2))
  expect_error(drdm(0.5, 1, drift = c(3, 1.5), gap = 1, ndt = 0.2,
                    sp = -0.1))
  expect_error(drdm(0.5, 1, drift = c(3, 1.5), gap = 1, ndt = -1))
})

# -----------------------------------------------------------------------------
# Stage 4 numerics: R vs Stan parity, integration and inversion identities
# -----------------------------------------------------------------------------

# Extends "drdm() integrates to one over responses and time when s != 1" to a
# shrinking start-point range: every sp below crosses the midpoint switch
# (A* = 1e-4 * s * sqrt(t), see .dwald_full()) somewhere inside [0.2, Inf), so
# a value discontinuity introduced there would still break normalization.
test_that("drdm() integrates to one across a shrinking start-point range", {
  density <- function(rt, response, sp) {
    vapply(rt, function(x) {
      drdm(x, response, drift = c(3, 1.5), gap = 0.7, sp = sp, ndt = 0.2)
    }, numeric(1))
  }
  for (sp in c(1e-6, 1e-3, 0.05)) {
    total <- sum(vapply(1:2, function(response) {
      stats::integrate(density, 0.2, Inf, response = response, sp = sp,
                       rel.tol = 1e-9)$value
    }, numeric(1)))
    expect_equal(total, 1, tolerance = 1e-6)
  }
})

# A start-point range of 1e-12 is smaller than any sp a user would fit but is
# not exactly zero, so it must still take the small-A midpoint branch to
# numerical precision; a regression in that branch's cutoff would show up
# here first.
test_that("drdm(sp = 1e-12) matches drdm(sp = 0) to 1e-10", {
  rt <- c(0.3, 0.55, 0.9, 1.6)
  response <- c(1L, 2L, 1L, 2L)
  d0 <- drdm(rt, response, drift = c(3, 1.5), gap = 0.8, ndt = 0.2, sp = 0,
             log = TRUE)
  d_tiny <- drdm(rt, response, drift = c(3, 1.5), gap = 0.8, ndt = 0.2,
                 sp = 1e-12, log = TRUE)
  expect_equal(d_tiny, d0, tolerance = 1e-10)
})

# .rdm_race() is checked against drdm() directly elsewhere; this exercises the
# user-facing rrdm()/qrdm() pair instead, so a bug introduced only in the ndt
# shift or the response coding applied after .rdm_race() returns (both are
# rrdm()-only code) would still be caught. The quantile tolerance is the
# order-statistic delta-method SE, sqrt(p(1-p)/n) / f(q), with f(q) the
# marginal RT density (summed over responses) at the theoretical quantile.
test_that("rrdm() reproduces drdm()'s choice probabilities and RT quantiles", {
  withr::local_seed(23)
  drift <- c(3, 1.5, 1)
  gap <- 0.8
  sp <- 0.3
  ndt <- 0.2
  s <- 0.9
  n <- 2e5
  dat <- rrdm(n, drift = drift, gap = gap, sp = sp, ndt = ndt, s = s)

  for (j in seq_along(drift)) {
    p_sim <- mean(dat$response == j)
    se_p <- sqrt(p_sim * (1 - p_sim) / n)
    p_exact <- stats::integrate(function(rt) {
      vapply(rt, function(x) {
        drdm(x, j, drift = drift, gap = gap, sp = sp, ndt = ndt, s = s)
      }, numeric(1))
    }, ndt, Inf, rel.tol = 1e-8)$value
    expect_lt(abs(p_sim - p_exact), 4 * se_p)
  }

  probs <- c(0.1, 0.5, 0.9)
  q_theory <- qrdm(probs, drift = drift, gap = gap, sp = sp, ndt = ndt, s = s)
  q_sim <- stats::quantile(dat$rt, probs = probs, type = 7, names = FALSE)
  f_at_q <- vapply(seq_along(probs), function(k) {
    sum(vapply(seq_along(drift), function(j) {
      drdm(q_theory[k], j, drift = drift, gap = gap, sp = sp, ndt = ndt, s = s)
    }, numeric(1)))
  }, numeric(1))
  se_q <- sqrt(probs * (1 - probs) / n) / f_at_q
  expect_lt(max(abs(q_sim - q_theory) / se_q), 4)
})

# prdm() is assembled from .pwald_full()'s survival antiderivative directly
# (Tillman et al.'s Appendix A), never from drdm(); this checks the two
# constructions agree by independent numerical integration.
test_that("prdm() equals drdm() integrated and summed over responses", {
  drift <- c(2.5, 1.3)
  gap <- 0.9
  sp <- 0.25
  ndt <- 0.15
  s <- 0.8
  for (q in c(0.4, 0.9, 2.0)) {
    direct <- sum(vapply(seq_along(drift), function(j) {
      stats::integrate(function(rt) {
        vapply(rt, function(x) {
          drdm(x, j, drift = drift, gap = gap, sp = sp, ndt = ndt, s = s)
        }, numeric(1))
      }, ndt, q, rel.tol = 1e-10)$value
    }, numeric(1)))
    expect_equal(
      prdm(q, drift = drift, gap = gap, sp = sp, ndt = ndt, s = s),
      direct, tolerance = 1e-7
    )
  }
})

test_that("qrdm(prdm(q)) recovers q", {
  drift <- c(2.5, 1.3)
  gap <- 0.9
  sp <- 0.25
  ndt <- 0.15
  s <- 0.8
  for (q in c(0.4, 0.9, 2.0)) {
    p <- prdm(q, drift = drift, gap = gap, sp = sp, ndt = ndt, s = s)
    q_back <- qrdm(p, drift = drift, gap = gap, sp = sp, ndt = ndt, s = s)
    expect_equal(q_back, q, tolerance = 1e-6)
  }
})

# Protects rdm_log_pdf()/rdm_log_surv() in inst/stan_chunks/rdm_functions.stan
# against the same log-space quadrature reference used for R's
# .dwald_full()/.pwald_full(); the reference functions are copied in (not
# sourced) because local/ is gitignored. cmdstanr::expose_functions() does not
# build on this machine, so a fixed_param generated-quantities program
# evaluates the compiled package chunk directly (see test-model_rdm-stability.R's
# header note and local/rdm_stan_grid.R for the full-size version of this grid).
# Mutant (applied locally, not left in the tree): flipping the sign of the
# "2 * drift * b / s2" term in rdm_log_surv's lEb left the density grid
# unaffected (lEb only feeds the survivor: max rel error 4.9e-10, unchanged)
# but broke the survival grid from max rel error 2.5e-10 to 0.53, far past the
# 1e-7 gate below.
test_that("Stan rdm_log_pdf()/rdm_log_surv() match the quadrature reference on a grid", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")

  lse <- function(a, b) {
    m <- pmax(a, b)
    out <- m + log(exp(a - m) + exp(b - m))
    out[is.infinite(m) & m < 0] <- -Inf
    out
  }
  lde <- function(a, b) {
    d <- b - a
    out <- ifelse(d > -0.6931472, a + log(-expm1(d)), a + log1p(-exp(d)))
    out[a <= b] <- -Inf
    out
  }
  log1mexp <- function(x) ifelse(x > -0.6931472, log(-expm1(x)), log1p(-exp(x)))
  lwald <- function(t, u, v, s) {
    log(u) - 0.5 * (log(2 * pi) + 2 * log(s) + 3 * log(t)) - (u - v * t)^2 / (2 * s^2 * t)
  }
  lwald_cdf <- function(t, u, v, s) {
    z1 <- (v * t - u) / (s * sqrt(t))
    z2 <- -(v * t + u) / (s * sqrt(t))
    lse(stats::pnorm(z1, log.p = TRUE), 2 * u * v / s^2 + stats::pnorm(z2, log.p = TRUE))
  }
  lwald_surv <- function(t, u, v, s) {
    z1 <- (v * t - u) / (s * sqrt(t))
    z2 <- -(v * t + u) / (s * sqrt(t))
    lc <- lwald_cdf(t, u, v, s)
    ifelse(
      lc < log(0.5), log1mexp(lc),
      lde(stats::pnorm(-z1, log.p = TRUE), 2 * u * v / s^2 + stats::pnorm(z2, log.p = TRUE))
    )
  }
  simpson_log <- function(lf, lo, hi, n = 4001L) {
    u <- seq(lo, hi, length.out = n)
    y <- lf(u)
    y[!is.finite(y)] <- -Inf
    m <- max(y)
    if (!is.finite(m)) return(-Inf)
    w <- c(1, rep(c(4, 2), (n - 3) / 2), 4, 1)
    h <- (hi - lo) / (n - 1)
    m + log(sum(w * exp(y - m)) * h / 3)
  }
  int_log <- function(lf, lo, hi) {
    y <- lf(seq(lo, hi, length.out = 401))
    y[!is.finite(y)] <- -Inf
    m <- max(y)
    if (!is.finite(m)) return(-Inf)
    from_lo <- y[1] >= y[length(y)]
    g <- function(x) {
      u <- if (from_lo) lo + exp(x) else hi - exp(x)
      z <- exp(lf(u) - m + x)
      z[!is.finite(z)] <- 0
      z
    }
    r <- tryCatch(
      stats::integrate(g, -60, log(hi - lo), rel.tol = 1e-12, abs.tol = 0,
                       subdivisions = 5000L, stop.on.error = FALSE),
      error = function(e) NULL
    )
    if (is.null(r) || r$message != "OK" || !is.finite(r$value) || r$value <= 0) {
      return(simpson_log(lf, lo, hi))
    }
    m + log(r$value)
  }
  ref_lpdf <- function(t, v, gap, A, s) {
    if (A == 0) return(lwald(t, gap, v, s))
    int_log(function(u) lwald(t, u, v, s), gap, gap + A) - log(A)
  }
  ref_lsurv <- function(t, v, gap, A, s) {
    if (A == 0) return(lwald_surv(t, gap, v, s))
    logC <- int_log(function(u) lwald_cdf(t, u, v, s), gap, gap + A) - log(A)
    if (is.finite(logC) && logC < log(0.5)) return(log1mexp(logC))
    int_log(function(u) lwald_surv(t, u, v, s), gap, gap + A) - log(A)
  }

  grid <- expand.grid(
    t = c(1e-4, 1e-2, 0.3, 2, 20),
    drift = c(0.5, 3, 8),
    A = c(0, 1e-5, 1e-3, 0.3, 1),
    s = c(0.5, 1, 1.5)
  )
  grid$gap <- 1
  n_grid <- nrow(grid)

  sc <- system.file("stan_chunks", package = "bmm")
  funs <- paste(
    read_lines2(file.path(sc, "cswald_helper_functions.stan")),
    read_lines2(file.path(sc, "rdm_functions.stan")),
    sep = "\n"
  )
  code <- paste0(
    "functions {\n", funs, "\n}\n",
    "data { int G; vector[G] t; vector[G] v; vector[G] gap; vector[G] A; vector[G] s; }\n",
    "generated quantities {\n  vector[G] ld; vector[G] ls;\n",
    "  for (i in 1:G) {\n",
    "    if (A[i] == 0) { ld[i] = swald_lpdf(t[i] | v[i], gap[i], 0, s[i]);",
    " ls[i] = swald_lccdf(t[i] | v[i], gap[i], 0, s[i]); }\n",
    "    else { ld[i] = rdm_log_pdf(t[i], v[i], gap[i], A[i], s[i]);",
    " ls[i] = rdm_log_surv(t[i], v[i], gap[i], A[i], s[i]); }\n",
    "  }\n}\n"
  )
  # one compiled program for this test file
  mod <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(code), quiet = TRUE)
  fit <- mod$sample(
    data = list(G = n_grid, t = grid$t, v = grid$drift, gap = grid$gap,
               A = grid$A, s = grid$s),
    fixed_param = TRUE, chains = 1, iter_sampling = 1, iter_warmup = 0,
    refresh = 0, show_messages = FALSE, seed = 1, sig_figs = 18
  )
  draws <- fit$draws(format = "matrix")
  nm <- colnames(draws)
  draws <- as.numeric(draws[1, ])
  ld <- draws[match(paste0("ld[", seq_len(n_grid), "]"), nm)]
  ls <- draws[match(paste0("ls[", seq_len(n_grid), "]"), nm)]

  ref_d <- mapply(ref_lpdf, grid$t, grid$drift, grid$gap, grid$A, grid$s)
  ref_s <- mapply(ref_lsurv, grid$t, grid$drift, grid$gap, grid$A, grid$s)

  expect_true(all(is.finite(ld)))
  expect_true(all(is.finite(ls)))
  expect_lt(max(abs(ld - ref_d) / pmax(1, abs(ref_d))), 1e-7)
  expect_lt(max(abs(ls - ref_s) / pmax(1, abs(ref_s))), 1e-7)
})
