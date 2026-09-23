# =============================================================================
# Tests for cswald model (model-specific tests)
# Distribution function tests are in test-distributions.R
# =============================================================================

# -----------------------------------------------------------------------------
# Model construction tests
# -----------------------------------------------------------------------------

test_that("cswald() creates model with correct structure", {
  model <- cswald(rt = "rt", response = "response", version = "simple")

  expect_s3_class(model, "bmmodel")
  expect_s3_class(model, "cswald")
  expect_s3_class(model, "cswald_simple")
  expect_equal(model$resp_vars$rt, "rt")
  expect_equal(model$resp_vars$response, "response")
  expect_equal(model$version, "simple")
})

test_that("cswald() creates crisk version correctly", {
  model <- cswald(rt = "rt", response = "response", version = "crisk")

  expect_s3_class(model, "cswald_crisk")
  expect_equal(model$version, "crisk")
  expect_true("zr" %in% names(model$parameters))
})

test_that("cswald simple version has correct parameters", {
  model <- cswald(rt = "rt", response = "response", version = "simple")

  expect_true(all(c("drift", "bound", "ndt", "s") %in% names(model$parameters)))
  expect_false("zr" %in% names(model$parameters))
})

test_that("cswald crisk version has correct parameters", {
  model <- cswald(rt = "rt", response = "response", version = "crisk")

  expect_true(all(c("drift", "bound", "ndt", "zr", "s") %in% names(model$parameters)))
})

test_that("cswald accepts custom links", {
  model <- cswald(rt = "rt", response = "response",
                  links = list(drift = "identity"), version = "simple")
  expect_equal(model$links$drift, "identity")
})

# -----------------------------------------------------------------------------
# Data validation tests (check_data.cswald)
# -----------------------------------------------------------------------------

test_that("check_data.cswald errors when required variables missing", {
  model <- cswald(rt = "rt", response = "response")

  expect_error(
    check_data(model, data.frame(x = 1), bmf(drift ~ 1)),
    "RT variable 'rt' is not present"
  )

  expect_error(
    check_data(model, data.frame(rt = 1), bmf(drift ~ 1)),
    "response variable 'response' is not present"
  )
})

test_that("check_data.cswald errors when RT contains NA", {
  model <- cswald(rt = "rt", response = "response")
  dat <- data.frame(rt = c(0.5, NA, 0.8), response = c(1, 1, 0))

  expect_error(
    check_data(model, dat, bmf(drift ~ 1)),
    "RT variable 'rt' contains.*NA"
  )
})

test_that("check_data.cswald errors when response contains NA", {
  model <- cswald(rt = "rt", response = "response")
  dat <- data.frame(rt = c(0.5, 0.6, 0.8), response = c(1, NA, 0))

  expect_error(
    check_data(model, dat, bmf(drift ~ 1)),
    "response variable 'response' contains.*NA"
  )
})

test_that("check_data.cswald errors when RT contains negative values", {
  model <- cswald(rt = "rt", response = "response")
  dat <- data.frame(rt = c(-0.5, 0.6, 0.8), response = c(1, 1, 0))

  expect_error(
    check_data(model, dat, bmf(drift ~ 1)),
    "reaction times are lower than zero"
  )
})

test_that("check_data.cswald warns when RT > 10 seconds", {
  # Use crisk version and large dataset to avoid other warnings
  model <- cswald(rt = "rt", response = "response", version = "crisk")
  dat <- data.frame(
    rt = c(runif(99, 0.4, 1.5), 15),
    response = rep(c(1, 0), 50)
  )

  expect_warning(
    check_data(model, dat, bmf(drift ~ 1)),
    "larger than 10 seconds"
  )
})

test_that("check_data.cswald warns when RT < 0.1 seconds", {
  # Use crisk version and large dataset to avoid other warnings
  model <- cswald(rt = "rt", response = "response", version = "crisk")
  dat <- data.frame(
    rt = c(0.05, runif(99, 0.4, 1.5)),
    response = rep(c(1, 0), 50)
  )

  expect_warning(
    check_data(model, dat, bmf(drift ~ 1)),
    "smaller than 0.100 seconds"
  )
})

test_that("check_data.cswald handles different response formats", {
  # Use crisk version to avoid error rate warning, and use larger datasets
  model <- cswald(rt = "rt", response = "response", version = "crisk")

  # integer 0/1 - should work silently
  dat_int <- data.frame(
    rt = runif(100, 0.4, 1.5),
    response = rep(c(1L, 0L), 50)
  )
  expect_silent(check_data(model, dat_int, bmf(drift ~ 1)))

  # logical - should warn and convert
  dat_logical <- data.frame(
    rt = runif(100, 0.4, 1.5),
    response = rep(c(TRUE, FALSE), 50)
  )
  expect_warning(
    check_data(model, dat_logical, bmf(drift ~ 1)),
    "boolean"
  )

  # character upper/lower - should warn and convert
  dat_char <- data.frame(
    rt = runif(100, 0.4, 1.5),
    response = rep(c("upper", "lower"), 50)
  )
  expect_warning(
    check_data(model, dat_char, bmf(drift ~ 1)),
    "character"
  )

  # factor - should warn and convert
  dat_factor <- data.frame(
    rt = runif(100, 0.4, 1.5),
    response = factor(rep(c("upper", "lower"), 50))
  )
  expect_warning(
    check_data(model, dat_factor, bmf(drift ~ 1)),
    "character"
  )
})

test_that("check_data.cswald errors on invalid response values", {
  # Use crisk version to avoid error rate warning
  model <- cswald(rt = "rt", response = "response", version = "crisk")

  # invalid integer values
  dat <- data.frame(
    rt = runif(100, 0.4, 1.5),
    response = c(rep(1, 50), rep(2, 50))
  )
  expect_error(
    check_data(model, dat, bmf(drift ~ 1)),
    "values other than 0 and 1"
  )

  # invalid character values
  dat_char <- data.frame(
    rt = runif(100, 0.4, 1.5),
    response = c(rep("upper", 50), rep("invalid", 50))
  )
  expect_error(
    check_data(model, dat_char, bmf(drift ~ 1)),
    "invalid character values"
  )
})

test_that("check_data.cswald warns about high error rate for simple version", {
  model <- cswald(rt = "rt", response = "response", version = "simple")

  # 30% error rate
  dat <- data.frame(
    rt = runif(100, 0.4, 1.5),
    response = c(rep(0, 30), rep(1, 70))
  )

  expect_warning(
    check_data(model, dat, bmf(drift ~ 1)),
    "error rate"
  )
})

test_that("check_data.cswald does not warn about error rate for crisk version", {
  model <- cswald(rt = "rt", response = "response", version = "crisk")

  # 30% error rate - no warning for crisk
  dat <- data.frame(
    rt = runif(100, 0.4, 1.5),
    response = c(rep(0, 30), rep(1, 70))
  )

  expect_silent(check_data(model, dat, bmf(drift ~ 1)))
})

test_that("check_data.cswald returns a data.frame", {
  # Use crisk version to avoid error rate warning
  model <- cswald(rt = "rt", response = "response", version = "crisk")
  dat <- data.frame(rt = runif(100, 0.4, 1.5), response = rep(c(0, 1), 50))

  result <- check_data(model, dat, bmf(drift ~ 1))
  expect_s3_class(result, "data.frame")
})

# -----------------------------------------------------------------------------
# Formula conversion tests (bmf2bf.cswald)
# -----------------------------------------------------------------------------

test_that("bmf2bf.cswald creates correct brms formula", {
  model <- cswald(rt = "rt", response = "response")
  formula <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1)

  bf <- bmf2bf(model, formula)

  expect_s3_class(bf, "brmsformula")
  # check that response includes dec() term
  expect_true(grepl("dec", deparse(bf$formula)))
})

# -----------------------------------------------------------------------------
# Model configuration tests
# -----------------------------------------------------------------------------

test_that("configure_model.cswald_simple returns correct components", {
  skip_on_cran()

  model <- cswald(rt = "rt", response = "response", version = "simple")
  dat <- data.frame(
    rt = runif(100, 0.4, 1.5),
    response = sample(c(0, 1), 100, replace = TRUE)
  )
  formula <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1)

  config <- configure_model(model, dat, formula)

  expect_true(all(c("formula", "data", "stanvars") %in% names(config)))
  expect_s3_class(config$formula, "brmsformula")
  expect_s3_class(config$formula$family, "customfamily")
  expect_equal(config$formula$family$name, "cswald")
})

test_that("configure_model.cswald_crisk returns correct components", {
  skip_on_cran()

  model <- cswald(rt = "rt", response = "response", version = "crisk")
  dat <- data.frame(
    rt = runif(100, 0.4, 1.5),
    response = sample(c(0, 1), 100, replace = TRUE)
  )
  formula <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1, zr ~ 1)

  config <- configure_model(model, dat, formula)

  expect_true(all(c("formula", "data", "stanvars") %in% names(config)))
  expect_s3_class(config$formula, "brmsformula")
  expect_s3_class(config$formula$family, "customfamily")
  expect_equal(config$formula$family$name, "cswald_crisk")
})

cswald_data <- function(n = 100) {
  data.frame(
    rt = runif(n, 0.4, 1.5),
    response = sample(c(0, 1), n, replace = TRUE)
  )
}

test_that("cswald slices its decision variable exactly when brms slices the data", {
  skip_on_cran()

  dat <- cswald_data()
  model <- cswald(rt = "rt", response = "response", version = "simple")
  formula <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1)

  expect_equal(configure_model(model, dat, formula)$formula$family$vars, "dec")

  withr::with_options(list(brms.threads = brms::threading(2)), {
    expect_equal(
      configure_model(model, dat, formula)$formula$family$vars,
      "dec[start:end]"
    )
  })

  # force = TRUE tells brms to compile with threading but leave the generated
  # code alone, so start/end are never defined and slicing would not compile
  withr::with_options(list(brms.threads = brms::threading(2, force = TRUE)), {
    expect_equal(configure_model(model, dat, formula)$formula$family$vars, "dec")
  })

  # brms accepts a bare number for the option
  withr::with_options(list(brms.threads = 2), {
    expect_equal(
      configure_model(model, dat, formula)$formula$family$vars,
      "dec[start:end]"
    )
  })
})

# stanc accepts a call with two vector dpars swapped, and the parity test builds
# its own call, so the whole generated call is pinned here for both versions
cswald_generated_calls <- list(
  simple = c(
    serial = "cswald_lpdf(Y | mu, drift, bound, ndt, s, sndt, dec)",
    threaded = "cswald_lpdf(Y[start:end] | mu, drift, bound, ndt, s, sndt, dec[start:end])"
  ),
  crisk = c(
    serial = "cswald_crisk_lpdf(Y | mu, drift, bound, ndt, zr, s, sndt, dec)",
    threaded = "cswald_crisk_lpdf(Y[start:end] | mu, drift, bound, ndt, zr, s, sndt, dec[start:end])"
  )
)

test_that("the generated cswald call pairs the sliced response with the sliced decisions", {
  skip_on_cran()

  dat <- cswald_data()
  formula <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1)

  for (version in names(cswald_generated_calls)) {
    model <- cswald(rt = "rt", response = "response", version = version)
    expect_match(
      suppressWarnings(stancode(formula, dat, model)),
      cswald_generated_calls[[version]][["serial"]],
      fixed = TRUE
    )
    expect_match(
      suppressWarnings(stancode(formula, dat, model, threads = brms::threading(2))),
      cswald_generated_calls[[version]][["threaded"]],
      fixed = TRUE
    )
  }
})

test_that("an explicit threads = NULL beats a global threading option", {
  skip_on_cran()

  dat <- cswald_data()
  formula <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1)

  # brms reads threads = NULL as "threading off" and generates serial code, so
  # the family must not slice its decisions even though the option is set
  withr::local_options(brms.threads = brms::threading(2))
  for (version in names(cswald_generated_calls)) {
    code <- suppressWarnings(stancode(
      formula, dat, cswald(rt = "rt", response = "response", version = version),
      threads = NULL
    ))
    expect_match(code, cswald_generated_calls[[version]][["serial"]], fixed = TRUE)
    expect_false(grepl("start:end", code, fixed = TRUE))
  }
})

test_that("threading a non-zero sndt indexes the decisions with the sliced row", {
  skip_on_cran()

  dat <- cswald_data()
  # a non-zero sndt drops the family back to loop = TRUE and vars = "dec[n]";
  # brms rewrites that n to the nn it also indexes Y with, so the two stay paired
  formula <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1, sndt ~ 1)

  for (version in names(cswald_generated_calls)) {
    code <- suppressWarnings(stancode(
      formula, dat, cswald(rt = "rt", response = "response", version = version),
      threads = brms::threading(2)
    ))
    expect_match(code, "Y[nn] | mu[n]", fixed = TRUE)
    expect_match(code, "dec[nn]);", fixed = TRUE)
    expect_false(grepl("dec[n])", code, fixed = TRUE))
  }
})

# -----------------------------------------------------------------------------
# Integration tests with mock backend
# -----------------------------------------------------------------------------

test_that("cswald simple version runs with mock backend", {
  skip_on_cran()

  dat <- rcswald(n = 100, drift = 2, bound = 1.5, ndt = 0.3)
  model <- cswald(rt = "rt", response = "response", version = "simple")
  formula <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1)

  expect_silent(
    bmm(formula, dat, model, backend = "mock", mock = 1, rename = FALSE)
  )
})

test_that("cswald crisk version runs with mock backend", {
  skip_on_cran()

  dat <- rcswald(n = 100, drift = 2, bound = 1.5, ndt = 0.3, zr = 0.5)
  model <- cswald(rt = "rt", response = "response", version = "crisk")
  formula <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1, zr ~ 1)

  expect_silent(
    bmm(formula, dat, model, backend = "mock", mock = 1, rename = FALSE)
  )
})

test_that("cswald simple version with predictors runs with mock backend", {
  skip_on_cran()

  dat <- rcswald(n = 200, drift = 2, bound = 1.5, ndt = 0.3)
  dat$condition <- rep(c("A", "B"), each = 100)
  model <- cswald(rt = "rt", response = "response", version = "simple")
  formula <- bmf(drift ~ condition, bound ~ 1, ndt ~ 1)

  expect_silent(
    bmm(formula, dat, model, backend = "mock", mock = 1, rename = FALSE)
  )
})

test_that("cswald crisk version allows negative drift (more lower responses)", {
  skip_on_cran()

  # Generate data with negative drift (biased toward lower boundary)
  dat <- rcswald(n = 100, drift = -1.5, bound = 1.5, ndt = 0.3, zr = 0.5)

  # Most responses should be lower (0) with negative drift
  expect_true(mean(dat$response == 0) > 0.5)

  model <- cswald(rt = "rt", response = "response", version = "crisk")
  formula <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1, zr ~ 1)

  # Model should run without error
  expect_silent(
    bmm(formula, dat, model, backend = "mock", mock = 1, rename = FALSE)
  )
})

test_that("cswald handles all-correct responses", {
  skip_on_cran()

  # Data with all correct responses
  dat <- data.frame(
    rt = runif(100, 0.4, 1.5),
    response = rep(1, 100)
  )

  model <- cswald(rt = "rt", response = "response", version = "simple")
  formula <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1)

  # Should run without error (though may warn about 0% error rate)
  expect_silent(
    bmm(formula, dat, model, backend = "mock", mock = 1, rename = FALSE)
  )
})

test_that("cswald handles high error rate data with crisk version", {
  skip_on_cran()

  # Data with 50% error rate
  dat <- data.frame(
    rt = runif(100, 0.4, 1.5),
    response = rep(c(0, 1), 50)
  )

  model <- cswald(rt = "rt", response = "response", version = "crisk")
  formula <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1, zr ~ 1)

  expect_silent(
    bmm(formula, dat, model, backend = "mock", mock = 1, rename = FALSE)
  )
})

# -----------------------------------------------------------------------------
# R <-> Stan parity of the likelihood
# -----------------------------------------------------------------------------

# Builds a standalone Stan program exposing both overloads of the cswald
# likelihood as generated quantities, so the vectorized (loop = FALSE) form that
# brms actually compiles can be compared against the scalar form and against the
# R implementation used by log_lik() and posterior_predict()
compile_cswald_parity_model <- function(version) {
  sc_path <- system.file("stan_chunks", package = "bmm")
  chunk <- function(f) read_lines2(file.path(sc_path, f))
  crisk <- version == "crisk"
  lpdf <- if (crisk) "cswald_crisk_lpdf" else "cswald_lpdf"
  scalar_args <- if (crisk) {
    "mu[n], drift[n], bound[n], ndt[n], zr[n], s[n], sndt[n], dec[n]"
  } else {
    "mu[n], drift[n], bound[n], ndt[n], s[n], sndt[n], dec[n]"
  }
  vector_args <- if (crisk) {
    "mu, drift, bound, ndt, zr, s, sndt, dec"
  } else {
    "mu, drift, bound, ndt, s, sndt, dec"
  }
  code <- glue(
    "functions {{\n{chunk('cswald_helper_functions.stan')}\n",
    "{chunk(paste0('cswald_', version, '_functions.stan'))}\n}}\n",
    "data {{\n  int N;\n  vector[N] rt;\n  array[N] int dec;\n  vector[N] mu;\n",
    "  vector[N] drift;\n  vector[N] bound;\n  vector[N] ndt;\n  vector[N] s;\n",
    "  vector<lower=0,upper=1>[N] zr;\n  vector<lower=0>[N] sndt;\n}}\n",
    "generated quantities {{\n",
    "  real lp_vector = {lpdf}(rt | {vector_args});\n",
    "  vector[N] lp_scalar;\n",
    "  for (n in 1:N) lp_scalar[n] = {lpdf}(rt[n] | {scalar_args});\n}}\n"
  )
  file <- file.path(tempdir(), glue("bmm_cswald_parity_{version}.stan"))
  writeLines(code, file)
  cmdstanr::cmdstan_model(file)
}

# Random parameters and observations for the parity test, plus 30 censored
# observations in the two regimes where the survivor leaves its vectorized
# probability-space path: 2*bound*drift/s^2 in the hundreds, and the deep tail
# where Phi(-z1) underflows to 0
cswald_parity_data <- function(n = 200, sndt = 0) {
  drift <- runif(n, 0.5, 4)
  bound <- runif(n, 0.5, 2)
  ndt <- runif(n, 0.05, 0.2)
  rt <- ndt + runif(n, 0.05, 2.5)

  drift_ext <- c(runif(15, 15, 25), runif(15, 20, 30))
  bound_ext <- c(runif(15, 15, 25), runif(15, 0.3, 0.6))
  ndt_ext <- runif(30, 0.05, 0.2)
  rt_ext <- ndt_ext + c(
    bound_ext[1:15] / drift_ext[1:15] * runif(15, 0.8, 1.2),
    runif(15, 2, 3)
  )

  list(
    N = n + 30,
    rt = c(rt, rt_ext),
    dec = c(sample(c(0L, 1L), n, replace = TRUE), rep(0L, 30)),
    mu = rep(0, n + 30),
    drift = c(drift, drift_ext),
    bound = c(bound, bound_ext),
    ndt = c(ndt, ndt_ext),
    s = c(runif(n, 0.7, 1.3), runif(30, 0.9, 1.1)),
    zr = c(runif(n, 0.2, 0.8), runif(30, 0.3, 0.7)),
    sndt = rep(sndt, n + 30)
  )
}

test_that("the vectorized cswald likelihood matches the scalar and R versions", {
  skip_on_cran()
  skip_if_not(requireNamespace("cmdstanr", quietly = TRUE), "cmdstanr not available")
  skip_if_not(nzchar(Sys.getenv("CMDSTAN", unset = "")) ||
    !is.null(tryCatch(cmdstanr::cmdstan_path(), error = function(e) NULL)),
  "CmdStan not installed"
  )

  # seeded so a failure reproduces. Under this seed the log-space fallback,
  # with both clamp conditions firing, takes 23 survivor terms for simple and 7
  # for crisk. Over seeds 1-200 the floor was 13 for simple and 6 for crisk
  # (below 13 in 170 of 200 seeds); neither version dropped to 0. Re-measuring
  # needs this loop's draw order, not two cswald_parity_data() calls in a row:
  # $sample() consumes two uniforms, so crisk draws from a shifted stream
  withr::local_seed(20260921)

  for (version in c("simple", "crisk")) {
    model <- compile_cswald_parity_model(version)
    # sndt = 0 is the closed form the vectorized overload evaluates; sndt > 0
    # exercises the convolution, including the strip where rt - ndt < sndt
    for (sndt in c(0, 0.1, 0.3)) {
      sdata <- cswald_parity_data(sndt = sndt)
      fit <- model$sample(
        data = sdata, chains = 1, iter_sampling = 1, fixed_param = TRUE,
        refresh = 0, show_messages = FALSE, sig_figs = 18, seed = 1
      )
      lp_vector <- as.numeric(fit$draws("lp_vector", format = "draws_matrix")[1, 1])
      lp_scalar <- as.numeric(fit$draws("lp_scalar", format = "draws_matrix")[1, ])
      lp_r <- with(sdata, .dcswald(rt, dec, drift, bound, ndt, zr, s, sndt,
        version = version, log = TRUE
      ))

      # the scalar overload is what brms compiles once sndt is in play; it must
      # match the R mirror used by log_lik()
      expect_equal(lp_scalar, lp_r, tolerance = 1e-8)
      # expect_equal() averages the relative difference over the differing
      # elements, so one bad element next to 229 good ones passes. The floor at
      # sndt = 0 is Stan's Phi(), accurate to ~1e-10 absolute against a 60-digit
      # reference. With sndt > 0 the censored terms are a difference quotient
      # of integrated survivors that swald_sndt_lccdf takes down to a relative
      # difference of 1e-8, so rounding in G is amplified up to 1e8: 1.1e-8
      # under this seed, 9.1e-8 over 60 seeds, all on censored deep-tail terms
      expect_lt(max(abs(lp_scalar - lp_r)), if (sndt == 0) 1e-8 else 1e-7)
      # the vectorized overload is what brms compiles at sndt = 0; it must
      # reproduce the scalar likelihood it replaced
      if (sndt == 0) expect_equal(lp_vector, sum(lp_scalar), tolerance = 1e-10)
    }
  }
})

# -----------------------------------------------------------------------------
# sndt: link resolution and family selection
# -----------------------------------------------------------------------------

test_that("cswald reports the link sndt is actually fixed on, and restores it", {
  model <- cswald(rt = "rt", response = "response")

  # a fixed value is a constant() on the link scale, so sndt = 0 is only
  # sndt = 0 under an identity link (s = 0 under its log link means s = 1)
  expect_equal(model$fixed_parameters$sndt, 0)
  expect_equal(model$links$sndt, "identity")
  expect_equal(model$links$s, "log")

  freed <- update_model_fixed_parameters(
    model, bmf(drift ~ 1, bound ~ 1, ndt ~ 1, sndt ~ 1)
  )
  expect_null(freed$fixed_parameters$sndt)
  expect_equal(freed$links$sndt, "log")

  refixed <- update_model_fixed_parameters(
    freed, bmf(drift ~ 1, bound ~ 1, ndt ~ 1, sndt = 0)
  )
  expect_equal(refixed$links$sndt, "identity")
})

test_that("a fixed sndt reaches the likelihood on the natural scale", {
  dat <- cswald_data()
  model <- cswald(rt = "rt", response = "response")
  formula <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1, sndt = 0.15)

  model <- check_model(model, dat, formula)
  config <- configure_model(model, dat, check_formula(model, dat, formula))
  prior <- configure_prior(model, dat, config$formula, NULL)

  # constant(0.15), not constant(log(0.15)) and not exp(0.15)
  expect_true(any(grepl("constant(0.15)", prior$prior, fixed = TRUE)))
  expect_equal(config$formula$family$link_sndt, "identity")
})

test_that("the vectorized family is chosen by the value of sndt, not its fixedness", {
  dat <- cswald_data()
  model <- cswald(rt = "rt", response = "response")
  configure <- function(f) {
    m <- check_model(model, dat, f)
    configure_model(m, dat, check_formula(m, dat, f))$formula$family
  }

  # only sndt == 0 collapses the likelihood to the closed form the vectorized
  # overload evaluates; a non-zero constant is still a convolution
  expect_false(configure(bmf(drift ~ 1, bound ~ 1, ndt ~ 1))$loop)
  expect_true(configure(bmf(drift ~ 1, bound ~ 1, ndt ~ 1, sndt = 0.15))$loop)
  expect_true(configure(bmf(drift ~ 1, bound ~ 1, ndt ~ 1, sndt ~ 1))$loop)
})

test_that("the cswald survivor keeps a finite gradient where Phi() nears underflow", {
  skip_on_cran()
  skip_if_not(requireNamespace("cmdstanr", quietly = TRUE), "cmdstanr not available")
  skip_if_not(nzchar(Sys.getenv("CMDSTAN", unset = "")) ||
    !is.null(tryCatch(cmdstanr::cmdstan_path(), error = function(e) NULL)),
  "CmdStan not installed"
  )

  sc_path <- system.file("stan_chunks", package = "bmm")
  file <- file.path(tempdir(), "bmm_cswald_surv_grad.stan")
  writeLines(c(
    "functions {", read_lines2(file.path(sc_path, "cswald_helper_functions.stan")), "}",
    "data { int N; vector[N] t; vector[N] drift0; vector[N] bound; vector[N] s; }",
    "parameters { real shift; }",
    "model { target += sum(swald_log_surv_vec(t, drift0 + shift, bound, s)); }"
  ), file)
  model <- cmdstanr::cmdstan_model(file)

  # Phi(-z1) sits just above its underflow at -37.5: 1 / Phi(z2) overflows in the
  # first case, and the survivor difference is subnormal in the second. Both go
  # through the recompute loop, so a shared parameter inherits their gradient
  sdata <- list(
    N = 3, t = c(2.7870782702230001, 0.5, 1),
    drift0 = c(22.120998953469002, 54.447222151364002, 2),
    bound = c(0.44479261357337002, 0.70710678118655002, 1),
    s = c(0.99227269766852, 1, 1)
  )
  # CmdStan rejects a point with a non-finite gradient, which makes cmdstanr
  # error; the huge `error` threshold keeps finite-difference noise in this deep
  # tail from failing the run on its own
  grad <- model$diagnose(
    data = sdata, init = list(list(shift = 0)), error = 1e6, seed = 1
  )$gradients()

  expect_true(all(is.finite(grad$model)))
})

test_that("swald_log_Phi keeps a finite gradient on its own where Phi() nears underflow", {
  skip_on_cran() # compiles a Stan program
  skip_if_not(requireNamespace("cmdstanr", quietly = TRUE), "cmdstanr not available")
  skip_if_not(nzchar(Sys.getenv("CMDSTAN", unset = "")) ||
    !is.null(tryCatch(cmdstanr::cmdstan_path(), error = function(e) NULL)),
  "CmdStan not installed"
  )

  sc_path <- system.file("stan_chunks", package = "bmm")
  file <- file.path(tempdir(), "bmm_cswald_log_phi_grad.stan")
  writeLines(c(
    "functions {", read_lines2(file.path(sc_path, "cswald_helper_functions.stan")), "}",
    "data { real z0; real weight; }",
    "parameters { real shift; }",
    "model { target += weight * swald_log_Phi(z0 + shift); }"
  ), file)
  model <- cmdstanr::cmdstan_model(file)

  # swald_log_surv redirects these points before swald_log_Phi sees them, so
  # the test above cannot tell whether swald_log_Phi's own crossover holds; the
  # callers that #406 adds are not behind that guard. Stan's Phi(-37.3) is
  # 8.2e-305, below the 1e-300 crossover and above the -37.5 cutoff where Phi()
  # returns 0. The weight stands in for the adjoint the callers pass down,
  # which swald_log_surv scales by up to 1 / 1e-300: with a unit adjoint,
  # 1 / Phi(z) still fits in a double and log(Phi(z)) would pass here. The
  # multiplier that overflows scales with Phi(z), so this z is a mild case:
  # 1.5e4 here against about 8 for a caller landing at the -37.5 cutoff
  grad <- model$diagnose(
    data = list(z0 = -37.3, weight = 1e6), # adjoint stand-in
    init = list(list(shift = 0)),
    error = 1e6, # finite-difference tolerance, unrelated to the weight
    seed = 1
  )$gradients()

  expect_true(all(is.finite(grad$model)))
})
