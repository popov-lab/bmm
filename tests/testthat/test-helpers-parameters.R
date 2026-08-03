test_that("k2sd works", {
  # Test format
  kappa <- runif(10)
  expect_length(k2sd(kappa), 10)
  expect_type(k2sd(kappa), "double")
  expect_length(k2sd(2), 1)
  expect_type(k2sd(2), "double")

  # Test edge cases
  expect_equal(k2sd(0), Inf)
  expect_equal(k2sd(Inf), 0)

  # Test vector of edge cases
  expect_equal(k2sd(c(0, Inf)), c(Inf, 0))

  # Test known values (compared to pre-computed results)
  expect_equal(k2sd(1), 1.270088, tolerance = 1e-6)
  expect_equal(k2sd(10), 0.3248638, tolerance = 1e-6)

  # Test NA handling
  expect_true(is.na(k2sd(NA)))
  expect_true(all(is.na(k2sd(c(1, NA, 3))[2])))

  # Test invalid inputs
  expect_error(k2sd("a"))
  expect_error(k2sd(NULL))
})

test_that("c parameter conversions work", {
  # Test basic conversion
  c_sqrtexp <- 4
  kappa <- 3
  c_bessel <- c_sqrtexp2bessel(c_sqrtexp, kappa)
  expect_equal(c_bessel2sqrtexp(c_bessel, kappa), c_sqrtexp)

  # Test vector inputs
  c_vec <- c(1, 2, 3)
  kappa_vec <- c(2, 3, 4)
  c_bessel_vec <- c_sqrtexp2bessel(c_vec, kappa_vec)
  expect_equal(c_bessel2sqrtexp(c_bessel_vec, kappa_vec), c_vec)

  # Test error handling
  expect_error(c_sqrtexp2bessel(-1, 2), "c must be non-negative")
  expect_error(c_bessel2sqrtexp(1, -2), "kappa must be non-negative")
})

test_that("identity is a no-op and round-trips", {
  x <- c(-2, -0.5, 0, 1.25, 3)
  eta <- link_transform(x, "identity", inverse = FALSE)
  back <- link_transform(eta, "identity", inverse = TRUE)
  expect_identical(eta, x)
  expect_identical(back, x)
})

test_that("log and inverse log round-trip on positive values", {
  x <- c(0.1, 0.5, 1, 2, 10)
  eta <- link_transform(x, "log", inverse = FALSE)
  back <- link_transform(eta, "log", inverse = TRUE)
  expect_equal(back, x, tolerance = 1e-12)
})

test_that("softplus and inverse softplus round-trip on positive values", {
  x <- c(0.5, 2, 5)
  eta <- link_transform(x, "softplus", inverse = FALSE)
  back <- link_transform(eta, "softplus", inverse = TRUE)
  expect_equal(back, x, tolerance = 1e-12)
})

test_that("softplus inverse maps reals to positive values", {
  eta <- c(-5, -1, 0, 2, 10)
  x <- link_transform(eta, "softplus", inverse = TRUE)
  expect_true(all(is.finite(x)))
  expect_true(all(x > 0))
})

test_that("log1p and inverse expm1 round-trip for x > -1", {
  x <- c(-0.5, -0.1, 0, 0.3, 5)
  eta <- link_transform(x, "log1p", inverse = FALSE)
  back <- link_transform(eta, "log1p", inverse = TRUE)
  expect_equal(back, x, tolerance = 1e-12)
})

test_that("logm1 (brms::logm1) round-trip for x > 1", {
  testthat::skip_if_not_installed("brms")
  x <- c(1.1, 2, 5, 10)
  eta <- link_transform(x, "logm1", inverse = FALSE)
  back <- link_transform(eta, "logm1", inverse = TRUE)
  expect_equal(back, x, tolerance = 1e-12)
})

test_that("inverse link is its own inverse (x != 0)", {
  x <- c(-3, -0.5, 0.2, 4)
  eta <- link_transform(x, "inverse", inverse = FALSE)
  back <- link_transform(eta, "inverse", inverse = TRUE)
  expect_equal(back, x, tolerance = 1e-12)
})

test_that("sqrt link round-trip for x >= 0", {
  x <- c(0, 0.01, 1, 2, 9)
  eta <- link_transform(x, "sqrt", inverse = FALSE)
  back <- link_transform(eta, "sqrt", inverse = TRUE)
  expect_equal(back, x, tolerance = 1e-12)
})

test_that("logit round-trip for probabilities in (0,1)", {
  p <- c(0.1, 0.25, 0.5, 0.75, 0.9)
  eta <- link_transform(p, "logit", inverse = FALSE)
  back <- link_transform(eta, "logit", inverse = TRUE)
  expect_equal(back, p, tolerance = 1e-12)
})

test_that("probit round-trip for probabilities in (0,1)", {
  p <- c(0.01, 0.2, 0.5, 0.8, 0.99)
  eta <- link_transform(p, "probit", inverse = FALSE)
  back <- link_transform(eta, "probit", inverse = TRUE)
  expect_equal(back, p, tolerance = 1e-12)
})

test_that("tan_half round-trip on a safe interval (-pi, pi)", {
  x <- c(-2, -1, 0, 1, 2)  # well within (-pi, pi)
  eta <- link_transform(x, "tan_half", inverse = FALSE)   # tan(x/2)
  back <- link_transform(eta, "tan_half", inverse = TRUE) # 2*atan(eta)
  expect_equal(back, x, tolerance = 1e-12)
})


test_that("cloglog round-trip for probabilities in (0,1)", {
  p <- seq(0.1,0.9,by = 0.1)
  eta <- link_transform(p, "cloglog", inverse = FALSE)
  back <- link_transform(eta, "cloglog", inverse = TRUE)
  expect_equal(back, p, tolerance = 1e-9)
})

test_that("loglog round-trip for probabilities in (0,1)", {
  p <- seq(0.1,0.9,by = 0.1)
  eta <- link_transform(p, "loglog", inverse = FALSE)     # log(-log(p))
  back <- link_transform(eta, "loglog", inverse = TRUE)   # exp(-exp(eta))
  expect_equal(back, p, tolerance = 1e-12)
})

test_that("loglog inverse maps reals to (0,1)", {
  eta <- c(-3, -2, 0, 1, 3)
  p <- link_transform(eta, "loglog", inverse = TRUE)
  expect_true(all(is.finite(p)))
  expect_true(all(p > 0 & p < 1))
})

test_that("loglog is monotone decreasing in p", {
  p <- c(0.1, 0.2, 0.4, 0.8)  # increasing p
  eta <- link_transform(p, "loglog", inverse = FALSE)
  # as p increases, eta decreases strictly
  expect_true(all(diff(eta) < 0))
})

test_that("loglog vectorization and NA propagation", {
  p <- c(0.2, NA_real_, 0.7)
  eta <- link_transform(p, "loglog", inverse = FALSE)
  back <- link_transform(eta, "loglog", inverse = TRUE)
  expect_true(is.na(eta[2]))
  expect_true(is.na(back[2]))
  expect_equal(back[c(1,3)], p[c(1,3)], tolerance = 1e-12)
})

test_that("vectorization works and NA positions are preserved", {
  p <- c(0.2, NA_real_, 0.8)
  eta <- link_transform(p, "logit", inverse = FALSE)
  back <- link_transform(eta, "logit", inverse = TRUE)
  expect_true(is.na(eta[2]))
  expect_true(is.na(back[2]))
  expect_equal(back[c(1,3)], p[c(1,3)], tolerance = 1e-12)
})

test_that("unknown link errors clearly", {
  expect_error(link_transform(1:3, "not_a_link"))
})

test_that("non-numeric values error", {
  expect_error(link_transform(c("a","b"), "log"))
})

test_that("NULL link is treated as identity", {
  x <- c(-2, -0.5, 0, 1.25, 3)
  eta <- link_transform(x, NULL, inverse = FALSE)
  back <- link_transform(eta, NULL, inverse = TRUE)
  expect_identical(eta, x)
  expect_identical(back, x)
})


# ===========================================================================
# .is_softmax_param()
# ===========================================================================

test_that(".is_softmax_param detects mixture3p softmax params", {
  mock_model <- structure(list(), class = c("mixture3p", "bmmodel"))
  expect_true(.is_softmax_param("thetat", mock_model))
  expect_true(.is_softmax_param("thetant", mock_model))
  expect_false(.is_softmax_param("kappa", mock_model))
})

test_that(".is_softmax_param returns FALSE for non-mixture3p models", {
  mock_model <- structure(list(), class = c("mixture2p", "bmmodel"))
  expect_false(.is_softmax_param("thetat", mock_model))

  mock_sdm <- structure(list(), class = c("sdm", "bmmodel"))
  expect_false(.is_softmax_param("kappa", mock_sdm))
})

# ===========================================================================
# .get_parameter_info()
# ===========================================================================

test_that(".get_parameter_info returns correct info for SDM params", {
  skip_on_cran()
  path <- test_path("assets/bmmfit_example1.rds")
  skip_if_not(file.exists(path), "SDM fixture not available (excluded by .Rbuildignore)")
  fit <- readRDS(path)

  info_c <- .get_parameter_info(fit, "c")
  expect_equal(info_c$type, "dpar")
  expect_equal(info_c$link, "log")
  expect_false(info_c$softmax)

  info_kappa <- .get_parameter_info(fit, "kappa")
  expect_equal(info_kappa$type, "dpar")
  expect_equal(info_kappa$link, "log")
  expect_false(info_kappa$softmax)
})

# ===========================================================================
# native_parameters() and native_transform()
# Tier 1: Unit tests (always run, no fitted model)
# Tier 2: Fixture-based integration tests (skip on CRAN)
# Model-fitting tests live in tests/internal/test-native_parameters.R
# ===========================================================================

load_np_sdm_fit <- function() {
  path <- test_path("assets/bmmfit_example1.rds")
  skip_if_not(file.exists(path), "SDM fixture not available (excluded by .Rbuildignore)")
  readRDS(path)
}

load_np_m3_fit <- function() {
  path <- test_path("assets/bmmfit_m3_ppcheck.rds")
  skip_if_not(file.exists(path), "m3 fixture not available (excluded by .Rbuildignore)")
  readRDS(path)
}

# ===========================================================================
# Tier 1: .np_softmax()
# ===========================================================================

test_that(".np_softmax matches softmax() with an implicit zero reference", {
  set.seed(20)
  mats <- list(a = matrix(rnorm(12), 4), b = matrix(rnorm(12), 4))
  out <- .np_softmax(mats)

  expected <- vapply(
    seq_along(mats$a),
    function(i) softmax(c(mats$a[i], mats$b[i], 0)),
    numeric(3)
  )
  expect_equal(as.vector(out$a), expected[1, ])
  expect_equal(as.vector(out$b), expected[2, ])
  expect_equal(as.vector(out$a + out$b) + expected[3, ], rep(1, 12))
})

test_that(".np_softmax with a single parameter reduces to the logit link", {
  eta <- matrix(c(1.2, -0.3, 0))
  expect_equal(as.vector(.np_softmax(list(p = eta))$p), plogis(as.vector(eta)))
})

test_that(".np_softmax is numerically stable at extreme values", {
  out <- .np_softmax(list(
    a = matrix(c(1000, -1000, 0)),
    b = matrix(c(-1000, 1000, 0))
  ))
  values <- unlist(out)
  expect_true(all(is.finite(values)))
  expect_true(all(values >= 0 & values <= 1))
  expect_equal(as.vector(out$a), c(1, 0, 1 / 3))
})

test_that(".np_softmax preserves names and dimensions", {
  mats <- list(a = matrix(0, 5, 3), b = matrix(0, 5, 3))
  out <- .np_softmax(mats)
  expect_named(out, c("a", "b"))
  expect_equal(lapply(out, dim), lapply(mats, dim))
})

# ===========================================================================
# Tier 1: native_transform()
# ===========================================================================

test_that("native_transform applies the elementwise inverse link of each parameter", {
  model <- .model_sdm()
  set.seed(21)
  linpred <- list(
    mu = matrix(rnorm(6), 3),
    c = matrix(rnorm(6), 3),
    kappa = matrix(rnorm(6), 3)
  )
  out <- native_transform(model, linpred, data.frame())

  expect_named(out, names(linpred))
  expect_equal(lapply(out, dim), lapply(linpred, dim))
  expect_equal(out$c, exp(linpred$c))
  expect_equal(out$kappa, exp(linpred$kappa))
  expect_equal(out$mu, 2 * atan(linpred$mu))
})

test_that("native_transform treats a parameter without a declared link as identity", {
  model <- .model_sdm()
  model$links$kappa <- NULL
  linpred <- list(kappa = matrix(c(1, 2)))
  expect_equal(native_transform(model, linpred, data.frame())$kappa, linpred$kappa)
})

test_that("native_transform softmaxes mixture3p weights without a model specific method", {
  model <- .model_mixture3p()
  set.seed(22)
  linpred <- list(
    kappa = matrix(rnorm(4), 2),
    thetat = matrix(rnorm(4), 2),
    thetant = matrix(rnorm(4), 2)
  )
  out <- native_transform(model, linpred, data.frame())

  expect_equal(out$kappa, exp(linpred$kappa))
  expect_true(all(out$thetat > 0 & out$thetat < 1))
  expect_true(all(out$thetat + out$thetant < 1))
  expect_equal(out[c("thetat", "thetant")], .np_softmax(linpred[c("thetat", "thetant")]))
})

test_that("native_transform is consistent between the logit and softmax mixtures", {
  eta <- matrix(c(0.7, -1.4))
  logit_model <- .model_mixture2p()
  softmax_model <- .model_mixture3p()
  softmax_model$links <- list(thetat = "softmax")

  expect_equal(
    native_transform(logit_model, list(thetat = eta), data.frame())$thetat,
    native_transform(softmax_model, list(thetat = eta), data.frame())$thetat
  )
})

test_that("native_transform errors informatively on an unsupported link", {
  model <- .model_sdm()
  model$links$kappa <- "banana"
  expect_error(
    native_transform(model, list(kappa = matrix(1)), data.frame()),
    "kappa"
  )
  expect_error(
    native_transform(model, list(kappa = matrix(1)), data.frame()),
    "native_transform"
  )
})

test_that("native_transform requires a bmmodel", {
  expect_error(native_transform(list(), list(), data.frame()), "bmmodel")
})

test_that("every link declared by a supported model can be transformed", {
  models <- supported_models(print_call = FALSE)
  links <- unlist(lapply(models, function(name) {
    constructor <- get_model(name)
    versions <- eval(formals(constructor)$version)
    versions <- if (is.null(versions)) list(NULL) else as.list(versions)
    unlist(lapply(versions, function(version) {
      model <- tryCatch(constructor(version = version), error = function(e) NULL)
      unlist(model$links)
    }))
  }))
  links <- unique(links)
  expect_gt(length(links), 0)

  for (link in links[links != "softmax"]) {
    expect_silent(link_transform(0.5, link, inverse = TRUE))
  }
})

# ===========================================================================
# Tier 1: .np_grid_vars() and .np_newdata()
# ===========================================================================

test_that(".np_grid_vars honours the three re_formula regimes", {
  fit <- load_np_m3_fit()
  pars <- names(fit$bmm$model$parameters)

  expect_setequal(.np_grid_vars(fit, pars, NULL), c("cond", "ID"))
  expect_setequal(.np_grid_vars(fit, pars, NA), "cond")
  expect_setequal(.np_grid_vars(fit, pars, ~ (1 | ID)), c("cond", "ID"))
  expect_setequal(.np_grid_vars(fit, pars, ~ (1 | other)), "cond")
})

test_that(".np_grid_vars skips parameters without a formula", {
  fit <- load_np_m3_fit()
  expect_length(.np_grid_vars(fit, "b", NULL), 0)
})

test_that(".np_grid_vars keeps a grouping variable that is also a predictor", {
  fit <- load_np_m3_fit()
  fit$bmm$user_formula$c <- stats::as.formula(c ~ 1 + ID + (1 | ID))
  expect_setequal(.np_grid_vars(fit, "c", NA), "ID")
})

test_that(".np_newdata returns the unique observed cells in a stable order", {
  fit <- load_np_sdm_fit()
  grid <- .np_newdata(fit, "set_size", NULL)
  expect_equal(nrow(grid), 4)
  expect_equal(as.character(grid$set_size), as.character(sort(unique(fit$data$set_size))))
})

test_that(".np_newdata returns a single row when there are no predictors", {
  fit <- load_np_sdm_fit()
  expect_equal(nrow(.np_newdata(fit, character(0), NULL)), 1)
})

test_that(".np_newdata completes user data including matrix columns", {
  fit <- load_np_m3_fit()
  newdata <- data.frame(cond = fit$data$cond[1:2], ID = c(1L, 2L))
  completed <- .np_newdata(fit, c("cond", "ID"), newdata)

  expect_setequal(names(completed), names(fit$data))
  expect_true(is.matrix(completed$Y))
  expect_equal(nrow(completed$Y), 2)
  expect_equal(completed$nTrials, rep(fit$data$nTrials[1], 2))
})

test_that(".np_newdata errors when user data omits a predictor", {
  fit <- load_np_m3_fit()
  expect_error(
    .np_newdata(fit, c("cond", "ID"), data.frame(cond = fit$data$cond[1])),
    "ID"
  )
})

test_that(".np_newdata errors when a grid variable is not in the data", {
  fit <- load_np_sdm_fit()
  expect_error(.np_newdata(fit, "nonexistent", NULL), "nonexistent")
})

# ===========================================================================
# Tier 1: .np_linpred() parameter resolution
# ===========================================================================

test_that(".np_linpred errors when a requested parameter is not in the formula", {
  fit <- load_np_sdm_fit()
  grid <- .np_newdata(fit, "set_size", NULL)
  expect_error(
    .np_linpred(fit, c("c", "ghost"), "ghost", grid, NA, 1:5, list()),
    "ghost"
  )
})

test_that(".np_linpred warns and drops an unresolvable bystander parameter", {
  fit <- load_np_sdm_fit()
  grid <- .np_newdata(fit, "set_size", NULL)
  expect_warning(
    out <- .np_linpred(fit, c("c", "ghost"), "c", grid, NA, 1:5, list()),
    "ghost"
  )
  expect_named(out, "c")
})

# ===========================================================================
# Tier 2: native_parameters() on the SDM fixture
# ===========================================================================

test_that("native_parameters returns one row per draw, cell and parameter", {
  skip_on_cran()
  fit <- load_np_sdm_fit()
  out <- native_parameters(fit)

  expect_named(out, c(".chain", ".iteration", ".draw", "set_size", "parameter", "value"))
  expect_equal(nrow(out), brms::ndraws(fit) * 4 * 3)
  expect_setequal(unique(out$parameter), c("mu", "c", "kappa"))
})

test_that("native_parameters applies the inverse links of the model", {
  skip_on_cran()
  fit <- load_np_sdm_fit()
  native <- native_parameters(fit, draw_ids = 1:20)
  sampling <- native_parameters(fit, draw_ids = 1:20, scale = "sampling")

  expect_equal(
    native$value[native$parameter == "c"],
    exp(sampling$value[sampling$parameter == "c"])
  )
  expect_equal(
    native$value[native$parameter == "mu"],
    2 * atan(sampling$value[sampling$parameter == "mu"])
  )
  expect_true(all(native$value[native$parameter == "kappa"] > 0))
})

test_that("native_parameters returns the untransformed linear predictor on the sampling scale", {
  skip_on_cran()
  fit <- load_np_sdm_fit()
  grid <- fit$data[!duplicated(fit$data$set_size), , drop = FALSE]
  grid <- grid[order(grid$set_size), , drop = FALSE]
  expected <- brms::posterior_linpred(
    fit,
    newdata = grid, dpar = "c", re_formula = NULL, draw_ids = 1:20
  )

  out <- native_parameters(fit, draw_ids = 1:20, scale = "sampling", pars = "c")
  expect_equal(out$value, as.vector(expected))
})

test_that("native_parameters reports fixed parameters at their constant", {
  skip_on_cran()
  fit <- load_np_sdm_fit()
  out <- native_parameters(fit, pars = "mu")
  expect_true(all(out$value == 0))
})

test_that("native_parameters draw indices match the draws of the fit", {
  skip_on_cran()
  fit <- load_np_sdm_fit()
  out <- native_parameters(fit, pars = "kappa")
  reference <- .np_draw_index(fit, seq_len(brms::ndraws(fit)))

  expect_equal(unique(out$.draw), reference$.draw)
  expect_equal(out$.chain[seq_len(nrow(reference))], reference$.chain)
  expect_equal(out$.iteration[seq_len(nrow(reference))], reference$.iteration)
})

test_that("native_parameters uses the same draws for every parameter", {
  skip_on_cran()
  fit <- load_np_sdm_fit()
  out <- native_parameters(fit, ndraws = 10)
  per_par <- tapply(out$.draw, out$parameter, function(ids) sort(unique(ids)))

  expect_equal(length(unique(per_par)), 1)
  expect_length(per_par[[1]], 10)
})

test_that("native_parameters is reproducible given draw_ids", {
  skip_on_cran()
  fit <- load_np_sdm_fit()
  expect_equal(
    native_parameters(fit, draw_ids = c(1, 5, 9)),
    native_parameters(fit, draw_ids = c(1, 5, 9))
  )
})

test_that("native_parameters summarises the transformed draws, not the reverse", {
  skip_on_cran()
  fit <- load_np_sdm_fit()
  grid <- fit$data[!duplicated(fit$data$set_size), , drop = FALSE]
  grid <- grid[order(grid$set_size), , drop = FALSE]
  linpred <- brms::posterior_linpred(fit, newdata = grid, dpar = "c", re_formula = NULL)

  out <- native_parameters(fit, summary = TRUE)
  estimate <- out$Estimate[out$parameter == "c"]

  expect_equal(estimate, unname(colMeans(exp(linpred))))
  expect_true(all(estimate > exp(colMeans(linpred))))
})

test_that("native_parameters summary has the expected shape and interval names", {
  skip_on_cran()
  fit <- load_np_sdm_fit()
  out <- native_parameters(fit, summary = TRUE)
  expect_named(out, c("set_size", "parameter", "Estimate", "Est.Error", "Q2.5", "Q97.5"))
  expect_equal(nrow(out), 4 * 3)

  narrow <- native_parameters(fit, summary = TRUE, prob = 0.5)
  expect_named(narrow, c("set_size", "parameter", "Estimate", "Est.Error", "Q25", "Q75"))
  expect_true(all(narrow$Q75 - narrow$Q25 <= out$Q97.5 - out$Q2.5))
})

test_that("native_parameters summary supports robust estimates", {
  skip_on_cran()
  fit <- load_np_sdm_fit()
  draws <- native_parameters(fit, pars = "c")
  robust <- native_parameters(fit, pars = "c", summary = TRUE, robust = TRUE)

  expect_equal(robust$Estimate, as.vector(tapply(draws$value, draws$set_size, stats::median)))
  expect_equal(robust$Est.Error, as.vector(tapply(draws$value, draws$set_size, stats::mad)))
})

test_that("native_parameters subsets the output without changing the values", {
  skip_on_cran()
  fit <- load_np_sdm_fit()
  all_pars <- native_parameters(fit, draw_ids = 1:20)
  kappa_only <- native_parameters(fit, draw_ids = 1:20, pars = "kappa")

  expect_equal(unique(kappa_only$parameter), "kappa")
  expect_equal(kappa_only, subset(all_pars, parameter == "kappa"), ignore_attr = TRUE)
})

test_that("native_parameters validates its inputs", {
  skip_on_cran()
  fit <- load_np_sdm_fit()
  expect_error(native_parameters(fit, pars = "nope"), "mu.*c.*kappa")
  expect_error(native_parameters(fit, scale = "raw"))
  expect_error(native_parameters(fit$data), "bmmfit")
})

test_that("native_parameters ignores re_formula when the model has no group-level effects", {
  skip_on_cran()
  fit <- load_np_sdm_fit()
  expect_equal(
    native_parameters(fit, re_formula = NA, draw_ids = 1:20),
    native_parameters(fit, re_formula = NULL, draw_ids = 1:20)
  )
})

# ===========================================================================
# Tier 2: native_parameters() on the m3 fixture
# ===========================================================================

test_that("native_parameters handles multinomial models with matrix columns", {
  skip_on_cran()
  fit <- load_np_m3_fit()
  out <- native_parameters(fit, draw_ids = 1:10)

  expect_named(out, c(".chain", ".iteration", ".draw", "cond", "ID", "parameter", "value"))
  expect_equal(nrow(out), 10 * 120 * 4)
  expect_false(any(c("Y", "nTrials", "Idx_corr", "n_corr") %in% names(out)))
})

test_that("native_parameters drops grouping variables when re_formula is NA", {
  skip_on_cran()
  fit <- load_np_m3_fit()
  out <- native_parameters(fit, re_formula = NA, draw_ids = 1:10)

  expect_false("ID" %in% names(out))
  expect_equal(nrow(out), 10 * 3 * 4)
})

test_that("native_parameters transforms m3 nlpars and keeps b at its constant", {
  skip_on_cran()
  fit <- load_np_m3_fit()
  native <- native_parameters(fit, re_formula = NA, draw_ids = 1:10)
  sampling <- native_parameters(fit, re_formula = NA, draw_ids = 1:10, scale = "sampling")

  expect_true(all(native$value[native$parameter == "b"] == 0.1))
  for (par in c("a", "c", "d")) {
    expect_equal(
      native$value[native$parameter == par],
      exp(sampling$value[sampling$parameter == par])
    )
  }
})

test_that("native_parameters completes partial user newdata", {
  skip_on_cran()
  fit <- load_np_m3_fit()
  newdata <- data.frame(
    cond = factor(levels(fit$data$cond), levels = levels(fit$data$cond)),
    ID = 1L
  )
  partial <- native_parameters(fit, newdata = newdata, draw_ids = 1:10, pars = "c")
  full <- native_parameters(fit, draw_ids = 1:10, pars = "c")

  expect_equal(nrow(partial), 10 * 3)
  expect_equal(sort(partial$value), sort(subset(full, ID == 1)$value))
})

test_that("native_parameters marginalises over new grouping levels when asked", {
  skip_on_cran()
  fit <- load_np_m3_fit()
  newdata <- data.frame(
    cond = factor(levels(fit$data$cond), levels = levels(fit$data$cond)),
    ID = 9999L
  )
  expect_error(native_parameters(fit, newdata = newdata, draw_ids = 1:5), "9999")

  out <- native_parameters(
    fit,
    newdata = newdata, draw_ids = 1:5,
    allow_new_levels = TRUE, sample_new_levels = "gaussian"
  )
  expect_true(all(is.finite(out$value)))
})
