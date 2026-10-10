# ===========================================================================
# native_parameters(population_summary = ) and its helpers
# Tier 1: Unit tests of the integrals (no fitted model)
# Tier 2: Fixture-based tests (skip when the fixtures are not available)
# Model-fitting tests live in tests/internal/test-native_parameters.R
# ===========================================================================

load_pop_m3_fit <- function() {
  path <- test_path("assets/bmmfit_m3_ppcheck.rds")
  skip_if_not(file.exists(path), "m3 fixture not available (excluded by .Rbuildignore)")
  restructure(readRDS(path))
}

load_pop_sdm_fit <- function() {
  path <- test_path("assets/bmmfit_example1.rds")
  skip_if_not(file.exists(path), "SDM fixture not available (excluded by .Rbuildignore)")
  restructure(readRDS(path))
}

# dnorm() underflows to zero beyond |z| = 38, where exp() of the integrand can
# still overflow
normal_expectation <- function(f, eta, s) {
  stats::integrate(
    function(z) f(eta + s * z) * stats::dnorm(z), -38, 38,
    rel.tol = 1e-13, abs.tol = 0, subdivisions = 2000L
  )$value
}

population_grid <- function() {
  grid <- expand.grid(eta = c(-3, -1, 0, 0.5, 2), s = c(0.1, 0.5, 1, 2, 3))
  list(eta = matrix(grid$eta, nrow = 1), s = matrix(grid$s, nrow = 1))
}

# ===========================================================================
# Tier 1: quadrature and closed forms
# ===========================================================================

test_that(".np_gh_nodes integrates polynomials of the standard normal exactly", {
  nodes <- .np_gh_nodes(100)
  expect_equal(sum(nodes$w), 1, tolerance = 1e-13)
  expect_equal(sum(nodes$w * nodes$x), 0, tolerance = 1e-13)
  expect_equal(sum(nodes$w * nodes$x^2), 1, tolerance = 1e-12)
  expect_equal(sum(nodes$w * nodes$x^4), 3, tolerance = 1e-11)
})

test_that("the probit population mean matches numerical integration to 1e-8", {
  g <- population_grid()
  out <- .np_population_moments(g$eta, g$s, "probit")
  for (i in seq_along(g$eta)) {
    mean <- normal_expectation(stats::pnorm, g$eta[i], g$s[i])
    second <- normal_expectation(function(v) stats::pnorm(v)^2, g$eta[i], g$s[i])
    expect_lt(abs(out$mean[i] - mean), 1e-8)
    expect_lt(abs(out$sd[i] - sqrt(second - mean^2)), 1e-6)
  }
})

test_that("the logit population mean and SD agree with Monte Carlo within its error", {
  withr::local_seed(531)
  z <- stats::rnorm(1e6)
  g <- population_grid()
  out <- .np_population_moments(g$eta, g$s, "logit")
  for (i in seq_along(g$eta)) {
    p <- stats::plogis(g$eta[i] + g$s[i] * z)
    se_mean <- stats::sd(p) / sqrt(length(p))
    expect_lt(abs(out$mean[i] - mean(p)), 4 * se_mean)
    expect_lt(abs(out$sd[i] - stats::sd(p)), 4 * se_mean + 1e-4)
  }
})

test_that("the logit population mean matches numerical integration at s = 3", {
  out <- .np_population_moments(matrix(c(-2, 0, 1.5), 1), matrix(3, 1, 3), "logit")
  reference <- vapply(c(-2, 0, 1.5), normal_expectation, numeric(1), f = stats::plogis, s = 3)
  expect_lt(max(abs(out$mean - reference)), 1e-8)
})

test_that("the log population mean is exp(eta + s^2 / 2)", {
  g <- population_grid()
  out <- .np_population_moments(g$eta, g$s, "log")
  expect_equal(out$mean, exp(g$eta + g$s^2 / 2))
  for (i in c(1, 7, 13, 25)) {
    mean <- normal_expectation(exp, g$eta[i], g$s[i])
    second <- normal_expectation(function(v) exp(2 * v), g$eta[i], g$s[i])
    expect_equal(out$mean[i], mean, tolerance = 1e-8)
    expect_equal(out$sd[i], sqrt(second - mean^2), tolerance = 1e-7)
  }
})

test_that("the identity link returns the latent mean and SD", {
  g <- population_grid()
  expect_equal(.np_population_moments(g$eta, g$s, "identity"), list(mean = g$eta, sd = g$s))
})

test_that("generic links are integrated by quadrature", {
  out <- .np_population_moments(matrix(0.3, 1), matrix(1.2, 1), "softplus")
  f <- function(v) log1p(exp(v))
  mean <- normal_expectation(f, 0.3, 1.2)
  expect_equal(out$mean[1], mean, tolerance = 1e-10)
  expect_equal(
    out$sd[1],
    sqrt(normal_expectation(function(v) (f(v) - mean)^2, 0.3, 1.2)),
    tolerance = 1e-9
  )
})

test_that("the population mean and median coincide with the median person when s = 0", {
  eta <- matrix(c(-1.2, 0, 0.7, 2), 2)
  s <- matrix(0, 2, 2)
  for (link in c("identity", "log", "logit", "probit", "softplus", "cloglog", "tan_half")) {
    native <- link_transform(eta, link, inverse = TRUE)
    moments <- .np_population_moments(eta, s, link)
    quartiles <- .np_population_quartiles(eta, s, link)
    expect_equal(moments$mean, native, info = link)
    expect_equal(moments$sd, matrix(0, 2, 2), info = link)
    expect_equal(quartiles$median, native, info = link)
    if (link != "tan_half") {
      expect_equal(quartiles$q25, native, info = link)
      expect_equal(quartiles$q75, native, info = link)
    }
  }
})

test_that("the median and quartiles are the inverse link of the latent quartiles", {
  withr::local_seed(531)
  draws <- stats::plogis(0.4 + 1.7 * stats::rnorm(2e5))
  out <- .np_population_quartiles(matrix(0.4, 1), matrix(1.7, 1), "logit")
  expected <- stats::quantile(draws, c(0.5, 0.25, 0.75), names = FALSE)
  expect_equal(c(out$median, out$q25, out$q75), expected, tolerance = 0.01)
})

test_that("the quartiles stay ordered under a decreasing inverse link", {
  out <- .np_population_quartiles(matrix(c(-1, 0.5), 1), matrix(c(0.8, 2), 1), "loglog")
  expect_true(all(out$q25 < out$median & out$median < out$q75))
})

test_that("the inverse link has no population mean", {
  expect_warning(
    out <- .np_population_moments(matrix(1, 1), matrix(1, 1), "inverse"),
    "does not exist"
  )
  expect_true(is.na(out$mean[1]))
})

test_that("the population statistics keep the dimensions of the draws", {
  eta <- matrix(stats::rnorm(6), 2, 3)
  s <- matrix(abs(stats::rnorm(6)), 2, 3)
  for (summary in c("mean", "median")) {
    for (link in c("log", "logit", "probit", "tan_half")) {
      out <- .np_population_normal(eta, s, link, summary)
      expect_true(all(vapply(out, function(m) identical(dim(m), c(2L, 3L)), logical(1))))
    }
  }
})

# ===========================================================================
# Tier 1: circular location parameters
# ===========================================================================

test_that("the circular mean and SD match numerical integration", {
  eta <- matrix(c(-0.5, 0.3, 1), 1)
  s <- matrix(c(0.4, 1, 1.5), 1)
  out <- .np_circular_moments(eta, s)
  for (i in seq_along(eta)) {
    cosine <- normal_expectation(function(v) cos(2 * atan(v)), eta[i], s[i])
    sine <- normal_expectation(function(v) sin(2 * atan(v)), eta[i], s[i])
    # 200 nodes: about 1e-12 at s = 1, 1e-5 at s = 2
    expect_equal(out$mean[i], atan2(sine, cosine), tolerance = 1e-6)
    expect_equal(out$sd[i], sqrt(-2 * log(sqrt(cosine^2 + sine^2))), tolerance = 1e-6)
  }
})

test_that("the circular median minimises the expected arc distance", {
  # the arc distance has a kink that integrate() cannot handle; a stratified
  # grid of normal quantiles approximates E d to about 1e-6
  z <- stats::qnorm((seq_len(2e4) - 0.5) / 2e4)
  cases <- list(c(0.3, 0.4), c(1, 1), c(-0.5, 2), c(2, 0.7))
  for (case in cases) {
    theta <- 2 * atan(case[1] + case[2] * z)
    expected_distance <- function(m) {
      d <- abs(theta - m)
      mean(pmin(d, 2 * pi - d))
    }
    candidates <- seq(-pi, pi, length.out = 721)
    coarse <- candidates[which.min(vapply(candidates, expected_distance, numeric(1)))]
    best <- stats::optimize(expected_distance, coarse + c(-0.02, 0.02), tol = 1e-8)$minimum
    median <- .np_circular_median(case[1], case[2])
    gap <- abs(atan2(sin(median - best), cos(median - best)))
    expect_lt(gap, 1e-3, label = paste(case, collapse = ", "))
  }
})

test_that("the circular median differs from 2 * atan(eta) once mass wraps around", {
  # eta = -0.5, s = 2 puts most mass beyond -pi/2
  expect_gt(abs(.np_circular_median(-0.5, 2) - 2 * atan(-0.5)), 1)
})

test_that("the circular median of a sample matches a brute-force search", {
  withr::local_seed(531)
  for (i in 1:5) {
    theta <- 2 * atan(stats::rnorm(1, 0, 1) + 1.5 * stats::rnorm(101))
    arc <- outer(theta, theta, function(a, b) {
      d <- abs(a - b)
      pmin(d, 2 * pi - d)
    })
    expect_equal(.np_circular_sample_median(theta), theta[which.min(colSums(arc))])
  }
})

test_that("the analytic and sample circular medians agree", {
  withr::local_seed(531)
  theta <- 2 * atan(-0.5 + 2 * stats::rnorm(2e4))
  sample_median <- .np_circular_sample_median(theta)
  analytic <- .np_circular_median(-0.5, 2)
  expect_lt(abs(atan2(sin(sample_median - analytic), cos(sample_median - analytic))), 0.05)
})

# ===========================================================================
# Tier 1: latent variance and Monte Carlo statistics
# ===========================================================================

test_that(".np_quadratic_form computes z' Sigma z per draw and cell", {
  z <- rbind(c(1, 0, 2), c(1, 1, -1))
  sds <- rbind(c(0.5, 1, 0.3), c(0.2, 0.4, 0.6))
  rho <- c(0.3, -0.5)
  out <- .np_quadratic_form(z, sds, list(list(j = 1, k = 2, rho = rho)))
  for (d in 1:2) {
    correlation <- diag(3)
    correlation[1, 2] <- correlation[2, 1] <- rho[d]
    sigma <- diag(sds[d, ]) %*% correlation %*% diag(sds[d, ])
    for (cell in 1:2) {
      expect_equal(out[d, cell], drop(t(z[cell, ]) %*% sigma %*% z[cell, ]))
    }
  }
})

test_that(".np_mc_stats summarises the levels of each cell separately", {
  values <- cbind(matrix(1:6, 2), matrix(10 * (1:6), 2))
  out <- .np_mc_stats(values, n_cells = 2, n_levels = 3, "mean", circular = FALSE)
  expect_equal(out$mean[, 1], rowMeans(values[, 1:3]))
  expect_equal(out$mean[, 2], rowMeans(values[, 4:6]))
  expect_equal(out$sd[, 2], apply(values[, 4:6], 1, stats::sd))

  out <- .np_mc_stats(values, 2, 3, "median", circular = FALSE)
  expect_named(out, c("median", "q25", "q75"))
  expect_equal(out$median[, 2], apply(values[, 4:6], 1, stats::median))
})

test_that(".np_check_mc_error warns when the Monte Carlo error dominates", {
  withr::local_seed(531)
  # posterior SD of the mean 0.01, within-draw SD 1: with 100 levels the Monte
  # Carlo SE is 0.1
  noisy <- list(par = list(
    mean = matrix(stats::rnorm(400, 0.5, sqrt(0.01^2 + 1 / 100)), 400),
    sd = matrix(1, 400)
  ))
  expect_warning(.np_check_mc_error(noisy, 100), "ndraws_population")
  precise <- list(par = list(
    mean = matrix(stats::rnorm(400, 0.5, 1), 400),
    sd = matrix(1, 400)
  ))
  expect_no_warning(.np_check_mc_error(precise, 1000))
})

test_that("parameters fixed to a constant in the predictions have no group-level terms", {
  prep <- list(
    dpars = list(mu2 = 0, kappa = list(re = list(Z = list(ID = matrix(1)))), theta1 = list(re = list())),
    nlpars = list(c = list(re = list(Z = list(item = matrix(1)))))
  )
  expect_null(.np_prep_term(prep, "mu2"))
  expect_equal(.np_prep_term(prep, "c"), prep$nlpars$c)
  expect_setequal(.np_prep_groups(prep), c("ID", "item"))
})

test_that(".np_reserved_names includes statistic only for population summaries", {
  expect_false("statistic" %in% .np_reserved_names(FALSE, 0.95))
  expect_true("statistic" %in% .np_reserved_names(FALSE, 0.95, population = TRUE))
  expect_true("statistic" %in% .np_reserved_names(TRUE, 0.95, population = TRUE))
})

# ===========================================================================
# Tier 2: native_parameters(population_summary = ) on the fixtures
# ===========================================================================

test_that("the population mean of a log-link parameter is exp(eta + s^2 / 2) on a fit", {
  skip_on_cran()
  fit <- load_pop_m3_fit()
  median_person <- native_parameters(fit, re_formula = NA, scale = "sampling", pars = "a")
  population <- native_parameters(fit, population_summary = "mean", pars = "a")
  draws <- posterior::as_draws_df(fit)[seq_len(brms::ndraws(fit)), ]
  # (1 + cond || ID): the intercept and the condition's own effect add up
  slope <- c(
    "new distractors" = 0,
    "old reordered" = 1,
    "old same" = 1
  )
  slope_sd <- list(
    "new distractors" = 0,
    "old reordered" = draws$sd_ID__a_condoldreordered,
    "old same" = draws$sd_ID__a_condoldsame
  )
  for (cell in levels(median_person$cond)) {
    eta <- median_person$value[median_person$cond == cell]
    s2 <- draws$sd_ID__a_Intercept^2 + slope[[cell]] * slope_sd[[cell]]^2
    rows <- population[population$cond == cell, ]
    expect_equal(rows$value[rows$statistic == "mean"], exp(eta + s2 / 2))
    expect_equal(rows$value[rows$statistic == "sd"], exp(eta + s2 / 2) * sqrt(expm1(s2)))
  }
})

test_that("the population median of a log-link parameter is the median person", {
  skip_on_cran()
  fit <- load_pop_m3_fit()
  median_person <- native_parameters(fit, re_formula = NA, pars = "c")
  population <- native_parameters(fit, population_summary = "median", pars = "c")
  expect_equal(population$value[population$statistic == "median"], median_person$value)
  expect_true(all(
    population$value[population$statistic == "q25"] < median_person$value &
      population$value[population$statistic == "q75"] > median_person$value
  ))
})

test_that("Monte Carlo over new levels agrees with the closed form on a fit", {
  skip_on_cran()
  fit <- load_pop_m3_fit()
  grid_vars <- .np_grid_vars(fit, names(fit$bmm$model$parameters), NA)
  newdata <- .np_newdata(fit, grid_vars, NULL)
  draw_ids <- 1:20
  links <- list(a = "log", b = "identity", c = "log", d = "log")
  withr::local_seed(531)
  mc <- .np_population_mc(
    fit, c("b", "a", "c", "d"), "a", newdata, grid_vars, NULL, draw_ids,
    "mean", native = TRUE, links = links, ndraws_population = 20000, dots = list()
  )
  prep <- .np_re_design(fit, newdata, NULL, 1)
  latent <- .np_latent_sd(fit, prep, "a", draw_ids)
  eta <- .np_linpred(fit, "a", "a", newdata, NA, draw_ids, list())$a
  exact <- .np_population_moments(eta, latent$sd, "log")
  se <- exact$sd / sqrt(20000)
  expect_lt(max(abs(mc$a$mean - exact$mean) / se), 5)
})

test_that("native_parameters returns a statistic column and the method", {
  skip_on_cran()
  fit <- load_pop_m3_fit()
  out <- native_parameters(fit, population_summary = "mean", draw_ids = 1:10)
  expect_named(
    out,
    c(".chain", ".iteration", ".draw", "cond", "parameter", "statistic", "value")
  )
  expect_setequal(unique(out$statistic), c("mean", "sd"))
  expect_false("ID" %in% names(out))
  method <- attr(out, "population_method")
  expect_named(method, c("parameter", "statistic", "link", "method", "n", "groups"))
  expect_equal(method$method[method$parameter == "a"], c("closed form", "closed form"))
  expect_equal(method$method[method$parameter == "b"], rep("no group-level effects", 2))

  summary <- native_parameters(fit, population_summary = "median", draw_ids = 1:10, summary = TRUE)
  expect_named(
    summary,
    c("cond", "parameter", "statistic", "Estimate", "Est.Error", "Q2.5", "Q97.5")
  )
  expect_setequal(unique(summary$statistic), c("median", "q25", "q75"))
  expect_false(is.null(attr(summary, "population_method")))
})

test_that("population summaries without group-level variance reduce to the median person", {
  skip_on_cran()
  fit <- load_pop_m3_fit()
  expect_message(
    out <- native_parameters(fit, re_formula = NA, population_summary = "mean", draw_ids = 1:10),
    "no group-level effects"
  )
  median_person <- native_parameters(fit, re_formula = NA, draw_ids = 1:10)
  expect_equal(out$value[out$statistic == "mean"], median_person$value)
  expect_true(all(out$value[out$statistic == "sd"] == 0))

  sdm <- load_pop_sdm_fit()
  expect_message(
    out <- native_parameters(sdm, population_summary = "median", pars = "c"),
    "no group-level effects"
  )
  expect_equal(
    out$value[out$statistic == "q25"],
    native_parameters(sdm, pars = "c")$value
  )
})

test_that("a re_formula restricts the terms integrated over", {
  skip_on_cran()
  fit <- load_pop_m3_fit()
  draws <- posterior::as_draws_df(fit)[1:10, ]
  intercept_only <- native_parameters(
    fit, population_summary = "mean", re_formula = ~ (1 | ID), pars = "a", draw_ids = 1:10
  )
  median_person <- native_parameters(fit, re_formula = NA, scale = "sampling", pars = "a", draw_ids = 1:10)
  rows <- intercept_only$statistic == "mean"
  expect_equal(
    intercept_only$value[rows],
    exp(median_person$value + rep(draws$sd_ID__a_Intercept^2 / 2, 3))
  )
})

test_that("population summaries refuse arguments that sample new levels", {
  skip_on_cran()
  fit <- load_pop_m3_fit()
  expect_error(
    native_parameters(fit, population_summary = "mean", allow_new_levels = TRUE),
    "allow_new_levels"
  )
  expect_error(
    native_parameters(fit, population_summary = "mean", ndraws_population = 1),
    "ndraws_population"
  )
  expect_error(native_parameters(fit, population_summary = "average"))
})

test_that("the sampling scale reports the latent normal distribution", {
  skip_on_cran()
  fit <- load_pop_m3_fit()
  latent <- native_parameters(
    fit, population_summary = "median", scale = "sampling", pars = "a", draw_ids = 1:10
  )
  eta <- native_parameters(fit, re_formula = NA, scale = "sampling", pars = "a", draw_ids = 1:10)
  expect_equal(latent$value[latent$statistic == "median"], eta$value)
  native <- native_parameters(fit, population_summary = "median", pars = "a", draw_ids = 1:10)
  expect_equal(
    native$value[native$statistic == "q75"],
    exp(latent$value[latent$statistic == "q75"])
  )
})
