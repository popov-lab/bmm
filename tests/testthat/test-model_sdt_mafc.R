# Test m-AFC SDT model specification and integration

############################################################################# !
# MODEL CONSTRUCTOR TESTS                                                ####
############################################################################# !

test_that("sdt_mafc model can be created with default arguments", {
  expect_silent(sdt_mafc("n_correct", "n_trials", m = 4))
})

test_that("sdt_mafc model can be created with all distribution options", {
  expect_silent(sdt_mafc("n_correct", "n_trials", m = 4, dist = "normal"))
  expect_silent(sdt_mafc("n_correct", "n_trials", m = 4, dist = "logistic"))
  expect_silent(sdt_mafc("n_correct", "n_trials", m = 4, dist = "gumbel_min"))
  expect_silent(sdt_mafc("n_correct", "n_trials", m = 4, dist = "gumbel_max"))
})

test_that("sdt_mafc model has correct class structure", {
  model <- sdt_mafc("n_correct", "n_trials", m = 4)
  expect_s3_class(model, "bmmodel")
  expect_s3_class(model, "sdt")
  expect_s3_class(model, "sdt_mafc")
})

test_that("sdt_mafc model parameters are correctly defined", {
  model <- sdt_mafc("n_correct", "n_trials", m = 4)
  expect_true("dprime" %in% names(model$parameters))
  expect_false("criterion" %in% names(model$parameters))
})

test_that("sdt_mafc model has correct default link functions", {
  model <- sdt_mafc("n_correct", "n_trials", m = 4)
  expect_equal(model$links$dprime, "identity")
})

test_that("sdt_mafc model accepts custom links", {
  model <- sdt_mafc("n_correct", "n_trials", m = 4, links = list(dprime = "log"))
  expect_equal(model$links$dprime, "log")
})

test_that("sdt_mafc model stores m and distribution info correctly", {
  model <- sdt_mafc("n_correct", "n_trials", m = 6, dist = "gumbel_min")
  expect_equal(model$other_vars$m, 6L)
  expect_equal(model$other_vars$dist, "gumbel_min")

  model2 <- sdt_mafc("n_correct", "n_trials", m = 3, dist = "logistic")
  expect_equal(model2$other_vars$dist, "logistic")
})

test_that("sdt_mafc accepts m as a constant or a column name", {
  m_const <- sdt_mafc("n_correct", "n_trials", m = 4)
  expect_equal(m_const$other_vars$m, 4L)

  m_col <- sdt_mafc("n_correct", "n_trials", m = "set_size")
  expect_identical(m_col$other_vars$m, "set_size")
})

test_that("sdt_mafc model has default priors and init_ranges", {
  model <- sdt_mafc("n_correct", "n_trials", m = 4)
  expect_true("dprime" %in% names(model$default_priors))
  expect_length(model$init_ranges$dprime, 2)
  expect_true(model$init_ranges$dprime[1] < model$init_ranges$dprime[2])
})

test_that("sdt_mafc requires response, n_trials, and m", {
  expect_error(sdt_mafc("n_correct", "n_trials"))
  expect_error(sdt_mafc("n_correct", m = 4))
})

test_that("sdt_mafc rejects m < 2 and invalid distributions", {
  expect_error(sdt_mafc("n_correct", "n_trials", m = 1), "m must be")
  expect_error(sdt_mafc("n_correct", "n_trials", m = 4, dist = "foo"))
})


############################################################################# !
# CHECK_DATA TESTS                                                       ####
############################################################################# !

test_that("sdt_mafc check_data adds m_afc and dist_type columns", {
  model <- sdt_mafc("n_correct", "n_trials", m = 4, dist = "logistic")
  formula <- bmf(dprime ~ 1)
  dat <- data.frame(n_correct = c(70, 65), n_trials = c(100, 100))

  result <- check_data(model, dat, formula)
  expect_true(all(c("m_afc", "dist_type") %in% colnames(result)))
  expect_equal(unique(result$m_afc), 4L)
  expect_equal(unique(result$dist_type), 4L)
})

test_that("sdt_mafc check_data resolves per-row set size from a data column", {
  model <- sdt_mafc("n_correct", "n_trials", m = "set_size")
  dat <- data.frame(
    set_size = c(2L, 4L, 6L),
    n_correct = c(70, 55, 40),
    n_trials = c(100, 100, 100)
  )

  result <- check_data(model, dat, bmf(dprime ~ 1))
  expect_equal(result$m_afc, c(2L, 4L, 6L))
})

test_that("sdt_mafc check_data errors when the set-size column is missing", {
  model <- sdt_mafc("n_correct", "n_trials", m = "set_size")
  dat <- data.frame(n_correct = c(70, 65), n_trials = c(100, 100))
  expect_error(check_data(model, dat, bmf(dprime ~ 1)), "Set-size column")
})

test_that("sdt_mafc check_data validates response counts", {
  model <- sdt_mafc("n_correct", "n_trials", m = 4)
  formula <- bmf(dprime ~ 1)

  expect_error(
    check_data(model, data.frame(n_correct = c(-1, 40), n_trials = c(50, 50)), formula),
    "non-negative"
  )
  expect_error(
    check_data(model, data.frame(n_correct = c(60, 40), n_trials = c(50, 50)), formula),
    "must not exceed"
  )
})


############################################################################# !
# DISTRIBUTION FUNCTION TESTS                                            ####
############################################################################# !

test_that("rsdt_mafc generates counts matching the design", {
  n_correct <- rsdt_mafc(5, 50, m = 4, dprime = 1.5)
  expect_length(n_correct, 5)
  expect_true(all(n_correct >= 0))
  expect_true(all(n_correct <= 50))
})

test_that("rsdt_mafc recycles vectorized parameters per observation", {
  n_correct <- rsdt_mafc(12, 100, m = 4, dprime = rnorm(12, 1.5, 0.3))
  expect_length(n_correct, 12)

  n2 <- rsdt_mafc(3, 100, m = c(2, 4, 8), dprime = 1.5)
  expect_length(n2, 3)
})

test_that("rsdt_mafc validates input", {
  expect_error(rsdt_mafc(c(2, 3), 100, m = 4, dprime = 1),
               "single positive integer")
  expect_error(rsdt_mafc(2, 100, m = 1, dprime = 1), "m must be")
})

test_that("dsdt_mafc recycles parameters across observations", {
  d <- dsdt_mafc(n_correct = c(60, 80), n_trials = 100,
                 dprime = c(1, 2), m = c(4, 4))
  expect_length(d, 2)
  expect_false(d[1] == d[2])
})

test_that("dsdt_mafc returns a valid binomial density and respects log", {
  d <- dsdt_mafc(n_correct = 80, n_trials = 100, dprime = 1.5, m = 4)
  expect_true(d > 0 && d <= 1)
  ld <- dsdt_mafc(n_correct = 80, n_trials = 100, dprime = 1.5, m = 4, log = TRUE)
  expect_equal(log(d), ld, tolerance = 1e-10)
})

test_that("dsdt_mafc is vectorized over observations and works for all dists", {
  d <- dsdt_mafc(n_correct = c(60, 80), n_trials = c(100, 100), dprime = 1.5, m = 4)
  expect_length(d, 2)
  for (di in c("normal", "logistic", "gumbel_min", "gumbel_max")) {
    expect_true(dsdt_mafc(70, 100, dprime = 1.5, m = 4, dist = di) > 0,
                info = paste("dist:", di))
  }
})

test_that("dsdt_mafc validates input", {
  expect_error(dsdt_mafc(80, 100, dprime = 1, m = 1), "m must be")
  expect_error(dsdt_mafc(120, 100, dprime = 1, m = 4), "must not exceed")
})

test_that("m-AFC probability correct equals chance (1/m) at d' = 0", {
  for (di in c("normal", "logistic", "gumbel_min", "gumbel_max")) {
    for (m in c(2L, 3L, 4L, 8L)) {
      expect_equal(.mafc_pc_r(0, m, di), 1 / m, tolerance = 1e-5,
                   info = paste(di, "m =", m))
    }
  }
})

test_that("m-AFC probability correct is monotone increasing in d'", {
  for (di in c("normal", "logistic", "gumbel_min", "gumbel_max")) {
    pcs <- vapply(c(0, 0.5, 1, 2, 3), .mafc_pc_r, numeric(1), m = 4L, dist = di)
    expect_true(all(diff(pcs) > 0), info = di)
    expect_true(all(pcs > 0 & pcs < 1), info = di)
  }
})

test_that("gumbel_max m-AFC matches the exact softmax closed form", {
  for (m in c(2L, 3L, 4L, 8L)) {
    for (d in c(0, 1, 2, 3)) {
      expect_equal(.mafc_pc_r(d, m, "gumbel_max"), 1 / (1 + (m - 1) * exp(-d)),
                   tolerance = 1e-12, info = paste("m =", m, "d =", d))
    }
  }
})

test_that("gumbel_min m-AFC matches the exact Gamma-ratio closed form", {
  for (m in c(2L, 3L, 4L, 8L)) {
    for (d in c(0, 1, 2, 3)) {
      expect_equal(
        .mafc_pc_r(d, m, "gumbel_min"),
        exp(lgamma(1 + exp(-d)) + lgamma(m) - lgamma(m + exp(-d))),
        tolerance = 1e-12, info = paste("m =", m, "d =", d)
      )
    }
  }
})

test_that("normal m-AFC matches Phi(d'/sqrt(2)) at m = 2", {
  for (d in c(0, 0.5, 1, 2)) {
    expect_equal(.mafc_pc_r(d, 2L, "normal"), pnorm(d / sqrt(2)),
                 tolerance = 1e-10, info = paste("d =", d))
  }
})

test_that(".mafc_pc_r vectorized path matches elementwise evaluation", {
  d <- c(0.3, 1.1, 2.2)
  mm <- c(2L, 4L, 6L)
  for (di in c("normal", "logistic", "gumbel_min", "gumbel_max")) {
    vec <- .mafc_pc_r(d, mm, di)
    ref <- vapply(seq_along(d),
                  function(i) .mafc_pc_r(d[i], mm[i], di), numeric(1))
    expect_equal(vec, ref, tolerance = 1e-12, info = di)
  }
})

test_that("m-AFC probability correct matches independent numerical oracles", {
  # Independent oracles: latent-space integration / binomial-sum closed form,
  # none of which reuse the implementation in .mafc_pc_r.
  o_norm <- function(d, m) {
    stats::integrate(function(z) dnorm(z) * pnorm(z + d)^(m - 1),
                     -Inf, Inf, rel.tol = 1e-12)$value
  }
  o_logis <- function(d, m) {
    stats::integrate(function(z) dlogis(z) * plogis(z + d)^(m - 1),
                     -Inf, Inf, rel.tol = 1e-12)$value
  }
  # largest extreme value: pdf exp(-z - exp(-z)), cdf exp(-exp(-z))
  o_gmax <- function(d, m) {
    stats::integrate(function(z) exp(-z - exp(-z)) * exp(-exp(-(z + d)))^(m - 1),
                     -Inf, Inf, rel.tol = 1e-12)$value
  }
  # smallest extreme value: binomial-sum closed form of the same integral
  o_gmin <- function(d, m) {
    a <- exp(d)
    k <- 0:(m - 1)
    sum(choose(m - 1, k) * (-1)^k / (a * k + 1))
  }
  for (m in c(3L, 4L, 6L)) {
    for (d in c(0.5, 1, 2)) {
      expect_equal(.mafc_pc_r(d, m, "normal"), o_norm(d, m), tolerance = 1e-5)
      expect_equal(.mafc_pc_r(d, m, "logistic"), o_logis(d, m), tolerance = 1e-7)
      expect_equal(.mafc_pc_r(d, m, "gumbel_max"), o_gmax(d, m), tolerance = 1e-6)
      expect_equal(.mafc_pc_r(d, m, "gumbel_min"), o_gmin(d, m), tolerance = 1e-9)
    }
  }
})


############################################################################# !
# FORMULA AND CONFIGURE_MODEL TESTS                                      ####
############################################################################# !

test_that("sdt_mafc produces valid stancode including its custom functions", {
  dat <- data.frame(n_correct = c(70, 65, 80, 60), n_trials = rep(100, 4),
                    id = c(1, 1, 2, 2))
  model <- sdt_mafc("n_correct", "n_trials", m = 4)
  code <- stancode(bmf(dprime ~ 1), data = dat, model = model)
  expect_true(grepl("sdt_mafc", code))
  expect_true(grepl("mafc_pc", code))
  expect_true(grepl("sdt_quantile", code))
})

test_that("sdt_mafc produces valid stancode for all distributions", {
  dat <- data.frame(n_correct = c(70, 65), n_trials = c(100, 100))
  for (d in c("normal", "logistic", "gumbel_min", "gumbel_max")) {
    model <- sdt_mafc("n_correct", "n_trials", m = 4, dist = d)
    code <- stancode(bmf(dprime ~ 1), data = dat, model = model)
    expect_true(nchar(code) > 0, info = paste("dist:", d))
  }
})

test_that("sdt_mafc handles predictors and random effects in the formula", {
  dat <- data.frame(n_correct = c(70, 65, 80, 60), n_trials = rep(100, 4),
                    condition = c("A", "B", "A", "B"), id = c(1, 1, 2, 2))
  model <- sdt_mafc("n_correct", "n_trials", m = 4)
  expect_true(nchar(stancode(bmf(dprime ~ condition), dat, model = model)) > 0)
  expect_true(nchar(stancode(bmf(dprime ~ 1 + (1 | id)), dat, model = model)) > 0)
})

test_that("sdt_mafc integrates with the bmm pipeline via mock backend", {
  dat <- data.frame(n_trials = rep(100L, 8))
  dat$n_correct <- rsdt_mafc(nrow(dat), dat$n_trials, m = 4, dprime = 1.2)
  for (d in c("normal", "logistic", "gumbel_min", "gumbel_max")) {
    model <- sdt_mafc("n_correct", "n_trials", m = 4, dist = d)
    expect_silent(
      bmm(bmf(dprime ~ 1), dat, model, backend = "mock", mock_fit = 1, rename = FALSE)
    )
  }
})

test_that("sdt_mafc fits mixed set sizes through the pipeline (mock)", {
  dat <- data.frame(set_size = rep(c(2L, 4L, 6L), each = 8), n_trials = 100L)
  dat$n_correct <- rsdt_mafc(nrow(dat), dat$n_trials, m = dat$set_size,
                             dprime = 1.2)
  model <- sdt_mafc("n_correct", "n_trials", m = "set_size")
  expect_silent(
    bmm(bmf(dprime ~ 1), dat, model, backend = "mock", mock_fit = 1, rename = FALSE)
  )
})

test_that("sdt_mafc default_prior returns a valid prior object", {
  dat <- data.frame(n_correct = c(70, 65), n_trials = c(100, 100))
  model <- sdt_mafc("n_correct", "n_trials", m = 4)
  prior <- default_prior(bmf(dprime ~ 1), data = dat, model = model)
  expect_s3_class(prior, "brmsprior")
})
