# Tier 1: Unit tests for .bmmfit_resolve_par() using mock fixtures
# These test the parameter routing logic without requiring posterior samples

load_mock_mixture2p <- function() {
  readRDS(system.file("extdata", "mock_bmmfit_mixture2p.rds", package = "bmm"))
}

load_mock_m3 <- function() {
  readRDS(system.file("extdata", "mock_bmmfit_m3.rds", package = "bmm"))
}

load_sdm_fit <- function() {
  readRDS(test_path("assets/bmmfit_example1.rds"))
}

test_that(".bmmfit_resolve_par routes dpar to nlpar for mixture2p", {
  mock <- load_mock_mixture2p()
  result <- .bmmfit_resolve_par(mock, dpar = "kappa", nlpar = NULL)
  expect_null(result$dpar)
  expect_equal(result$nlpar, "kappa")
})

test_that(".bmmfit_resolve_par routes dpar to nlpar for mixture2p thetat", {
  mock <- load_mock_mixture2p()
  result <- .bmmfit_resolve_par(mock, dpar = "thetat", nlpar = NULL)
  expect_null(result$dpar)
  expect_equal(result$nlpar, "thetat")
})

test_that(".bmmfit_resolve_par routes dpar to nlpar for m3", {
  mock <- load_mock_m3()
  result <- .bmmfit_resolve_par(mock, dpar = "c", nlpar = NULL)
  expect_null(result$dpar)
  expect_equal(result$nlpar, "c")
})

test_that(".bmmfit_resolve_par leaves dpar unchanged for SDM", {
  fit <- load_sdm_fit()
  result <- .bmmfit_resolve_par(fit, dpar = "c", nlpar = NULL)
  expect_equal(result$dpar, "c")
  expect_null(result$nlpar)
})

test_that(".bmmfit_resolve_par leaves dpar unchanged for SDM kappa", {
  fit <- load_sdm_fit()
  result <- .bmmfit_resolve_par(fit, dpar = "kappa", nlpar = NULL)
  expect_equal(result$dpar, "kappa")
  expect_null(result$nlpar)
})

test_that(".bmmfit_resolve_par respects explicit nlpar argument", {
  mock <- load_mock_mixture2p()
  result <- .bmmfit_resolve_par(mock, dpar = NULL, nlpar = "kappa")
  expect_null(result$dpar)
  expect_equal(result$nlpar, "kappa")
})

test_that(".bmmfit_resolve_par does not reroute when both dpar and nlpar given", {
  mock <- load_mock_mixture2p()
  result <- .bmmfit_resolve_par(mock, dpar = "kappa", nlpar = "thetat")
  expect_equal(result$dpar, "kappa")
  expect_equal(result$nlpar, "thetat")
})

test_that(".bmmfit_resolve_par passes through when both are NULL", {
  mock <- load_mock_mixture2p()
  result <- .bmmfit_resolve_par(mock, dpar = NULL, nlpar = NULL)
  expect_null(result$dpar)
  expect_null(result$nlpar)
})

test_that(".bmmfit_resolve_par errors on invalid dpar name", {
  mock <- load_mock_mixture2p()
  expect_error(
    .bmmfit_resolve_par(mock, dpar = "not_a_par", nlpar = NULL),
    "not found in model"
  )
})


# Tier 2: Fixture-based integration tests (SDM fixture — has posterior samples)

test_that("emmeans() returns valid emmGrid for SDM dpar", {
  skip_if_not_installed("emmeans")
  fit <- load_sdm_fit()
  em <- emmeans::emmeans(fit, ~ set_size, dpar = "c")
  expect_s4_class(em, "emmGrid")
  expect_true(nrow(as.data.frame(em)) > 0)
})

test_that("pairs() produces valid contrast from emmGrid", {
  skip_if_not_installed("emmeans")
  fit <- load_sdm_fit()
  em <- emmeans::emmeans(fit, ~ set_size, dpar = "c")
  contrasts <- pairs(em)
  expect_s4_class(contrasts, "emmGrid")
  expect_true(nrow(as.data.frame(contrasts)) > 0)
})
