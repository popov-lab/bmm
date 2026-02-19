# Tier 1: Unit tests for .bmmfit_resolve_par() using mock fixtures
# These test the parameter routing logic without requiring posterior samples

test_that(".bmmfit_resolve_par routes dpar to nlpar for mixture2p", {
  mock <- readRDS(test_path("assets/mock_bmmfit_mixture2p.rds"))
  result <- .bmmfit_resolve_par(mock, dpar = "kappa", nlpar = NULL)
  expect_null(result$dpar)
  expect_equal(result$nlpar, "kappa")
})

test_that(".bmmfit_resolve_par routes dpar to nlpar for mixture2p thetat", {
  mock <- readRDS(test_path("assets/mock_bmmfit_mixture2p.rds"))
  result <- .bmmfit_resolve_par(mock, dpar = "thetat", nlpar = NULL)
  expect_null(result$dpar)
  expect_equal(result$nlpar, "thetat")
})

test_that(".bmmfit_resolve_par routes dpar to nlpar for m3", {
  mock <- readRDS(test_path("assets/mock_bmmfit_m3.rds"))
  result <- .bmmfit_resolve_par(mock, dpar = "c", nlpar = NULL)
  expect_null(result$dpar)
  expect_equal(result$nlpar, "c")
})

test_that(".bmmfit_resolve_par leaves dpar unchanged for SDM", {
  fit <- readRDS(test_path("assets/bmmfit_example1.rds"))
  result <- .bmmfit_resolve_par(fit, dpar = "c", nlpar = NULL)
  expect_equal(result$dpar, "c")
  expect_null(result$nlpar)
})

test_that(".bmmfit_resolve_par leaves dpar unchanged for SDM kappa", {
  fit <- readRDS(test_path("assets/bmmfit_example1.rds"))
  result <- .bmmfit_resolve_par(fit, dpar = "kappa", nlpar = NULL)
  expect_equal(result$dpar, "kappa")
  expect_null(result$nlpar)
})

test_that(".bmmfit_resolve_par respects explicit nlpar argument", {
  mock <- readRDS(test_path("assets/mock_bmmfit_mixture2p.rds"))
  result <- .bmmfit_resolve_par(mock, dpar = NULL, nlpar = "kappa")
  expect_null(result$dpar)
  expect_equal(result$nlpar, "kappa")
})

test_that(".bmmfit_resolve_par does not reroute when both dpar and nlpar given", {
  mock <- readRDS(test_path("assets/mock_bmmfit_mixture2p.rds"))
  result <- .bmmfit_resolve_par(mock, dpar = "kappa", nlpar = "thetat")
  expect_equal(result$dpar, "kappa")
  expect_equal(result$nlpar, "thetat")
})

test_that(".bmmfit_resolve_par passes through when both are NULL", {
  mock <- readRDS(test_path("assets/mock_bmmfit_mixture2p.rds"))
  result <- .bmmfit_resolve_par(mock, dpar = NULL, nlpar = NULL)
  expect_null(result$dpar)
  expect_null(result$nlpar)
})


# Tier 2: Fixture-based integration tests (SDM fixture — has posterior samples)

test_that("emmeans() returns valid emmGrid for SDM dpar", {
  skip_if_not_installed("emmeans")
  fit <- readRDS(test_path("assets/bmmfit_example1.rds"))
  em <- emmeans::emmeans(fit, ~ set_size, dpar = "c")
  expect_s4_class(em, "emmGrid")
  expect_true(nrow(as.data.frame(em)) > 0)
})

test_that("pairs() produces valid contrast from emmGrid", {
  skip_if_not_installed("emmeans")
  fit <- readRDS(test_path("assets/bmmfit_example1.rds"))
  em <- emmeans::emmeans(fit, ~ set_size, dpar = "c")
  contrasts <- pairs(em)
  expect_s4_class(contrasts, "emmGrid")
  expect_true(nrow(as.data.frame(contrasts)) > 0)
})
