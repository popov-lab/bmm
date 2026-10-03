# Tier 1: Unit tests for .bmmfit_resolve_par() using mock fixtures
# These test the parameter routing logic without requiring posterior samples

load_mock_mixture2p <- function() {
  path <- test_path("assets/mock_bmmfit_mixture2p.rds")
  skip_if_not(file.exists(path), "fixture not available (excluded by .Rbuildignore)")
  readRDS(path)
}

load_mock_m3 <- function() {
  path <- test_path("assets/mock_bmmfit_m3.rds")
  skip_if_not(file.exists(path), "fixture not available (excluded by .Rbuildignore)")
  readRDS(path)
}

load_sdm_fit <- function() {
  path <- test_path("assets/bmmfit_example1.rds")
  skip_if_not(file.exists(path), "SDM fixture not available (excluded by .Rbuildignore)")
  readRDS(path)
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

# posterior_epred (#475) ------------------------------------------------------

test_that("emmeans(epred = TRUE) works on a fit saved without posterior_epred", {
  skip_if_not_installed("emmeans")
  skip_on_cran()
  fit <- load_fixture_fit("bmmfit_ddm_ppcheck.rds")
  em <- as.data.frame(emmeans::emmeans(fit, ~1, epred = TRUE))
  expect_true(is.finite(em$emmean))
  expect_equal(em$emmean, mean(brms::posterior_epred(fit)), tolerance = 0.05)
})

test_that("emmeans(epred = TRUE) refuses a circular model", {
  skip_if_not_installed("emmeans")
  fit <- load_sdm_fit()
  class(fit$bmm$model) <- c("bmmodel", "circular", "sdm", "sdm_simple")
  # a stand-in for the number brms would compute for a native family, so that
  # only the refusal before brms can produce the expected message
  fit$formula$family$posterior_epred <- function(prep) stop("reached brms")
  expect_error(emmeans::emmeans(fit, ~set_size, epred = TRUE),
               "not defined for the sdm model")
})
