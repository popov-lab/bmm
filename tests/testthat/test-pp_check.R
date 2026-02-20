# Tests for pp_check.bmmfit() — multinomial model support

load_m3_fit <- function() {
  readRDS(system.file("extdata", "bmmfit_m3_ppcheck.rds", package = "bmm"))
}

load_sdm_fit <- function() {
  readRDS(system.file("extdata", "bmmfit_example1.rds", package = "bmm"))
}

test_that("pp_check() returns ggplot for m3 model", {
  skip_if_not_installed("ggplot2")
  fit <- load_m3_fit()
  p <- pp_check(fit, ndraws = 5)
  expect_s3_class(p, "ggplot")
  expect_length(p$layers, 2L)  # geom_col (observed) + geom_pointrange (predicted)
})

test_that("pp_check() aggregates over all predictors by default", {
  skip_if_not_installed("ggplot2")
  fit <- load_m3_fit()
  p <- pp_check(fit, ndraws = 5)
  expect_s3_class(p$facet, "FacetNull")
})

test_that("pp_check() respects ndraws and probs arguments", {
  skip_if_not_installed("ggplot2")
  fit <- load_m3_fit()
  p <- pp_check(fit, ndraws = 8, probs = c(0.10, 0.90))
  expect_match(p$labels$subtitle, "8 draws")
  expect_match(p$labels$subtitle, "10\u201390% CrI")
})

test_that("pp_check() group argument facets by that predictor", {
  skip_if_not_installed("ggplot2")
  fit <- load_m3_fit()
  p <- pp_check(fit, ndraws = 5, group = "cond")
  expect_s3_class(p, "ggplot")
  expect_s3_class(p$facet, "FacetWrap")
})

test_that("pp_check() group with non-predictor column warns", {
  skip_if_not_installed("ggplot2")
  fit <- load_m3_fit()
  # ID is a random-effect grouping var, not a predictor
  expect_warning(
    p <- pp_check(fit, ndraws = 5, group = "ID"),
    "not a predictor"
  )
  expect_s3_class(p, "ggplot")
  expect_s3_class(p$facet, "FacetNull")
})

test_that("pp_check() group with unknown column warns", {
  skip_if_not_installed("ggplot2")
  fit <- load_m3_fit()
  expect_warning(
    p <- pp_check(fit, ndraws = 5, group = "no_such_col"),
    "not a predictor"
  )
  expect_s3_class(p, "ggplot")
  expect_s3_class(p$facet, "FacetNull")
})

test_that("pp_check() warns when type is specified for multinomial model", {
  skip_if_not_installed("ggplot2")
  fit <- load_m3_fit()
  expect_warning(
    pp_check(fit, type = "hist", ndraws = 5),
    "type.*ignored"
  )
})

test_that("pp_check() uses bayesplot theme_default", {
  skip_if_not_installed("ggplot2")
  fit <- load_m3_fit()
  p <- pp_check(fit, ndraws = 5)
  expect_equal(p$theme$text$family, "serif")
})

test_that("pp_check() uses response category names on x-axis", {
  skip_if_not_installed("ggplot2")
  fit <- load_m3_fit()
  p <- pp_check(fit, ndraws = 5)
  pdata <- p$data
  expect_true("category" %in% names(pdata))
  expect_equal(levels(pdata$category), c("corr", "other", "dist", "npl"))
})

test_that("pp_check() delegates to brms for non-multinomial bmmfit", {
  skip_if_not_installed("ggplot2")
  fit <- load_sdm_fit()
  p <- pp_check(fit, ndraws = 5)
  expect_s3_class(p, "ggplot")
})

test_that("pp_check() accepts re_formula for multinomial model", {
  skip_if_not_installed("ggplot2")
  fit <- load_m3_fit()
  p <- pp_check(fit, ndraws = 5, re_formula = NA)
  expect_s3_class(p, "ggplot")
})

# Auto-grouped type selection for non-multinomial models

test_that("pp_check() auto-selects grouped type when group is provided", {
  skip_if_not_installed("ggplot2")
  fit <- load_sdm_fit()
  # bayesplot warns "group unrecognized" for custom families — upstream issue
  p <- suppressWarnings(pp_check(fit, group = "set_size", ndraws = 5))
  expect_s3_class(p, "ggplot")
})

test_that("pp_check() auto-converts explicit type to grouped variant", {
  skip_if_not_installed("ggplot2")
  fit <- load_sdm_fit()
  p <- pp_check(fit, type = "stat", group = "set_size", ndraws = 5)
  expect_s3_class(p, "ggplot")
})

test_that(".auto_grouped_type() appends _grouped when variant exists", {
  expect_equal(.auto_grouped_type("dens_overlay"), "dens_overlay_grouped")
  expect_equal(.auto_grouped_type("stat"), "stat_grouped")
  expect_equal(.auto_grouped_type("violin"), "violin_grouped")
})

test_that(".auto_grouped_type() does not double-append _grouped", {
  expect_equal(.auto_grouped_type("dens_overlay_grouped"), "dens_overlay_grouped")
})

test_that(".auto_grouped_type() returns type unchanged when no grouped variant", {
  expect_equal(.auto_grouped_type("hist"), "hist")
  expect_equal(.auto_grouped_type("loo_pit"), "loo_pit")
})

test_that(".resolve_pp_conditions() excludes infrastructure columns", {
  fit <- load_m3_fit()
  conds <- .resolve_pp_conditions(fit)
  expect_true("cond" %in% conds)
  expect_false("nTrials" %in% conds)
  expect_false("ID" %in% conds)
  expect_false(any(grepl("^Idx_", conds)))
  expect_false(any(grepl("^n_", conds)))
})
