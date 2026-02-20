# Tests for pp_check.bmmfit() — multinomial model support

load_m3_fit <- function() {
  readRDS(test_path("../../vignettes/articles/assets/bmmfit_m3_vignette.rds"))
}

test_that("pp_check() returns ggplot for m3 model", {
  skip_if_not_installed("ggplot2")
  fit <- load_m3_fit()
  p <- pp_check(fit, ndraws = 5)
  expect_s3_class(p, "ggplot")
  # 4 layers: ribbon, predicted line, observed points, observed dashed line
  expect_length(p$layers, 4L)
})

test_that("pp_check() facets by condition for m3 model", {
  skip_if_not_installed("ggplot2")
  fit <- load_m3_fit()
  p <- pp_check(fit, ndraws = 5)
  expect_s3_class(p$facet, "FacetWrap")
})

test_that("pp_check() respects ndraws and probs arguments", {
  skip_if_not_installed("ggplot2")
  fit <- load_m3_fit()
  p <- pp_check(fit, ndraws = 8, probs = c(0.10, 0.90))
  expect_match(p$labels$subtitle, "8 draws")
  expect_match(p$labels$subtitle, "10\u201390%")
})

test_that("pp_check() group argument adds facet dimension", {
  skip_if_not_installed("ggplot2")
  fit <- load_m3_fit()
  p <- pp_check(fit, ndraws = 5, group = "ID")
  expect_s3_class(p, "ggplot")
  expect_s3_class(p$facet, "FacetWrap")
})

test_that("pp_check() group with unknown column warns and falls back", {
  skip_if_not_installed("ggplot2")
  fit <- load_m3_fit()
  expect_warning(
    p <- pp_check(fit, ndraws = 5, group = "no_such_col"),
    "not found in model data"
  )
  expect_s3_class(p, "ggplot")
})

test_that("pp_check() warns when type is specified for multinomial model", {
  skip_if_not_installed("ggplot2")
  fit <- load_m3_fit()
  expect_warning(
    pp_check(fit, type = "dens_overlay", ndraws = 5),
    "type.*ignored"
  )
})

test_that("pp_check() uses response category names on x-axis", {
  skip_if_not_installed("ggplot2")
  fit <- load_m3_fit()
  p <- pp_check(fit, ndraws = 5)
  plot_data <- ggplot2::ggplot_build(p)$layout$layout
  # Check that category is a factor with the right levels in the plot data
  pdata <- p$data
  expect_true("category" %in% names(pdata))
  expect_equal(levels(pdata$category), c("corr", "other", "dist", "npl"))
})

test_that("pp_check() delegates to brms for non-multinomial bmmfit", {
  skip_if_not_installed("ggplot2")
  fit <- readRDS(test_path("assets/bmmfit_example1.rds"))
  p <- pp_check(fit, ndraws = 5)
  expect_s3_class(p, "ggplot")
})

test_that(".resolve_pp_conditions() excludes infrastructure columns", {
  fit <- load_m3_fit()
  conds <- .resolve_pp_conditions(fit)
  # Should include cond but not Y, nTrials, Idx_*, n_*, or ID (random effect)
  expect_true("cond" %in% conds)
  expect_false("nTrials" %in% conds)
  expect_false("ID" %in% conds)
  expect_false(any(grepl("^Idx_", conds)))
  expect_false(any(grepl("^n_", conds)))
})
