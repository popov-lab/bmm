# Tests for pp_check.bmmfit() — multinomial model support

load_m3_fit <- function() {
  path <- test_path("assets/bmmfit_m3_ppcheck.rds")
  skip_if_not(file.exists(path), "M3 fixture not available (excluded by .Rbuildignore)")
  readRDS(path)
}

load_sdm_fit <- function() {
  path <- test_path("assets/bmmfit_example1.rds")
  skip_if_not(file.exists(path), "SDM fixture not available (excluded by .Rbuildignore)")
  readRDS(path)
}

test_that("pp_check() returns ggplot for m3 model", {
  fit <- load_m3_fit()
  p <- pp_check(fit, ndraws = 5)
  expect_s3_class(p, "ggplot")
  expect_length(p$layers, 2L)  # geom_col (observed) + geom_pointrange (predicted)
})

test_that("pp_check() aggregates over all predictors by default", {
  fit <- load_m3_fit()
  p <- pp_check(fit, ndraws = 5)
  expect_s3_class(p$facet, "FacetNull")
})

test_that("pp_check() respects ndraws and probs arguments", {
  fit <- load_m3_fit()
  p <- pp_check(fit, ndraws = 8, probs = c(0.10, 0.90))
  expect_match(p$labels$subtitle, "8 draws")
  expect_match(p$labels$subtitle, "10\u201390% CrI")
})

test_that("pp_check() group argument facets by that predictor", {
  fit <- load_m3_fit()
  p <- pp_check(fit, ndraws = 5, group = "cond")
  expect_s3_class(p, "ggplot")
  expect_s3_class(p$facet, "FacetWrap")
})

test_that("pp_check() group with non-predictor column warns", {
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
  fit <- load_m3_fit()
  expect_warning(
    p <- pp_check(fit, ndraws = 5, group = "no_such_col"),
    "not a predictor"
  )
  expect_s3_class(p, "ggplot")
  expect_s3_class(p$facet, "FacetNull")
})

test_that("pp_check() warns when type is specified for multinomial model", {
  fit <- load_m3_fit()
  expect_warning(
    pp_check(fit, type = "hist", ndraws = 5),
    "type.*ignored"
  )
})

test_that("pp_check() uses bayesplot theme_default", {
  fit <- load_m3_fit()
  p <- pp_check(fit, ndraws = 5)
  expect_equal(p$theme$text$family, "serif")
})

test_that("pp_check() uses response category names on x-axis", {
  fit <- load_m3_fit()
  p <- pp_check(fit, ndraws = 5)
  pdata <- p$data
  expect_true("category" %in% names(pdata))
  expect_equal(levels(pdata$category), c("corr", "other", "dist", "npl"))
})

test_that("pp_check() delegates to brms for non-multinomial bmmfit", {
  fit <- load_sdm_fit()
  p <- pp_check(fit, ndraws = 5)
  expect_s3_class(p, "ggplot")
})

test_that("pp_check() accepts re_formula for multinomial model", {
  fit <- load_m3_fit()
  p <- pp_check(fit, ndraws = 5, re_formula = NA)
  expect_s3_class(p, "ggplot")
})

# Auto-grouped type selection for non-multinomial models

test_that("pp_check() auto-selects grouped type when group is provided", {
  fit <- load_sdm_fit()
  # bayesplot warns "group unrecognized" for custom families — upstream issue
  p <- suppressWarnings(pp_check(fit, group = "set_size", ndraws = 5))
  expect_s3_class(p, "ggplot")
})

test_that("pp_check() auto-converts explicit type to grouped variant", {
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


# Multi-observable checks: pp_check(fit, resp_var = ...) (#401)

load_ppcheck_fit <- function(name) {
  path <- test_path("assets", name)
  skip_if_not(file.exists(path), "fixture not available (excluded by .Rbuildignore)")
  readRDS(path)
}

test_that("pp_check(resp_var = 'rt') checks the observed response times", {
  fit <- load_ppcheck_fit("bmmfit_ddm_ppcheck.rds")
  p <- pp_check(fit, resp_var = "rt", ndraws = 5)
  expect_s3_class(p, "ggplot")
  expect_setequal(p$data$value[p$data$is_y_label == "italic(y)"], fit$data$rt)
})

test_that("pp_check(resp_var = 'signed_rt') signs RTs by the response", {
  fit <- load_ppcheck_fit("bmmfit_ddm_ppcheck.rds")
  p <- pp_check(fit, resp_var = "signed_rt", ndraws = 5)
  expect_s3_class(p, "ggplot")
  y <- p$data$value[p$data$is_y_label == "italic(y)"]
  expect_setequal(y, fit$data$rt * (2 * fit$data$response - 1))
  expect_true(any(y < 0))
})

test_that("pp_check(negative_rt = TRUE) redirects to the signed_rt check", {
  fit <- load_ppcheck_fit("bmmfit_ddm_ppcheck.rds")
  set.seed(99)
  expect_message(p1 <- pp_check(fit, negative_rt = TRUE, ndraws = 5),
                 "signed_rt")
  set.seed(99)
  p2 <- pp_check(fit, resp_var = "signed_rt", ndraws = 5)
  expect_equal(p1$data, p2$data)
})

test_that("pp_check(resp_var = 'all') panels all checks", {
  fit <- load_ppcheck_fit("bmmfit_ddm_ppcheck.rds")
  p <- pp_check(fit, resp_var = "all", ndraws = 5)
  expect_s3_class(p, "bayesplot_grid")
  expect_length(p$bayesplots, 3L)
  expect_warning(pp_check(fit, resp_var = "all", type = "hist", ndraws = 5),
                 "ignored")
})

test_that("pp_check(resp_var) supports grouping and type overrides", {
  fit <- load_ppcheck_fit("bmmfit_ddm_ppcheck.rds")
  p <- pp_check(fit, resp_var = "rt", group = "cond", ndraws = 5)
  expect_s3_class(p, "ggplot")
  expect_s3_class(p$facet, "FacetWrap")
  expect_s3_class(pp_check(fit, resp_var = "rt", type = "hist", ndraws = 5),
                  "ggplot")
  expect_error(pp_check(fit, resp_var = "rt", type = "no_such_type"),
               "not a supported")
  expect_error(pp_check(fit, resp_var = "rt", type = "data"),
               "not a supported")
  expect_error(pp_check(fit, resp_var = "rt", group = "no_such_col"),
               "must name a column")
})

test_that("pp_check(resp_var) accepts draw_ids and rejects newdata", {
  fit <- load_ppcheck_fit("bmmfit_ddm_ppcheck.rds")
  expect_s3_class(pp_check(fit, resp_var = "rt", draw_ids = 1:3), "ggplot")
  expect_error(pp_check(fit, resp_var = "rt", newdata = fit$data),
               "not supported")
})

test_that("pp_check(resp_var) works for the 3par ezdm model", {
  fit <- load_ppcheck_fit("bmmfit_ezdm3_ppcheck.rds")
  p <- pp_check(fit, resp_var = "mean_pc", ndraws = 5)
  expect_s3_class(p, "ggplot")
  y <- p$data$value[p$data$is_y_label == "italic(y)"]
  expect_setequal(y, fit$data$n_upper / fit$data$n_trials)
  expect_s3_class(pp_check(fit, resp_var = "var_rt", ndraws = 5), "ggplot")
})

test_that("pp_check(resp_var) drops undefined 4par ezdm cells with a warning", {
  fit <- load_ppcheck_fit("bmmfit_ezdm4_ppcheck.rds")
  expect_warning(p <- pp_check(fit, resp_var = "mean_rt_lower", ndraws = 20),
                 "Dropped")
  expect_s3_class(p, "ggplot")
  expect_s3_class(suppressWarnings(pp_check(fit, resp_var = "all", ndraws = 5)),
                  "bayesplot_grid")
})
