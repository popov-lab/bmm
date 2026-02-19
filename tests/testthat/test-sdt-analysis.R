# Tests for roc_sdt() and auc_sdt() in helpers-analysis.R

load_sdt_fit <- function(name) {
  # devtools::test() sets wd to tests/testthat/, so navigate up to package root
  path <- test_path("..", "..", "vignettes", "articles", "assets",
                    paste0(name, ".rds"))
  if (!file.exists(path)) skip("SDT fit file not found")
  readRDS(path)
}


############################################################################# !
# INPUT VALIDATION                                                       ####
############################################################################# !

test_that("roc_sdt() errors on non-bmmfit input", {
  expect_error(roc_sdt(list()), "bmmfit object")
})

test_that("roc_sdt() errors on non-SDT bmmfit", {
  fit <- load_sdt_fit("bmmfit_sdt_binary")
  # Temporarily patch model class to simulate a non-SDT model
  fake_fit <- fit
  class(fake_fit$bmm$model) <- c("bmmodel", "sdm")
  expect_error(roc_sdt(fake_fit), "only available for SDT")
})

test_that("roc_sdt() errors for m-AFC and ranking models", {
  fit <- load_sdt_fit("bmmfit_sdt_binary")

  fake_mafc    <- fit
  fake_mafc$bmm$model$version <- "mafc"
  expect_error(roc_sdt(fake_mafc), "not defined")

  fake_ranking <- fit
  fake_ranking$bmm$model$version <- "ranking"
  expect_error(roc_sdt(fake_ranking), "not defined")
})

test_that("auc_sdt() errors for m-AFC and ranking models", {
  fit <- load_sdt_fit("bmmfit_sdt_binary")

  fake_mafc <- fit
  fake_mafc$bmm$model$version <- "mafc"
  expect_error(auc_sdt(fake_mafc), "not defined")
})


############################################################################# !
# ROC — BINARY MODEL                                                     ####
############################################################################# !

test_that("roc_sdt() returns bmm_sdt_roc for binary fit", {
  fit <- load_sdt_fit("bmmfit_sdt_binary")
  roc <- roc_sdt(fit, n_points = 20)
  expect_s3_class(roc, "bmm_sdt_roc")
  expect_s3_class(roc, "data.frame")
})

test_that("roc_sdt() binary output has correct columns", {
  fit <- load_sdt_fit("bmmfit_sdt_binary")
  roc <- roc_sdt(fit, n_points = 20)
  expect_true(all(c("FA", "Hit", ".draw") %in% names(roc)))
})

test_that("roc_sdt() binary FA and Hit are in [0, 1]", {
  fit <- load_sdt_fit("bmmfit_sdt_binary")
  roc <- roc_sdt(fit, n_points = 20)
  expect_true(all(roc$FA  >= 0 & roc$FA  <= 1))
  expect_true(all(roc$Hit >= 0 & roc$Hit <= 1))
})

test_that("roc_sdt() binary includes trivial endpoints (0,0) and (1,1)", {
  fit <- load_sdt_fit("bmmfit_sdt_binary")
  roc <- roc_sdt(fit, n_points = 20)
  # Every draw should start at (0,0) and end at (1,1)
  draw1 <- roc[roc$.draw == 1L, ]
  expect_true(any(draw1$FA == 0 & draw1$Hit == 0))
  expect_true(any(draw1$FA == 1 & draw1$Hit == 1))
})

test_that("roc_sdt() binary has n_points + 2 rows per draw × condition", {
  fit    <- load_sdt_fit("bmmfit_sdt_binary")
  n_pts  <- 20L
  roc    <- roc_sdt(fit, n_points = n_pts)
  n_draws <- length(unique(roc$.draw))
  n_cond  <- if ("condition" %in% names(roc)) length(unique(roc$condition)) else 1L
  expect_equal(nrow(roc), n_draws * n_cond * (n_pts + 2L))
})

test_that("roc_sdt() binary summary attribute exists and is valid", {
  fit  <- load_sdt_fit("bmmfit_sdt_binary")
  roc  <- roc_sdt(fit, n_points = 20)
  summ <- attr(roc, "summary")
  expect_true(is.data.frame(summ))
  expect_true(all(c("FA", "Hit_mean", "Hit_lower", "Hit_upper") %in% names(summ)))
  expect_true(all(summ$Hit_lower <= summ$Hit_mean))
  expect_true(all(summ$Hit_mean  <= summ$Hit_upper))
})

test_that("roc_sdt() binary respects conditions argument", {
  fit  <- load_sdt_fit("bmmfit_sdt_binary")
  cond <- data.frame(condition = "easy")
  roc  <- roc_sdt(fit, conditions = cond, n_points = 10)
  expect_true("condition" %in% names(roc))
  expect_equal(unique(roc$condition), "easy")
})

test_that("roc_sdt() binary: Hit increases monotonically with FA (median)", {
  fit  <- load_sdt_fit("bmmfit_sdt_binary")
  roc  <- roc_sdt(fit, n_points = 30)
  summ <- attr(roc, "summary")
  # For each condition group, median Hit should be non-decreasing in FA
  cond_groups <- if ("condition" %in% names(summ)) unique(summ$condition) else "all"
  for (cg in cond_groups) {
    if ("condition" %in% names(summ)) {
      sub <- summ[summ$condition == cg, ]
    } else {
      sub <- summ
    }
    sub <- sub[order(sub$FA), ]
    expect_true(all(diff(sub$Hit_mean) >= -1e-6),
                info = paste("Monotonicity violated for condition:", cg))
  }
})


############################################################################# !
# ROC — RATING MODEL                                                     ####
############################################################################# !

test_that("roc_sdt() returns bmm_sdt_roc for rating fit", {
  fit <- load_sdt_fit("bmmfit_sdt_rating")
  roc <- roc_sdt(fit)
  expect_s3_class(roc, "bmm_sdt_roc")
})

test_that("roc_sdt() rating has K+1 points per draw × condition", {
  fit       <- load_sdt_fit("bmmfit_sdt_rating")
  n_ratings <- fit$bmm$model$other_vars$n_ratings  # 6
  roc       <- roc_sdt(fit)
  n_draws   <- length(unique(roc$.draw))
  n_cond    <- max(1L, nrow(.sdt_resolve_conditions(fit, NULL)))
  expected_rows_per_draw <- (n_ratings + 1L) * n_cond
  expect_equal(nrow(roc), n_draws * expected_rows_per_draw)
})

test_that("roc_sdt() rating FA and Hit are in [0, 1]", {
  fit <- load_sdt_fit("bmmfit_sdt_rating")
  roc <- roc_sdt(fit)
  expect_true(all(roc$FA  >= 0 & roc$FA  <= 1))
  expect_true(all(roc$Hit >= 0 & roc$Hit <= 1))
})

test_that("roc_sdt() rating includes endpoints (0,0) and (1,1)", {
  fit <- load_sdt_fit("bmmfit_sdt_rating")
  roc <- roc_sdt(fit)
  draw1 <- roc[roc$.draw == 1L, ]
  expect_true(any(abs(draw1$FA) < 1e-8 & abs(draw1$Hit) < 1e-8))
  expect_true(any(abs(draw1$FA - 1) < 1e-8 & abs(draw1$Hit - 1) < 1e-8))
})

test_that("roc_sdt() gumbel rating works without error", {
  fit <- load_sdt_fit("bmmfit_sdt_gumbel")
  expect_s3_class(roc_sdt(fit), "bmm_sdt_roc")
})


############################################################################# !
# AUC — BINARY MODEL (ANALYTICAL)                                        ####
############################################################################# !

test_that("auc_sdt() returns bmm_sdt_auc for binary fit", {
  fit <- load_sdt_fit("bmmfit_sdt_binary")
  auc <- auc_sdt(fit)
  expect_s3_class(auc, "bmm_sdt_auc")
  expect_s3_class(auc, "data.frame")
})

test_that("auc_sdt() binary output has correct columns", {
  fit <- load_sdt_fit("bmmfit_sdt_binary")
  auc <- auc_sdt(fit)
  expect_true(all(c("AUC", ".draw") %in% names(auc)))
})

test_that("auc_sdt() binary AUC values are in (0.5, 1)", {
  fit <- load_sdt_fit("bmmfit_sdt_binary")
  auc <- auc_sdt(fit)
  # dprime is positive (easy + hard conditions both > 0) so AUC > 0.5
  expect_true(all(auc$AUC > 0.5))
  expect_true(all(auc$AUC < 1))
})

test_that("auc_sdt() normal binary uses analytical formula Phi(d/sqrt(2))", {
  fit <- load_sdt_fit("bmmfit_sdt_binary")
  # Analytical AUC for normal EV-SDT = Phi(dprime / sqrt(2))
  auc  <- auc_sdt(fit, conditions = data.frame(condition = "easy"))
  mean_auc <- mean(auc$AUC)
  # Get corresponding dprime draws using a complete newdata row
  easy_row <- fit$data[fit$data$condition == "easy", ][1L, ]
  dp_mat <- brms::posterior_linpred(fit, dpar = "dprime",
                                    newdata = easy_row,
                                    re_formula = NA)
  expected_auc <- mean(pnorm(dp_mat[, 1L] / sqrt(2)))
  expect_equal(mean_auc, expected_auc, tolerance = 1e-10)
})

test_that("auc_sdt() binary summary attribute is valid", {
  fit  <- load_sdt_fit("bmmfit_sdt_binary")
  auc  <- auc_sdt(fit)
  summ <- attr(auc, "summary")
  expect_true(is.data.frame(summ))
  expect_true(all(c("AUC_mean", "AUC_lower", "AUC_upper") %in% names(summ)))
  expect_true(all(summ$AUC_lower <= summ$AUC_mean))
  expect_true(all(summ$AUC_mean  <= summ$AUC_upper))
})


############################################################################# !
# AUC — RATING MODEL (NUMERICAL)                                         ####
############################################################################# !

test_that("auc_sdt() works for rating fit (numerical path)", {
  fit <- load_sdt_fit("bmmfit_sdt_rating")
  auc <- auc_sdt(fit)
  expect_s3_class(auc, "bmm_sdt_auc")
  expect_true(all(auc$AUC > 0.5 & auc$AUC < 1))
})

test_that("auc_sdt() gumbel rating: AUC > 0.5 (positive dprime)", {
  fit <- load_sdt_fit("bmmfit_sdt_gumbel")
  auc <- auc_sdt(fit)
  expect_true(mean(auc$AUC) > 0.5)
})


############################################################################# !
# PRINT AND PLOT METHODS                                                 ####
############################################################################# !

test_that("print.bmm_sdt_roc() works without error", {
  fit <- load_sdt_fit("bmmfit_sdt_binary")
  roc <- roc_sdt(fit, n_points = 10)
  expect_output(print(roc), "SDT ROC curve")
})

test_that("print.bmm_sdt_auc() works without error", {
  fit <- load_sdt_fit("bmmfit_sdt_binary")
  auc <- auc_sdt(fit)
  expect_output(print(auc), "SDT AUC")
})

test_that("plot.bmm_sdt_roc() returns a ggplot object", {
  skip_if_not_installed("ggplot2")
  fit <- load_sdt_fit("bmmfit_sdt_binary")
  roc <- roc_sdt(fit, n_points = 10)
  p <- plot(roc)
  expect_s3_class(p, "ggplot")
})

test_that("plot.bmm_sdt_auc() returns a ggplot object", {
  skip_if_not_installed("ggplot2")
  fit <- load_sdt_fit("bmmfit_sdt_binary")
  auc <- auc_sdt(fit)
  p <- plot(auc)
  expect_s3_class(p, "ggplot")
})

test_that("plot.bmm_sdt_roc() works for rating fit", {
  skip_if_not_installed("ggplot2")
  fit <- load_sdt_fit("bmmfit_sdt_rating")
  roc <- roc_sdt(fit)
  p   <- plot(roc)
  expect_s3_class(p, "ggplot")
})


# ---------------------------------------------------------------------------
# pp_check
# ---------------------------------------------------------------------------

test_that("pp_check.bmmfit() returns ggplot for rating SDT model", {
  skip_if_not_installed("ggplot2")
  fit <- load_sdt_fit("bmmfit_sdt_rating")
  p <- pp_check(fit, ndraws = 5)
  expect_s3_class(p, "ggplot")
  # 4 layers: ribbon, line (predicted), point + line (observed)
  expect_length(p$layers, 4L)
})

test_that("pp_check.bmmfit() facets by stimulus for intercept-only rating fit", {
  skip_if_not_installed("ggplot2")
  fit <- load_sdt_fit("bmmfit_sdt_rating")
  p <- pp_check(fit, ndraws = 5)
  expect_s3_class(p$facet, "FacetWrap")
})

test_that("pp_check.bmmfit() respects ndraws and probs arguments", {
  skip_if_not_installed("ggplot2")
  fit <- load_sdt_fit("bmmfit_sdt_rating")
  p <- pp_check(fit, ndraws = 8, probs = c(0.10, 0.90))
  expect_match(p$labels$subtitle, "8 draws")
  expect_match(p$labels$subtitle, "10\u201390%")
})

test_that("pp_check.bmmfit() group adds facet dimension for rating model", {
  skip_if_not_installed("ggplot2")
  fit <- load_sdt_fit("bmmfit_sdt_rating")
  # id is a random-effect grouping var excluded by default; group = "id" re-adds it
  p <- pp_check(fit, ndraws = 5, group = "id")
  expect_s3_class(p, "ggplot")
  # With group, faceting must span stimulus × id → FacetGrid
  expect_s3_class(p$facet, "FacetGrid")
})

test_that("pp_check.bmmfit() group = column already in default grouping is a no-op", {
  skip_if_not_installed("ggplot2")
  fit <- load_sdt_fit("bmmfit_sdt_rating")
  p_default <- pp_check(fit, ndraws = 5)
  p_group   <- pp_check(fit, ndraws = 5, group = "stimulus")
  # stimulus is already the row dimension; passing it again should not change structure
  expect_s3_class(p_group$facet, class(p_default$facet)[1L])
})

test_that("pp_check.bmmfit() group with unknown column warns and falls back", {
  skip_if_not_installed("ggplot2")
  fit <- load_sdt_fit("bmmfit_sdt_rating")
  expect_warning(
    p <- pp_check(fit, ndraws = 5, group = "no_such_col"),
    "not found in model data"
  )
  expect_s3_class(p, "ggplot")
})

test_that("pp_check.bmmfit() warns when type is specified for rating model", {
  skip_if_not_installed("ggplot2")
  fit <- load_sdt_fit("bmmfit_sdt_rating")
  expect_warning(
    pp_check(fit, type = "dens_overlay", ndraws = 5),
    "type.*ignored"
  )
})

test_that("pp_check.bmmfit() delegates to brms for binary SDT model", {
  skip_if_not_installed("ggplot2")
  fit <- load_sdt_fit("bmmfit_sdt_binary")
  p <- pp_check(fit, ndraws = 5)
  expect_s3_class(p, "ggplot")
})
