# Needs sampling: run by hand with
#   NOT_CRAN=true testthat::test_file("tests/internal/test-kfold-update.R")
# kfold() and reloo() refit with update(newdata = <slice of fit$data>)

skip_if_not_installed("cmdstanr")
skip_if(is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)), "CmdStan is not installed")

fit_counts <- function(...) {
  suppressWarnings(bmm(
    ..., backend = "cmdstanr", chains = 1, iter = 300, warmup = 150,
    refresh = 0, silent = 2
  ))
}

rating_data <- data.frame(
  stimulus = rep(c(0L, 1L), 6), id = factor(rep(1:6, each = 2)),
  r1 = c(30, 8, 26, 10, 33, 6, 28, 9, 31, 7, 27, 11),
  r2 = c(25, 12, 27, 14, 22, 11, 24, 13, 23, 12, 26, 15),
  r3 = c(20, 15, 21, 16, 19, 17, 22, 15, 20, 16, 21, 14),
  r4 = c(15, 25, 16, 22, 17, 28, 14, 26, 16, 24, 15, 23),
  r5 = c(10, 40, 10, 38, 9, 38, 12, 37, 10, 41, 11, 37)
)
rating_fit <- fit_counts(
  bmf(d ~ 1 + (1 | id), criterion ~ 1, spacing ~ 1, sdratio ~ 1),
  rating_data, sdt_rating(paste0("r", 1:5), "stimulus")
)

m3_cats <- c("corr", "other", "dist", "npl")
m3_ids <- unique(oberauer_lewandowsky_2019_e1$ID)[1:6]
m3_fit <- fit_counts(
  bmf(corr ~ b + a + c, other ~ b + a, dist ~ b + d, npl ~ b,
      c ~ 1 + (1 | ID), a ~ 1, d ~ 1),
  oberauer_lewandowsky_2019_e1[oberauer_lewandowsky_2019_e1$ID %in% m3_ids, ],
  m3(resp_cats = m3_cats, num_options = paste0("n_", m3_cats), choice_rule = "simple",
     links = list(c = "log", a = "log", d = "log"))
)

fits <- list(sdt_rating = list(fit = rating_fit, group = "id"), m3 = list(fit = m3_fit, group = "ID"))

for (name in names(fits)) {
  fit <- fits[[name]]$fit

  test_that(paste("update(newdata = fit$data[-1, ]) refits a", name, "fit"), {
    up <- suppressWarnings(update(fit, newdata = fit$data[-1, ], refresh = 0, silent = 2))
    expect_s3_class(up, "bmmfit")
    expect_equal(nrow(up$data), nrow(fit$data) - 1)
  })

  test_that(paste("kfold() with grouped folds runs on a", name, "fit"), {
    kf <- suppressWarnings(brms::kfold(
      fit, folds = "grouped", group = fits[[name]]$group, K = 2, cores = 1,
      refresh = 0, silent = 2
    ))
    expect_s3_class(kf, "kfold")
    expect_true(is.finite(kf$estimates["elpd_kfold", "Estimate"]))
    expect_equal(length(kf$pointwise[, "elpd_kfold"]), nrow(fit$data))
  })
}
