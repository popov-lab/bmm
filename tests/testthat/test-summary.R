test_that("summary has reasonable outputs", {
  skip_on_cran()
  path <- test_path("assets/bmmfit_example1.rds")
  skip_if_not(file.exists(path), "SDM fixture not available (excluded by .Rbuildignore)")
  fit <- readRDS(path)
  summary1 <- suppressWarnings(summary(fit))
  expect_true(is.data.frame(summary1$fixed))
  expect_equal(
    rownames(summary1$fixed),
    c("mu_Intercept", "kappa_Intercept", "c_set_size1", "c_set_size2", "c_set_size3", "c_set_size4")
  )
  expect_equal(
    colnames(summary1$fixed),
    c("Estimate", "Est.Error", "l-95% CI", "u-95% CI", "Rhat", "Bulk_ESS", "Tail_ESS")
  )
  expect_output(print(summary1), "Constant Parameters:")
  expect_output(print(summary1), "Model: sdm")
  expect_output(print(summary1), "Links: mu = tan_half; c = log; kappa = log")
  expect_output(print(summary1), "Formula: mu = 0")
})

test_that(".summary_fixed_rows keeps a single fixed-effect row", {
  # gumbel-min sdt_ranking has one population coefficient (dprime_Intercept) but
  # two printed parameters (dprime, sdratio); the old sapply+apply errored here.
  one_row <- data.frame(Estimate = 0.6, Rhat = 1, row.names = "dprime_Intercept")
  out <- .summary_fixed_rows(one_row, c("dprime", "sdratio"))
  expect_s3_class(out, "data.frame")
  expect_identical(rownames(out), "dprime_Intercept")
})

test_that(".summary_fixed_rows selects all rows matching the printed parameters", {
  fixed <- data.frame(
    Estimate = 1:3, Rhat = c(1, NA, 1),
    row.names = c("dprime_Intercept", "sdratio_Intercept", "nuisance_Intercept")
  )
  out <- .summary_fixed_rows(fixed, c("dprime", "sdratio"))
  expect_identical(sort(rownames(out)), c("dprime_Intercept", "sdratio_Intercept"))
})
