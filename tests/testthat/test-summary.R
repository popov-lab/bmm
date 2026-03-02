test_that("summary has reasonable outputs", {
  skip_if_not(interactive())
  fit <- readRDS(test_path("assets/bmmfit_example1.rds"))
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
