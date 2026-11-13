fit <- readRDS(here::here("tests/testthat/assets/bmmfit_mixture3p.rds"))

test_that("summary returns brmssummary with backend brms", {
  out <- summary(fit, backend = "brms")
  expect_true(inherits(out, "brmssummary"))
})

test_that("summary returns bmmsummary with backend bmm", {
  out <- summary(fit, backend = "bmm")
  expect_true(inherits(out, "bmmsummary"))
})
