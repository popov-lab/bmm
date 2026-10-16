fits <- readRDS(here::here('tests/internal/ref_fits0301/fit_objects_0.4.1.9000_seed-365_20240215.rds'))

test_that("summary returns brmssummary with backend brms", {
  out <- summary(fits[[1]], backend = "brms")
  expect_true(inherits(out, "brmssummary"))
})


test_that("summary returns bmmsummary with backend bmm", {
  out <- summary(fits[[1]], backend = "bmm")
  expect_true(inherits(out, "bmmsummary"))
})




