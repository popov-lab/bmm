simulate_sdm_smoke_data <- function() {
  set.seed(123)
  condition <- factor(rep(paste0("condition", 1:2), each = 50))
  data.frame(
    y = rsdm(100),
    condition = condition
  )
}

test_that("SDM benchmark simulator output works with current SDM workflow", {
  dat <- simulate_sdm_smoke_data()
  formula <- bmf(c ~ 0 + condition, kappa ~ 0 + condition)

  sdata <- standata(formula, data = dat, model = sdm(resp_error = "y"))

  expect_equal(nrow(dat), 100L)
  expect_equal(sdata$N, nrow(dat))
  expect_equal(sdata$G_sdm_runs, 2L)
  expect_equal(as.integer(sdata$sdm_run_count), c(50L, 50L))
  expect_true(all(dat$y >= -pi & dat$y <= pi))
})

test_that("SDM generated Stan code includes custom SDM chunks", {
  dat <- simulate_sdm_smoke_data()
  formula <- bmf(c ~ 0 + condition, kappa ~ 0 + condition)

  code <- stancode(formula, data = dat, model = sdm(resp_error = "y"))

  expect_match(code, "sdm_simple_ldenom_chquad_adaptive", fixed = TRUE)
  expect_match(code, "sdm_simple_run_ldenom", fixed = TRUE)
  expect_match(code, "target += sdm_simple_run_ldenom", fixed = TRUE)
  expect_false(grepl("c[n] != c[n-1]", code, fixed = TRUE))
  expect_match(code, "COSN", fixed = TRUE)
})

test_that("SDM threaded Stan code slices denominator runs correctly", {
  dat <- simulate_sdm_smoke_data()
  formula <- bmf(c ~ 0 + condition, kappa ~ 0 + condition)

  code <- stancode(
    formula,
    data = dat,
    model = sdm(resp_error = "y"),
    threads = brms::threading(2, grainsize = 1)
  )

  expect_match(code, "target += reduce_sum", fixed = TRUE)
  expect_match(code, "sdm_simple_run_ldenom_slice", fixed = TRUE)
  expect_match(code, "run_start[g] >= start", fixed = TRUE)
  expect_match(code, "run_start[g] <= end", fixed = TRUE)
  expect_match(code, "run_start[g] - start + 1", fixed = TRUE)
})
