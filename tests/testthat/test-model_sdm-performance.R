source(test_path("..", "..", "inst", "benchmarks", "sdm", "simulate_sdm_data.R"))

test_that("SDM benchmark simulator output works with current SDM workflow", {
  sim <- simulate_sdm_benchmark_data(config = sdm_benchmark_configs()$smoke)
  formula <- bmf(c ~ 0 + condition, kappa ~ 0 + condition)

  sdata <- standata(formula, data = sim$data, model = sdm(resp_error = "y"))

  expect_equal(nrow(sim$data), 100L)
  expect_equal(nrow(sim$true_parameters), 10L)
  expect_equal(sdata$N, nrow(sim$data))
  expect_equal(sdata$G_sdm_runs, 2L)
  expect_equal(as.integer(sdata$sdm_run_count), c(50L, 50L))
  expect_true(all(sim$data$y >= -pi & sim$data$y <= pi))
})

test_that("SDM generated Stan code includes custom SDM chunks", {
  sim <- simulate_sdm_benchmark_data(config = sdm_benchmark_configs()$smoke)
  formula <- bmf(c ~ 0 + condition, kappa ~ 0 + condition)

  code <- stancode(formula, data = sim$data, model = sdm(resp_error = "y"))

  expect_match(code, "sdm_simple_ldenom_chquad_adaptive", fixed = TRUE)
  expect_match(code, "for (g in 1:G_sdm_runs)", fixed = TRUE)
  expect_match(code, "target += sdm_run_count[g] * z", fixed = TRUE)
  expect_false(grepl("c[n] != c[n-1]", code, fixed = TRUE))
  expect_match(code, "COSN", fixed = TRUE)
})
