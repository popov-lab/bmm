source(test_path("..", "..", "inst", "benchmarks", "sdm", "simulate_sdm_data.R"))

test_that("SDM benchmark simulator output works with current SDM workflow", {
  sim <- simulate_sdm_benchmark_data(config = sdm_benchmark_configs()$smoke)
  formula <- bmf(c ~ 0 + condition, kappa ~ 0 + condition)

  sdata <- standata(formula, data = sim$data, model = sdm(resp_error = "y"))

  expect_equal(nrow(sim$data), 100L)
  expect_equal(nrow(sim$true_parameters), 10L)
  expect_equal(sdata$N, nrow(sim$data))
  expect_true(all(sim$data$y >= -pi & sim$data$y <= pi))
})

test_that("SDM generated Stan code includes custom SDM chunks", {
  sim <- simulate_sdm_benchmark_data(config = sdm_benchmark_configs()$smoke)
  formula <- bmf(c ~ 0 + condition, kappa ~ 0 + condition)

  code <- stancode(formula, data = sim$data, model = sdm(resp_error = "y"))

  expect_match(code, "sdm_simple_ldenom_chquad_adaptive", fixed = TRUE)
  expect_match(code, "target += z", fixed = TRUE)
  expect_match(code, "COSN", fixed = TRUE)
})
