simulate_sdm_smoke_data <- function() {
  set.seed(123)
  condition <- factor(rep(paste0("condition", 1:2), each = 50))
  data.frame(
    y = rsdm(100),
    condition = condition
  )
}

test_that("SDM smoke data works with current SDM workflow", {
  dat <- simulate_sdm_smoke_data()
  formula <- bmf(c ~ 0 + condition, kappa ~ 0 + condition)

  sdata <- standata(formula, data = dat, model = sdm(resp_error = "y"))

  expect_equal(nrow(dat), 100L)
  expect_equal(sdata$N, nrow(dat))
  expect_equal(sdata$G_sdm_runs, 2L)
  expect_equal(as.integer(sdata$sdm_run_count), c(50L, 50L))
  expect_true(all(dat$y >= -pi & dat$y <= pi))
})

test_that("SDM run metadata works with common predictor formulas", {
  dat <- data.frame(
    y = rsdm(12),
    task = factor(rep(c("task1", "task2"), each = 6)),
    condition = factor(rep(rep(c("easy", "hard"), each = 3), times = 2))
  )

  sdata <- standata(bmf(c ~ 1, kappa ~ 1),
    data = dat,
    model = sdm(resp_error = "y")
  )
  expect_equal(sdata$G_sdm_runs, 1L)
  expect_equal(as.integer(sdata$sdm_run_start), 1L)
  expect_equal(as.integer(sdata$sdm_run_count), nrow(dat))
  expect_true(is.array(sdata$sdm_run_start))
  expect_true(is.array(sdata$sdm_run_count))

  sdata <- standata(bmf(c ~ 0 + task, kappa ~ 1),
    data = dat,
    model = sdm(resp_error = "y")
  )
  expect_equal(sdata$G_sdm_runs, 2L)
  expect_equal(as.integer(sdata$sdm_run_count), c(6L, 6L))

  sdata <- standata(bmf(c ~ 0 + task + task:condition, kappa ~ 1),
    data = dat,
    model = sdm(resp_error = "y")
  )
  expect_equal(sdata$G_sdm_runs, 4L)
  expect_equal(as.integer(sdata$sdm_run_count), rep(3L, 4))

  sdata <- standata(bmf(c ~ 0 + task, kappa ~ 0 + condition),
    data = dat,
    model = sdm(resp_error = "y")
  )
  expect_equal(sdata$G_sdm_runs, 4L)
  expect_equal(as.integer(sdata$sdm_run_count), rep(3L, 4))
})

test_that("SDM generated Stan code includes run-level denominator chunks", {
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

test_that("SDM accepts numeric threads shorthand", {
  dat <- simulate_sdm_smoke_data()
  formula <- bmf(c ~ 0 + condition, kappa ~ 0 + condition)

  code <- stancode(
    formula,
    data = dat,
    model = sdm(resp_error = "y"),
    threads = 2
  )

  expect_match(code, "target += reduce_sum", fixed = TRUE)
  expect_match(code, "sdm_simple_run_ldenom_slice", fixed = TRUE)
})
