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

test_that("SDM run metadata aligns with brms after NA rows are dropped", {
  set.seed(1)
  dat <- data.frame(
    y = rsdm(20),
    condition = factor(rep(c("a", "b"), each = 10))
  )
  dat$y[3] <- NA

  sdata <- suppressWarnings(standata(bmf(c ~ 0 + condition, kappa ~ 1),
    data = dat,
    model = sdm(resp_error = "y")
  ))
  expect_equal(sdata$N, 19)
  expect_equal(sum(sdata$sdm_run_count), sdata$N)
  expect_equal(as.integer(sdata$sdm_run_start), c(1L, 10L))
  expect_equal(as.integer(sdata$sdm_run_count), c(9L, 10L))

  # NA in a predictor column (data gets sorted, NA row sorts last and is dropped)
  withr::local_options(bmm.sort_data = TRUE)
  dat2 <- data.frame(
    y = rsdm(20),
    condition = factor(rep(c("a", "b"), each = 10))
  )
  dat2$condition[15] <- NA
  sdata2 <- suppressWarnings(standata(bmf(c ~ 0 + condition, kappa ~ 1),
    data = dat2,
    model = sdm(resp_error = "y")
  ))
  expect_equal(sdata2$N, 19)
  expect_equal(sum(sdata2$sdm_run_count), sdata2$N)
  expect_equal(as.integer(sdata2$sdm_run_start), c(1L, 11L))
  expect_equal(as.integer(sdata2$sdm_run_count), c(10L, 9L))
})

test_that("numeric brms.threads option selects the threaded SDM chunk", {
  dat <- simulate_sdm_smoke_data()
  formula <- bmf(c ~ 0 + condition, kappa ~ 0 + condition)

  # brms accepts a bare number in this option, so bmm must thread the chunk too
  withr::local_options(brms.threads = 2)
  code <- stancode(formula, data = dat, model = sdm(resp_error = "y"))

  expect_match(code, "target += reduce_sum", fixed = TRUE)
  expect_match(code, "+= sdm_simple_run_ldenom_slice(", fixed = TRUE)
  expect_false(grepl("+= sdm_simple_run_ldenom(c", code, fixed = TRUE))
})

test_that("SDM Stan code validates run metadata in transformed data", {
  dat <- simulate_sdm_smoke_data()
  formula <- bmf(c ~ 0 + condition, kappa ~ 0 + condition)

  code <- stancode(formula, data = dat, model = sdm(resp_error = "y"))

  expect_match(code, "sdm_run_total != N", fixed = TRUE)
  expect_match(code, "reject(\"bmm error", fixed = TRUE)
})

test_that("update() refreshes SDM run metadata for new data", {
  set.seed(42)
  dat20 <- data.frame(y = rsdm(20))
  dat40 <- data.frame(y = rsdm(40))

  # a mock stanfit with just enough structure for brms::update.brmsfit to
  # reach its data-revalidation path without any compilation or sampling
  methods::setClass("bmm_mock_stanfit", representation(sim = "list"))
  mockfit <- methods::new("bmm_mock_stanfit", sim = list(
    warmup = 1000, iter = 2000, chains = 1, thin = 1,
    samples = list(structure(list(), args = list(control = list())))
  ))

  fit <- suppressMessages(bmm(bmf(c ~ 1, kappa ~ 1),
    data = dat20, model = sdm(resp_error = "y"),
    backend = "mock", mock_fit = mockfit, rename = FALSE
  ))
  expect_equal(sum(fit$stanvars$sdm_run_count$sdata), 20L)

  fit2 <- suppressMessages(
    update(fit, newdata = dat40, testmode = TRUE, recompile = FALSE)
  )
  expect_equal(sum(fit2$stanvars$sdm_run_count$sdata), 40L)
  expect_equal(fit2$stanvars$G_sdm_runs$sdata, 1L)
})
