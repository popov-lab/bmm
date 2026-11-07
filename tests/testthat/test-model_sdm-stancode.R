simulate_sdm_smoke_data <- function() {
  withr::local_seed(123)
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

  expect_match(code, "real sdm_simple_ldenom(real c, real kappa", fixed = TRUE)
  expect_match(code, "target += sdm_simple_run_ldenom(c, kappa, sdm_gl_x, sdm_gl_w, sdm_wide_cosm1,", fixed = TRUE)
  expect_match(code, "vector[16] sdm_gl_x;", fixed = TRUE)
  expect_match(code, "vector[16] sdm_gl_w;", fixed = TRUE)
  expect_match(code, "vector[32] sdm_wide_cosm1;", fixed = TRUE)
  expect_match(code, "vector[32] sdm_wide_w;", fixed = TRUE)
  expect_false(grepl("c[n] != c[n-1]", code, fixed = TRUE))
  expect_false(grepl("COSN", code, fixed = TRUE))
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
  expect_match(code, "data vector sdm_gl_x, data vector sdm_gl_w, data vector sdm_wide_cosm1, data vector sdm_wide_w", fixed = TRUE)
  expect_match(code, "sdm_simple_run_ldenom_slice(c, kappa, sdm_gl_x, sdm_gl_w,", fixed = TRUE)
  expect_match(code, "sdm_wide_cosm1, sdm_wide_w, start,", fixed = TRUE)
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
  withr::local_seed(1)
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

# Cells whose values differ but whose labels pasted with "." coincide:
# (1, 5.2) and (1.5, 2) both read "1.5.2", as do ("a.b", "c") and ("a", "b.c").
# Every row of a run gets the normalising constant of the run's first row
test_that("SDM runs split cells whose pasted labels coincide", {
  model <- sdm(resp_error = "y")
  runs <- function(data, formula) {
    out <- sdm_run_metadata(data, formula, model)
    c(G = out$G_sdm_runs, start = out$sdm_run_start, count = out$sdm_run_count)
  }
  numeric_formula <- bmf(c ~ a * b, kappa ~ 1)

  numeric_pair <- data.frame(
    y = 0,
    a = c(1, 1, 1, 1, 1.5, 1.5, 1.5, 1.5),
    b = c(2, 2, 5.2, 5.2, 2, 2, 5.2, 5.2)
  )
  expect_equal(runs(numeric_pair, numeric_formula),
               c(G = 4, start1 = 1, start2 = 3, start3 = 5, start4 = 7,
                 count1 = 2, count2 = 2, count3 = 2, count4 = 2))

  factor_pair <- data.frame(
    y = 0,
    f1 = factor(c("a.b", "a.b", "a", "a", "a")),
    f2 = factor(c("c", "c", "b.c", "b.c", "b.c"))
  )
  expect_equal(runs(factor_pair, bmf(c ~ f1 + f2, kappa ~ 1)),
               c(G = 2, start1 = 1, start2 = 3, count1 = 2, count2 = 3))

  # A A B A: a cell that recurs after another starts a run of its own
  recurring <- data.frame(y = 0, a = c(1, 1, 1.5, 1), b = c(5.2, 5.2, 2, 5.2))
  expect_equal(runs(recurring, numeric_formula),
               c(G = 3, start1 = 1, start2 = 3, start3 = 4,
                 count1 = 2, count2 = 1, count3 = 1))

  one_predictor <- data.frame(y = 0, x = factor(c("u", "u", "v", "u", "w", "w")))
  expect_equal(runs(one_predictor, bmf(c ~ x, kappa ~ 1)),
               c(G = 4, start1 = 1, start2 = 3, start3 = 4, start4 = 5,
                 count1 = 2, count2 = 1, count3 = 1, count4 = 2))

  # rows 2 and 5 are dropped, as brms drops them; row 5 sat inside the (1.5, 2) run
  with_na <- data.frame(
    y = c(0, NA, 0, 0, 0, 0, 0),
    a = c(1, 1, 1, 1.5, 1.5, 1.5, 1.5),
    b = c(2, 2, 5.2, 2, NA, 2, 5.2)
  )
  expect_equal(runs(with_na, numeric_formula),
               c(G = 4, start1 = 1, start2 = 2, start3 = 3, start4 = 5,
                 count1 = 1, count2 = 1, count3 = 2, count4 = 1))
})

test_that("standata() of sorted data carries one SDM run per colliding cell", {
  withr::local_options(bmm.sort_data = TRUE)
  dat <- data.frame(
    y = c(-1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5),
    a = c(1.5, 1, 1.5, 1, 1, 1.5, 1, 1.5),
    b = c(2, 5.2, 5.2, 2, 5.2, 2, 2, 5.2)
  )
  sdata <- suppressMessages(standata(bmf(c ~ a * b, kappa ~ a * b),
    data = dat, model = sdm(resp_error = "y")
  ))
  expect_equal(sdata$G_sdm_runs, 4L)
  expect_equal(as.integer(sdata$sdm_run_start), c(1L, 3L, 5L, 7L))
  expect_equal(as.integer(sdata$sdm_run_count), c(2L, 2L, 2L, 2L))
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

test_that("threading(force = TRUE) selects the serial SDM chunk", {
  dat <- simulate_sdm_smoke_data()
  formula <- bmf(c ~ 0 + condition, kappa ~ 0 + condition)

  # brms compiles with threads under force = TRUE but leaves the code unsliced,
  # so the threaded chunk would reference start/end outside partial_log_lik
  withr::local_options(brms.threads = brms::threading(2, force = TRUE))
  code <- stancode(formula, data = dat, model = sdm(resp_error = "y"))

  expect_false(grepl("reduce_sum", code, fixed = TRUE))
  expect_match(code, "+= sdm_simple_run_ldenom(c", fixed = TRUE)
  expect_false(grepl("+= sdm_simple_run_ldenom_slice(", code, fixed = TRUE))
})

test_that("an explicit threads = NULL beats a global threading option", {
  dat <- simulate_sdm_smoke_data()
  formula <- bmf(c ~ 0 + condition, kappa ~ 0 + condition)

  # brms reads threads = NULL as "threading off" and generates serial code, so
  # the sliced chunk would reference start/end outside partial_log_lik
  withr::local_options(brms.threads = brms::threading(2))
  code <- stancode(formula, data = dat, model = sdm(resp_error = "y"), threads = NULL)

  expect_false(grepl("reduce_sum", code, fixed = TRUE))
  expect_match(code, "+= sdm_simple_run_ldenom(c", fixed = TRUE)
  expect_false(grepl("+= sdm_simple_run_ldenom_slice(", code, fixed = TRUE))
})

test_that("SDM Stan code validates run metadata in transformed data", {
  dat <- simulate_sdm_smoke_data()
  formula <- bmf(c ~ 0 + condition, kappa ~ 0 + condition)

  code <- stancode(formula, data = dat, model = sdm(resp_error = "y"))

  expect_match(code, "sdm_run_total != N", fixed = TRUE)
  expect_match(code, "reject(\"bmm error", fixed = TRUE)
})

test_that("update() refreshes SDM run metadata for new data", {
  skip_on_cran() # over 0.4 s in CRAN mode (ADR 0004, R10)
  withr::local_seed(42)
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

test_that("the deprecated sdmSimple() warns and builds the model sdm() builds", {
  expect_warning(old <- sdmSimple("dev_rad"), "`sdmSimple()` is deprecated", fixed = TRUE)
  current <- sdm("dev_rad")
  expect_identical(class(old), class(current))
  expect_identical(old$resp_vars, current$resp_vars)
})
