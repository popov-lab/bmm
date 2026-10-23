load_info_fit <- function(name) restructure(load_fixture_fit(name))

combine_twice <- function(fit) brms::combine_models(fit, fit)

test_that("fit_info() of a single cmdstanr fit reports the stored metadata", {
  fit <- load_info_fit("bmmfit_m3_ppcheck.rds")

  expect_equal(
    fit_info(fit, "time"),
    data.frame(chain_id = 1, warmup = 0.428, sampling = 0.279, total = 0.707)
  )
  expect_equal(
    fit_info(fit, "time_mean"),
    c(warmup = 0.428, sampling = 0.279, total = 0.707)
  )
  expect_length(fit_info(fit, "inv_metric"), 1)
  expect_length(fit_info(fit, "inv_metric")[[1]], 84)
  expect_equal(
    head(fit_info(fit, "inv_metric")[[1]], 3),
    c(0.0298361, 0.0410393, 0.0343983)
  )
  expect_equal(unlist(fit_info(fit, "step_size"), use.names = FALSE), 0.0981616)
})

# the fixture was trimmed to one chain but kept the metadata of all four
test_that("fit_info() reports the chains the stanfit has, not those of leftover metadata", {
  fit <- load_info_fit("bmmfit_legacy_sdm_v0.2.1.rds")

  expect_equal(brms::nchains(fit), 1)
  expect_equal(
    fit_info(fit, "time"),
    data.frame(chain_id = 1, warmup = 6.165, sampling = 7.183, total = 13.348)
  )
  expect_equal(
    fit_info(fit, "inv_metric"),
    list(`1` = c(0.00273441, 0.00153782, 0.00133636, 0.0014172, 0.00173459))
  )
  expect_equal(fit_info(fit, "step_size"), list(`1` = 0.313261))
})

test_that("fit_info() of a single cmdstanr fit takes the stored metadata over the stanfit", {
  fit <- load_info_fit("bmmfit_m3_ppcheck.rds")
  metadata <- attr(fit$fit, "metadata")
  metadata$inv_metric <- list(`1` = c(1, 2))
  metadata$step_size <- list(`1` = 0.5)
  metadata$time$chains$warmup <- 9
  attr(fit$fit, "metadata") <- metadata

  expect_equal(fit_info(fit, "inv_metric"), list(`1` = c(1, 2)))
  expect_equal(fit_info(fit, "step_size"), list(`1` = 0.5))
  expect_equal(fit_info(fit, "time")$warmup, 9)
})

test_that("fit_info() of a single cmdstanr unit_e fit returns NULL for the metric", {
  fit <- load_info_fit("bmmfit_m3_ppcheck.rds")
  attr(fit$fit, "metadata")$inv_metric <- list()

  expect_null(fit_info(fit, "inv_metric"))
})

test_that("fit_info() of a fit without chains returns NULL despite stored metadata", {
  fit <- load_info_fit("bmmfit_m3_ppcheck.rds")
  fit$fit@sim$samples <- list()

  for (what in c("time", "time_mean", "inv_metric", "step_size")) {
    expect_null(fit_info(fit, what))
  }
})

test_that("fit_info() of a combined cmdstanr fit reports every chain", {
  fit <- load_info_fit("bmmfit_m3_ppcheck.rds")
  cmb <- combine_twice(fit)

  expect_equal(brms::nchains(cmb), 2)
  expect_equal(fit_info(cmb, "time")$chain_id, 1:2)
  expect_type(fit_info(cmb, "time")$chain_id, "double")
  expect_equal(fit_info(cmb, "time")[2, -1], fit_info(fit, "time")[1, -1], ignore_attr = TRUE)
  expect_equal(fit_info(cmb, "time_mean")[["total"]], 0.707)
  expect_named(fit_info(cmb, "time_mean"), c("warmup", "sampling", "total"))

  inv_metric <- fit_info(cmb, "inv_metric")
  expect_named(inv_metric, c("1", "2"))
  expect_named(fit_info(cmb, "step_size"), c("1", "2"))
  expect_equal(inv_metric[[2]], fit_info(fit, "inv_metric")[[1]], tolerance = 1e-5)
  expect_equal(
    unlist(fit_info(cmb, "step_size"), use.names = FALSE),
    c(0.0981616, 0.0981616)
  )
})

adaptation_text <- function(metric, step_size) {
  if (is.matrix(metric)) {
    header <- "# Elements of inverse mass matrix:"
    rows <- apply(metric, 1, paste, collapse = ", ")
  } else {
    header <- "# Diagonal elements of inverse mass matrix:"
    rows <- paste(metric, collapse = ", ")
  }
  paste0("# Step size = ", step_size, "\n", header, "\n", paste0("# ", rows, collapse = "\n"), "\n# ")
}

# a second run of the same model with its own timing and adapted metric, the
# values that an update() refit stores
refit_with <- function(fit, metric, step_size, elapsed) {
  chain <- fit$fit@sim$samples[[1]]
  attr(chain, "elapsed_time") <- elapsed
  attr(chain, "adaptation_info") <- adaptation_text(metric, step_size)
  fit$fit@sim$samples[[1]] <- chain
  metadata <- attr(fit$fit, "metadata")
  metadata$inv_metric <- list(`1` = metric)
  metadata$step_size <- list(`1` = step_size)
  metadata$time$chains[, c("warmup", "sampling", "total")] <-
    list(elapsed[[1]], elapsed[[2]], sum(elapsed))
  attr(fit$fit, "metadata") <- metadata
  fit
}

test_that("fit_info() of fits combined one after the other reports each chain's own values", {
  fit1 <- load_info_fit("bmmfit_m3_ppcheck.rds")
  metric2 <- seq(0.5, by = 0.01, length.out = 84)
  metric3 <- seq(2, by = 0.02, length.out = 84)
  fit2 <- refit_with(fit1, metric2, 0.2, c(warmup = 1.5, sample = 0.75))
  fit3 <- refit_with(fit1, metric3, 0.3, c(warmup = 2.5, sample = 1.25))

  cmb <- brms::combine_models(brms::combine_models(fit1, fit2), fit3)

  expect_equal(brms::nchains(cmb), 3)
  expect_equal(fit_info(cmb, "time")$chain_id, 1:3)
  expect_equal(fit_info(cmb, "time")$warmup, c(0.428, 1.5, 2.5))
  expect_equal(fit_info(cmb, "time")$sampling, c(0.279, 0.75, 1.25))
  expect_equal(fit_info(cmb, "time")$total, c(0.707, 2.25, 3.75))
  expect_equal(unlist(fit_info(cmb, "step_size"), use.names = FALSE), c(0.0981616, 0.2, 0.3))
  expect_named(fit_info(cmb, "step_size"), c("1", "2", "3"))
  inv_metric <- fit_info(cmb, "inv_metric")
  expect_named(inv_metric, c("1", "2", "3"))
  expect_equal(inv_metric[[1]], fit_info(fit1, "inv_metric")[[1]], tolerance = 1e-5)
  expect_equal(inv_metric[[2]], metric2, tolerance = 1e-5)
  expect_equal(inv_metric[[3]], metric3, tolerance = 1e-5)
})

test_that("fit_info() of a combined fit reads a dense metric as a matrix", {
  fit1 <- load_info_fit("bmmfit_m3_ppcheck.rds")
  dense <- matrix(c(0.0087, -0.0003, -0.0003, 0.0123), 2)
  fit2 <- refit_with(fit1, dense, 1.09, c(warmup = 0.5, sample = 0.25))

  inv_metric <- fit_info(brms::combine_models(fit1, fit2), "inv_metric")

  expect_length(inv_metric[[1]], 84)
  expect_equal(inv_metric[[2]], dense, tolerance = 1e-5)
})

test_that("fit_info() of a combined fit takes the step size of a unit_e chain from its draws", {
  fit1 <- load_info_fit("bmmfit_m3_ppcheck.rds")
  chain <- fit1$fit@sim$samples[[1]]
  attr(chain, "adaptation_info") <- " "
  n <- length(attr(chain, "sampler_params")$stepsize__)
  attr(chain, "sampler_params")$stepsize__ <- c(rep(1, n - 1), 0.117)
  fit1$fit@sim$samples[[1]] <- chain
  attr(fit1$fit, "metadata") <- NULL

  cmb <- combine_twice(fit1)

  expect_null(fit_info(cmb, "inv_metric"))
  expect_equal(unlist(fit_info(cmb, "step_size"), use.names = FALSE), c(0.117, 0.117))
})

test_that("fit_info() reads the stanfit of an rstan fit", {
  fit <- load_info_fit("bmmfit_example1.rds")

  expect_equal(
    fit_info(fit, "time"),
    data.frame(chain_id = 1, warmup = 1.133, sampling = 1.396, total = 2.529)
  )
  expect_equal(
    fit_info(fit, "inv_metric")[[1]],
    c(0.0160426, 0.00597945, 0.00830777, 0.00713964, 0.0100673)
  )
  expect_equal(unlist(fit_info(fit, "step_size"), use.names = FALSE), 0.369539)

  cmb <- combine_twice(fit)
  expect_equal(fit_info(cmb, "time")$chain_id, 1:2)
  expect_length(fit_info(cmb, "inv_metric"), 2)
})

test_that("fit_info() returns NULL for what the sampler did not record", {
  fit <- load_info_fit("bmmfit_m3_ppcheck.rds")
  chain <- fit$fit@sim$samples[[1]]
  for (name in c("elapsed_time", "adaptation_info", "sampler_params")) {
    attr(chain, name) <- NULL
  }
  fit$fit@sim$samples[[1]] <- chain
  attr(fit$fit, "metadata") <- NULL

  for (what in c("time", "time_mean", "inv_metric", "step_size")) {
    expect_null(fit_info(fit, what))
  }

  fit$fit@sim$samples <- list()
  expect_null(fit_info(fit, "time"))
  expect_null(fit_info(fit, "inv_metric"))
})

test_that("fit_info() returns NULL for a fit without a stanfit", {
  fit <- load_fixture_fit("mock_bmmfit_m3.rds")

  expect_null(fit_info(fit, "time"))
  expect_null(fit_info(fit, "time_mean"))
  expect_null(fit_info(fit, "inv_metric"))
  expect_null(fit_info(fit, "step_size"))
})

test_that("parse_adaptation_info() reads non-finite elements as C++ prints them", {
  info <- paste(
    "# Step size = 0.5",
    "# Diagonal elements of inverse mass matrix:",
    "# NaN, Inf, 0.3, -inf, nan, -nan",
    "# Elapsed Time: 0.1 seconds (Warm-up)",
    sep = "\n"
  )

  expect_equal(
    bmm:::parse_adaptation_info(info),
    list(step_size = 0.5, inv_metric = c(NaN, Inf, 0.3, -Inf, NaN, NaN))
  )
})

test_that("parse_adaptation_info() keeps a dense row with a non-finite element", {
  info <- paste(
    "# Step size = 0.5",
    "# Elements of inverse mass matrix:",
    "# 1, 0.1",
    "# nan, 2",
    "# Elapsed Time: 0.1 seconds (Warm-up)",
    sep = "\n"
  )

  expect_equal(
    bmm:::parse_adaptation_info(info)$inv_metric,
    matrix(c(1, NaN, 0.1, 2), 2)
  )
})

test_that("parse_adaptation_info() returns no metric for a header without rows", {
  info <- "# Step size = 0.5\n# Diagonal elements of inverse mass matrix:\n# "

  expect_equal(
    bmm:::parse_adaptation_info(info),
    list(step_size = 0.5, inv_metric = NULL)
  )
})

test_that("fit_info() of a combined fit returns no time when one chain has none", {
  fit <- load_info_fit("bmmfit_m3_ppcheck.rds")
  untimed <- fit
  attr(untimed$fit@sim$samples[[1]], "elapsed_time") <- NULL

  cmb <- brms::combine_models(untimed, fit)

  expect_null(fit_info(cmb, "time"))
  expect_null(fit_info(cmb, "time_mean"))
})

test_that("fit_info() of a combined fit gives a unit_e chain a NULL metric element", {
  fit <- load_info_fit("bmmfit_m3_ppcheck.rds")
  unit_e <- fit
  attr(unit_e$fit@sim$samples[[1]], "adaptation_info") <- " "

  inv_metric <- fit_info(brms::combine_models(unit_e, fit), "inv_metric")

  expect_named(inv_metric, c("1", "2"))
  expect_null(inv_metric[[1]])
  expect_length(inv_metric[[2]], 84)
})

test_that("parse_adaptation_info() rejects a malformed number without a warning", {
  for (row in c("# 1..2, 3", "# 1e, 3", "# 1e+, 3")) {
    info <- paste(
      "# Step size = 0.5",
      "# Diagonal elements of inverse mass matrix:",
      row,
      sep = "\n"
    )

    expect_no_warning(parsed <- bmm:::parse_adaptation_info(info))
    expect_null(parsed$inv_metric)
  }
})

test_that("parse_adaptation_info() reads the MSVC spellings of NaN and Inf", {
  info <- paste(
    "# Step size = 0.5",
    "# Diagonal elements of inverse mass matrix:",
    "# -nan(ind), infinity, -Infinity, 0.3, .5",
    sep = "\n"
  )

  expect_no_warning(inv_metric <- bmm:::parse_adaptation_info(info)$inv_metric)
  # expect_equal() does not tell NaN from NA
  expect_identical(is.nan(inv_metric), c(TRUE, FALSE, FALSE, FALSE, FALSE))
  expect_equal(inv_metric[-1], c(Inf, -Inf, 0.3, 0.5))
})
