# report and performance read a bmm fit through the methods of
# R/easystats-report.R, on top of the insight methods. The fixture tests run on
# every stored fit of the model-cases registry. Each report() call runs loo()
# and the parameter summaries (up to 18 s for the sdm fixture), so each fixture
# gets a test block of its own (R6) and the report is cached per file. The
# model_performance() block reads that cache; testthat runs blocks in file
# order, so it is placed after the per-fixture blocks and finds them filled.

report_cache <- new.env(parent = emptyenv())

fixture_report <- function(file) {
  if (is.null(report_cache[[file]])) {
    fit <- load_fixture_fit(file)
    report_cache[[file]] <- list(
      fit = fit, report = suppressWarnings(suppressMessages(report::report(fit)))
    )
  }
  report_cache[[file]]
}

report_fit_label <- function(fit) paste(class(fit$bmm$model)[-1], collapse = "/")

for (case in model_test_cases("fixture")) {
  test_that(paste0("report() describes the fit in ", case$file, " in bmm's terms"), {
    skip_if_not_installed("report")
    skip_if_not_installed("performance")
    cached <- fixture_report(case$file)
    label <- report_fit_label(cached$fit)
    text <- as.character(cached$report)
    expect_s3_class(cached$report, "report")
    expect_length(text, 1)
    expect_identical(
      lengths(regmatches(text, gregexpr("We fitted a", text, fixed = TRUE))), 1L,
      info = label
    )
    expect_false(grepl("intercept-only", text, fixed = TRUE), info = label)
    expect_false(grepl("general linear", text, fixed = TRUE), info = label)
    expect_false(grepl("custom family", text, fixed = TRUE), info = label)
    expect_false(grepl("\\bb_[A-Za-z]", text), info = label)
    expect_false(grepl("R2", text, fixed = TRUE), info = label)
    expect_false(grepl("..", text, fixed = TRUE), info = label)
  })
}

test_that("format_model() names the bmm model and its version", {
  skip_if_not_installed("report")
  sdm_fit <- load_fixture_fit("bmmfit_example1.rds")
  expect_identical(
    report::format_model(sdm_fit),
    "Bayesian Signal Discrimination Model (SDM) by Oberauer (2023) (simple)"
  )
  ddm_fit <- load_fixture_fit("bmmfit_ddm_ppcheck.rds")
  expect_identical(report::format_model(ddm_fit), "Bayesian Diffusion Decision Model")
  m3_fit <- load_fixture_fit("bmmfit_m3_ppcheck.rds")
  expect_identical(
    report::format_model(m3_fit),
    "Bayesian Multinomial / Memory Measurement Model (custom)"
  )
  mixture2p_fit <- readRDS(test_path("assets", "mock_bmmfit_mixture2p.rds"))
  expect_identical(
    report::format_model(mixture2p_fit),
    "Bayesian Two-parameter mixture model by Zhang and Luck (2008)"
  )
})

test_that("report_model() gives the formula of each model parameter", {
  skip_if_not_installed("report")
  fit <- load_fixture_fit("bmmfit_example1.rds")
  model_report <- suppressMessages(report::report_model(fit))
  expect_match(
    model_report,
    "to predict dev_rad with set_size (formula: mu = 0; c ~ 0 + set_size; kappa ~ 1)",
    fixed = TRUE
  )
  expect_false(grepl("random effect", model_report, fixed = TRUE))
  expect_identical(
    as.character(summary(model_report)),
    "Bayesian Signal Discrimination Model (SDM) by Oberauer (2023) (simple) to predict dev_rad with set_size"
  )
})

test_that("report_model() lists all response variables in one sentence", {
  skip_if_not_installed("report")
  fit <- load_fixture_fit("bmmfit_ezdm4_ppcheck.rds")
  text <- suppressMessages(report::report_model(fit))
  expect_length(text, 1)
  expect_match(
    text,
    "to predict mean_rt_upper, mean_rt_lower, var_rt_upper, var_rt_lower, and n_upper",
    fixed = TRUE
  )
})

test_that("report_model() names the grouping factors of the random effects", {
  skip_if_not_installed("report")
  fit <- load_fixture_fit("bmmfit_m3_ppcheck.rds")
  expect_match(
    suppressMessages(report::report_model(fit)),
    "The model included ID as random effect",
    fixed = TRUE
  )
})

test_that("report_priors() gives the sentences of prior_info()", {
  skip_if_not_installed("report")
  for (case in model_test_cases("fixture")) {
    fit <- load_fixture_fit(case$file)
    expect_identical(
      as.character(report::report_priors(fit)),
      prior_report_text(prior_info(fit)),
      info = report_fit_label(fit)
    )
  }
})

# reads the reports cached by the per-fixture blocks above
test_that("model_performance() returns only the indices defined for bmm models", {
  skip_if_not_installed("report")
  skip_if_not_installed("performance")
  for (case in model_test_cases("fixture")) {
    cached <- fixture_report(case$file)
    perf <- attr(attr(cached$report, "table"), "performance")
    label <- report_fit_label(cached$fit)
    expect_true(all(c("LOOIC", "WAIC") %in% names(perf)), info = label)
    expect_false(any(c("R2", "RMSE") %in% names(perf)), info = label)
  }
  fit <- load_fixture_fit("bmmfit_example1.rds")
  expect_message(
    perf <- performance::model_performance(fit, metrics = c("WAIC", "R2", "RMSE")),
    "R2 and RMSE"
  )
  expect_identical(names(perf), "WAIC")
  expect_message(
    performance::model_performance(fit, metrics = c("ELPD", "R2_adjusted")),
    "^R2_adjusted is not defined"
  )
  expect_no_message(perf <- performance::model_performance(fit, metrics = "all"))
  expect_identical(names(perf), c("ELPD", "ELPD_SE", "LOOIC", "LOOIC_SE", "WAIC"))
  expect_no_message(perf <- performance::model_performance(fit, c("WAIC", "R2"), FALSE))
  expect_identical(names(perf), "WAIC")
})

# the brmsfit method ignores verbose for a sampled fit and report() calls the
# generic without metrics, so what the brmsfit method receives is only visible
# in its arguments
test_that("model_performance() hands metrics and verbose on as given", {
  skip_if_not_installed("report")
  skip_if_not_installed("performance")
  fit <- load_fixture_fit("bmmfit_example1.rds")
  received <- NULL
  local_mocked_bindings(
    model_performance.brmsfit = function(model, metrics = "all", verbose = TRUE, ...) {
      received <<- list(metrics = metrics, verbose = verbose, dots = list(...))
      NULL
    },
    .package = "performance"
  )
  performance::model_performance(fit)
  expect_identical(received, list(metrics = c("LOOIC", "WAIC"), verbose = TRUE, dots = list()))
  performance::model_performance(fit, verbose = FALSE)
  expect_identical(received, list(metrics = c("LOOIC", "WAIC"), verbose = FALSE, dots = list()))
  performance::model_performance(fit, c("WAIC", "R2"), FALSE)
  expect_identical(received, list(metrics = "WAIC", verbose = FALSE, dots = list()))
  expect_message(performance::model_performance(fit, "r2"), "r2 is not defined")
  expect_identical(received, list(metrics = c("LOOIC", "WAIC"), verbose = TRUE, dots = list()))
})
