# report and performance read a bmm fit through the methods of
# R/easystats-report.R, on top of the insight methods. The fixture tests run on
# every stored fit of the model-cases registry; report() is computed once per
# fit, because each call runs loo() and the parameter summaries.

report_fixtures <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      cache <<- lapply(model_test_cases("fixture"), function(case) {
        fit <- load_fixture_fit(case$file)
        list(fit = fit, report = suppressWarnings(suppressMessages(report::report(fit))))
      })
    }
    cache
  }
})

report_fit_label <- function(fit) paste(class(fit$bmm$model)[-1], collapse = "/")

test_that("report() describes each fit once, in bmm's terms", {
  skip_if_not_installed("report")
  skip_if_not_installed("performance")
  for (case in report_fixtures()) {
    label <- report_fit_label(case$fit)
    text <- as.character(case$report)
    expect_s3_class(case$report, "report")
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
  }
})

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
  expect_match(
    suppressMessages(report::report_model(fit)),
    "to predict dev_rad with set_size (formula: mu = 0; c ~ 0 + set_size; kappa ~ 1)",
    fixed = TRUE
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
  for (case in report_fixtures()) {
    expect_identical(
      as.character(report::report_priors(case$fit)),
      prior_report_text(prior_info(case$fit)),
      info = report_fit_label(case$fit)
    )
  }
})

test_that("model_performance() returns only the indices defined for bmm models", {
  skip_if_not_installed("report")
  skip_if_not_installed("performance")
  for (case in report_fixtures()) {
    perf <- attr(attr(case$report, "table"), "performance")
    label <- report_fit_label(case$fit)
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
})
