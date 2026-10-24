# insight reads a bmm fit through the bmmfit methods of R/easystats-insight.R,
# and parameters, bayestestR and modelbased read it through insight. The
# fixture tests run on every stored fit of the model-cases registry.

fixture_fits <- function() {
  lapply(model_test_cases("fixture"), function(case) load_fixture_fit(case$file))
}

# the population-level predictors of the parameters the user wrote formulas for
user_predictors <- function(fit) {
  uf <- fit$bmm$user_formula
  pars <- names(fit$bmm$model$parameters)
  setdiff(unique(unlist(.population_vars(uf[intersect(names(uf), pars)]))), pars)
}

fit_label <- function(fit) paste(class(fit$bmm$model)[-1], collapse = "/")

technical_parameters <- function(fit) {
  setdiff(names(fit$bmm$model$fixed_parameters), names(fit$bmm$model$parameters))
}

test_that("find_formula() gives one formula per model parameter and the user's predictors", {
  skip_if_not_installed("insight")
  for (fit in fixture_fits()) {
    f <- insight::find_formula(fit)
    label <- fit_label(fit)
    expect_s3_class(f, "insight_formula")
    expect_identical(names(f)[1], "conditional")
    expect_true(all(setdiff(names(f), "conditional") %in%
                      c(names(fit$bmm$model$parameters), paste0(names(fit$bmm$model$parameters), "_random"))),
                info = label)
    expect_setequal(insight::find_predictors(fit, flatten = TRUE) %||% character(0), user_predictors(fit))
    expect_identical(insight::is_nullmodel(fit), length(user_predictors(fit)) == 0, info = label)
    expect_setequal(insight::find_random(fit, flatten = TRUE) %||% character(0),
                    re_group_vars(fit$bmm$user_formula))
    expect_identical(insight::is_mixed_model(fit), length(re_group_vars(fit$bmm$user_formula)) > 0,
                     info = label)
    expect_identical(insight::find_response(fit), unlist(fit$bmm$model$resp_vars, use.names = FALSE))
  }
})

test_that("find_formula() names the predictors on the sdm and m3 fixtures", {
  skip_if_not_installed("insight")
  sdm_fit <- load_fixture_fit("bmmfit_example1.rds")
  expect_identical(insight::find_predictors(sdm_fit, flatten = TRUE), "set_size")
  expect_identical(names(insight::find_formula(sdm_fit)), c("conditional", "mu", "c", "kappa"))

  m3_fit <- load_fixture_fit("bmmfit_m3_ppcheck.rds")
  expect_identical(insight::find_predictors(m3_fit, flatten = TRUE), "cond")
  expect_identical(insight::find_predictors(m3_fit, effects = "random", flatten = TRUE), "ID")
  expect_identical(names(insight::find_random(m3_fit)), c("a_random", "c_random", "d_random"))
})

test_that("find_formula() keeps the group-level terms of a parameter apart", {
  skip_on_cran()
  skip_if_not_installed("insight")
  case <- model_test_cases("sd_default")$sdm
  fit <- suppressWarnings(suppressMessages(
    bmm(case$formula, case$data, case$model, backend = "mock", mock_fit = 1, rename = FALSE)
  ))
  f <- insight::find_formula(fit)
  expect_identical(all.vars(f$kappa), character(0))
  expect_identical(all.vars(f$kappa_random), "ID")
  expect_true(insight::is_mixed_model(fit))
  expect_identical(insight::find_random(fit, flatten = TRUE), "ID")
})
