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

test_that("model_info() describes a bmm model as nonlinear without a dispersion parameter", {
  skip_if_not_installed("insight")
  for (fit in fixture_fits()) {
    info <- insight::model_info(fit)
    label <- fit_label(fit)
    expect_false(info$is_linear, label = label)
    expect_false(info$is_dispersion, label = label)
    expect_true(info$is_bayesian, label = label)
    expect_identical(info$family, fit$family$family, info = label)
    expect_identical(info$is_mixed, length(re_group_vars(fit$bmm$user_formula)) > 0, info = label)
  }
})

test_that("get_parameters() returns every draw of the model parameters and no technical parameter", {
  skip_if_not_installed("insight")
  for (fit in fixture_fits()) {
    label <- fit_label(fit)
    draws <- insight::get_parameters(fit)
    expect_equal(nrow(draws), brms::ndraws(fit), info = label)
    expect_identical(names(draws), insight::find_parameters(fit, effects = "fixed", flatten = TRUE))
    pars <- insight::find_parameters(fit, effects = "full", flatten = TRUE)
    cleaned <- insight::clean_parameters(fit)
    # the main dpar's intercept is b_Intercept: a parameter only where the model declares mu
    expect_identical("b_Intercept" %in% pars, "mu" %in% names(fit$bmm$model$parameters), info = label)
    expect_identical("b_Intercept" %in% cleaned$Parameter, "mu" %in% names(fit$bmm$model$parameters))
    for (tp in technical_parameters(fit)) {
      expect_false(any(grepl(paste0("^b_", tp, "_"), pars)), info = paste(label, tp))
    }
    expect_true(all(cleaned$Component %in% names(fit$bmm$model$parameters)), info = label)
  }
})

test_that("find_parameters() has one component per model parameter", {
  skip_if_not_installed("insight")
  for (fit in fixture_fits()) {
    expected <- intersect(names(fit$bmm$model$parameters), names(fit$formula$pforms))
    expect_identical(names(insight::find_parameters(fit, effects = "fixed")), expected,
                     info = fit_label(fit))
    expect_identical(names(insight::find_parameters(fit, effects = "fixed", component = "conditional")),
                     expected, info = fit_label(fit))
  }
  expect_identical(names(insight::find_parameters(load_fixture_fit("bmmfit_cswald_ppcheck.rds"))),
                   c("drift", "bound", "ndt", "s"))
})

test_that("find_parameters() honours effects, component, parameters and flatten", {
  skip_if_not_installed("insight")
  fit <- load_fixture_fit("bmmfit_m3_ppcheck.rds")
  all_pars <- insight::find_parameters(fit)
  expect_type(all_pars, "list")
  expect_identical(insight::find_parameters(fit, flatten = TRUE), unlist(all_pars, use.names = FALSE))
  expect_identical(names(insight::find_parameters(fit, effects = "random")),
                   c("a_random", "c_random", "d_random"))
  random <- unlist(insight::find_parameters(fit, effects = "random"), use.names = FALSE)
  expect_true(all(grepl("^(sd|r)_ID__", random)))
  expect_true(any(startsWith(random, "r_")))
  expect_false(any(startsWith(insight::find_parameters(fit, effects = "random_variance", flatten = TRUE), "r_")))
  expect_true(all(startsWith(insight::find_parameters(fit, effects = "grouplevel", flatten = TRUE), "r_")))
  expect_identical(insight::find_parameters(fit, effects = "fixed", component = "c"),
                   list(c = c("b_c_Intercept", "b_c_condoldreordered", "b_c_condoldsame")))
  expect_identical(insight::find_parameters(fit, effects = "fixed", parameters = "Intercept$", flatten = TRUE),
                   c("b_b_Intercept", "b_a_Intercept", "b_c_Intercept", "b_d_Intercept"))
  expect_null(insight::find_parameters(fit, component = "zi", flatten = TRUE))
})

test_that("get_priors() gives the prior of every coefficient, constants with their value", {
  skip_if_not_installed("insight")
  for (fit in fixture_fits()) {
    label <- fit_label(fit)
    priors <- insight::get_priors(fit)
    expect_identical(priors$Parameter, insight::find_parameters(fit, flatten = TRUE), info = label)
    provenance <- prior_provenance(fit)
    expect_setequal(priors$Distribution, sub("\\(.*", "", provenance$prior))
  }

  sdm_priors <- insight::get_priors(load_fixture_fit("bmmfit_example1.rds"))
  expect_identical(sdm_priors, data.frame(
    Parameter = c("b_Intercept", paste0("b_c_set_size", 1:4), "b_kappa_Intercept"),
    Distribution = c("constant", rep("student_t", 5)),
    Location = c(0, rep(2, 4), 1.75),
    Scale = c(NA, rep(0.75, 5)),
    df = c(NA, rep(5, 5))
  ))

  m3_priors <- insight::get_priors(load_fixture_fit("bmmfit_m3_ppcheck.rds"))
  prior_of <- function(variable) {
    as.list(m3_priors[m3_priors$Parameter == variable, c("Distribution", "Location", "Scale")])
  }
  expect_identical(prior_of("b_a_Intercept"), list(Distribution = "normal", Location = 1, Scale = 1))
  expect_identical(prior_of("b_a_condoldsame"), list(Distribution = "normal", Location = 0, Scale = 0.5))
  expect_identical(prior_of("b_b_Intercept"), list(Distribution = "constant", Location = 0.1, Scale = NA_real_))
  expect_identical(prior_of("sd_ID__d_Intercept")$Distribution, "exponential")
})

test_that("get_data() returns the data the user passed, response categories as columns", {
  skip_if_not_installed("insight")
  for (fit in fixture_fits()) {
    label <- fit_label(fit)
    data <- insight::get_data(fit)
    expect_identical(nrow(data), nrow(fit$data), info = label)
    expect_true(all(insight::find_response(fit) %in% names(data)), info = label)
    expect_true(all(insight::find_predictors(fit, effects = "all", flatten = TRUE) %in% names(data)),
                info = label)
    expect_false(any(built_data_columns(fit$bmm$model) %in% names(data)), info = label)
  }

  m3_fit <- load_fixture_fit("bmmfit_m3_ppcheck.rds")
  data <- insight::get_data(m3_fit)
  expect_false("Y" %in% names(data))
  expect_identical(data[c("corr", "other", "dist", "npl")],
                   as.data.frame(m3_fit$data$Y)[c("corr", "other", "dist", "npl")])
  expect_identical(names(insight::get_data(m3_fit, effects = "random")), "ID")
  expect_false("ID" %in% names(insight::get_data(m3_fit, effects = "fixed")))
  expect_identical(names(insight::get_response(m3_fit)), c("corr", "other", "dist", "npl"))
})

test_that("link_function() and link_inverse() use the first parameter's link of every model", {
  skip_if_not_installed("insight")
  for (case in model_test_cases("links")) {
    fit <- structure(list(bmm = list(model = case$model)), class = c("bmmfit", "brmsfit"))
    label <- class(case$model)[length(class(case$model))]
    eta <- insight::link_function(fit)(insight::link_inverse(fit)(0.3))
    expect_equal(eta, 0.3, info = label)
  }
  fit <- structure(list(bmm = list(model = sdm("dev_rad"))), class = c("bmmfit", "brmsfit"))
  expect_identical(insight::link_inverse(fit)(0.3), 2 * atan(0.3))
})

test_that("parameters and bayestestR summarise every model parameter on its own", {
  skip_if_not_installed("insight")
  skip_if_not_installed("parameters")
  skip_if_not_installed("bayestestR")
  for (fit in fixture_fits()) {
    label <- fit_label(fit)
    mp <- suppressWarnings(suppressMessages(parameters::model_parameters(fit)))
    expect_identical(unique(mp$Component), names(insight::find_parameters(fit, effects = "fixed")),
                     info = label)
    dp <- suppressWarnings(suppressMessages(bayestestR::describe_posterior(fit)))
    draws <- insight::get_parameters(fit)
    varying <- names(draws)[vapply(draws, function(d) stats::sd(d) > 0, logical(1))]
    expect_true(all(varying %in% dp$Parameter), info = label)
  }
})

test_that("estimate_means() returns a model parameter on its native scale", {
  skip_if_not_installed("insight")
  skip_if_not_installed("modelbased")
  fit <- load_fixture_fit("bmmfit_example1.rds")
  means <- suppressWarnings(suppressMessages(
    modelbased::estimate_means(fit, by = "set_size", predict = "c")
  ))
  expect_identical(as.character(means$set_size), c("1", "2", "3", "4"))
  expect_identical(fit$bmm$model$links$c, "log")
  draws <- as.data.frame(fit, variable = paste0("b_c_set_size", 1:4))
  expect_equal(means$Median, unname(vapply(draws, function(d) stats::median(exp(d)), numeric(1))),
               tolerance = 1e-6)
})
