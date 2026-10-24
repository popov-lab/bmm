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
  }
})

test_that("find_response() names the response columns of the data", {
  skip_if_not_installed("insight")
  expected <- list(
    "bmmfit_example1.rds" = "dev_rad",
    "bmmfit_ddm_ppcheck.rds" = c("rt", "response"),
    "bmmfit_cswald_ppcheck.rds" = c("rt", "response"),
    "bmmfit_ezdm3_ppcheck.rds" = c("mean_rt", "var_rt", "n_upper"),
    "bmmfit_ezdm4_ppcheck.rds" = c("mean_rt_upper", "mean_rt_lower", "var_rt_upper", "var_rt_lower", "n_upper"),
    "bmmfit_m3_ppcheck.rds" = c("corr", "other", "dist", "npl"),
    "bmmfit_mpt_ppcheck.rds" = c("old", "new")
  )
  for (file in names(expected)) {
    expect_identical(insight::find_response(load_fixture_fit(file)), expected[[file]], info = file)
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

test_that("find_predictors() returns the user's variables for every effects and component", {
  skip_if_not_installed("insight")
  expected <- list(
    "bmmfit_m3_ppcheck.rds" = list(fixed = "cond", random = "ID", all = c("cond", "ID")),
    "bmmfit_mpt_ppcheck.rds" = list(fixed = NULL, random = "id", all = "id")
  )
  for (file in names(expected)) {
    fit <- load_fixture_fit(file)
    for (component in c("all", "conditional", "location", "distributional", "auxiliary")) {
      for (effects in names(expected[[file]])) {
        expect_identical(
          insight::find_predictors(fit, effects = effects, component = component, flatten = TRUE),
          expected[[file]][[effects]],
          info = paste(file, component, effects)
        )
      }
    }
  }
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

# formulas the registry does not cover: group-level terms that brms expands
# before insight reads them, and non-linear formulas with helper parameters
formula_cases <- function() {
  data <- small_data()
  data$ssn <- as.numeric(data$set_size)
  list(
    nested_groups = list(formula = bmf(c ~ 1 + (1 | ID/session), kappa ~ 1)),
    interaction_group = list(formula = bmf(c ~ 1 + (1 | ID:session), kappa ~ 1)),
    helper_parameters = list(formula = bmf(c ~ c0 + c1 * ssn, c0 ~ 1 + (1 | ID), c1 ~ 1, kappa ~ 1)),
    helper_function = list(
      formula = bmf(c ~ en + (st - en) * exp(-R * ssn), st ~ 1 + (1 | ID), en ~ 1, R ~ 1, kappa ~ 1)
    )
  ) |>
    lapply(function(case) c(case, list(data = data, model = sdm("dev_rad"))))
}

mock_formula_cases <- function() {
  lapply(formula_cases(), function(case) suppressWarnings(suppressMessages(
    bmm(case$formula, case$data, case$model, backend = "mock", mock_fit = 1, rename = FALSE)
  )))
}

# the stored-frame case of every model, the default-prior case of the model
# versions that have no stored-frame case, and the formula cases
draw_free_cases <- local({
  cases <- NULL
  function() {
    if (is.null(cases)) {
      stored <- model_test_cases("stored_frame")
      defaults <- model_test_cases("sd_default")
      version <- function(case) class(case$model)[length(class(case$model))]
      cases <<- c(stored, defaults[!vapply(defaults, version, "") %in% vapply(stored, version, "")],
                  formula_cases())
    }
    cases
  }
})

for (case_name in names(draw_free_cases())) {
  test_that(paste0("the methods that need no draws read the model (", case_name, ")"), {
    skip_on_cran()
    skip_if_not_installed("insight")
    case <- draw_free_cases()[[case_name]]
    fit <- suppressWarnings(suppressMessages(
      bmm(case$formula, case$data, case$model, backend = "mock", mock_fit = 1, rename = FALSE)
    ))
    response <- insight::find_response(fit)
    columns <- names(insight::get_data(fit))
    expect_true(length(response) > 0 && all(response %in% columns))
    expect_true(all(insight::find_predictors(fit, effects = "all", flatten = TRUE) %in% columns))
    expect_identical(length(insight::find_formula(fit)$conditional), 3L)
    expect_no_error(insight::model_info(fit))
    expect_no_error(insight::get_response(fit))
  })
}

test_that("find_random() reads nested and interaction grouping terms as brms does", {
  skip_on_cran()
  skip_if_not_installed("insight")
  fits <- mock_formula_cases()
  expect_identical(insight::find_random(fits$nested_groups), list(c_random = c("ID", "ID:session")))
  expect_identical(insight::find_random(fits$nested_groups, split_nested = TRUE, flatten = TRUE),
                   c("ID", "session"))
  expect_identical(insight::find_random(fits$interaction_group), list(c_random = "ID:session"))
  expect_identical(names(insight::get_data(fits$nested_groups, effects = "random")), c("ID", "session"))
})

test_that("the helper parameters of a non-linear formula are parameters, not predictors", {
  skip_on_cran()
  skip_if_not_installed("insight")
  fits <- mock_formula_cases()
  helpers <- fits$helper_parameters
  expect_identical(names(insight::find_formula(helpers)),
                   c("conditional", "mu", "c", "kappa", "c0", "c1", "c0_random"))
  expect_identical(insight::find_predictors(helpers, flatten = TRUE), "ssn")
  expect_identical(insight::find_predictors(helpers, effects = "all", flatten = TRUE), c("ssn", "ID"))
  expect_identical(insight::find_random(helpers), list(c0_random = "ID"))
  expect_true(insight::is_mixed_model(helpers))
  expect_identical(names(insight::get_data(helpers)), c("dev_rad", "ssn", "ID"))

  function_of_helpers <- fits$helper_function
  expect_identical(insight::find_predictors(function_of_helpers, effects = "all", flatten = TRUE),
                   c("ssn", "ID"))
  expect_identical(insight::find_random(function_of_helpers), list(st_random = "ID"))
})

test_that("find_response() names the count columns of sdt_cdp(), whatever their prefix", {
  skip_on_cran()
  skip_if_not_installed("insight")
  case <- model_test_cases("sd_default")$sdt_cdp
  fit <- suppressWarnings(suppressMessages(
    bmm(case$formula, case$data, case$model, backend = "mock", mock_fit = 1, rename = FALSE)
  ))
  expect_identical(insight::find_response(fit), paste0("cdp", 1:6))
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

test_that("get_priors() reads a brms mixture, whose first component's coefficients carry a prefix", {
  skip_on_cran()
  skip_if_not_installed("insight")
  case <- model_test_cases("stored_frame")$mixture2p
  fit <- suppressWarnings(suppressMessages(
    bmm(case$formula, case$data, case$model, backend = "mock", mock_fit = 1, rename = FALSE)
  ))
  skip_if_not(identical(fit$family$family, "mixture"), "mixture2p is not a brms mixture here")
  # the draws a sampled mixture2p fit of this formula has
  local_mocked_bindings(
    variables = function(x, ...) {
      c("b_mu1_Intercept", "b_mu2_Intercept", "b_kappa2_Intercept", "b_kappa_Intercept",
        "b_thetat_Intercept", "Intercept_mu1", "lprior", "lp__")
    },
    .package = "brms"
  )
  expect_identical(insight::get_priors(fit), data.frame(
    Parameter = c("b_mu1_Intercept", "b_kappa_Intercept", "b_thetat_Intercept"),
    Distribution = c("constant", "normal", "logistic"),
    Location = c(0, 2, 0),
    Scale = c(NA, 1, 1),
    df = NA_real_
  ))
})

test_that("a coefficient of mu whose predictor starts with another parameter's name belongs to mu", {
  skip_on_cran()
  skip_if_not_installed("insight")
  data <- small_data()[1:200, ]
  data$c_load <- as.numeric(data$set_size)
  prior <- brms::set_prior("normal(0, 0.3)", class = "b") +
    brms::set_prior("normal(0, 2)", class = "b", dpar = "c")
  fit <- suppressWarnings(suppressMessages(
    bmm(bmf(mu ~ 1 + c_load, c ~ 1 + c_load, kappa ~ 1), data, sdm("dev_rad"),
        prior = prior, backend = "mock", mock_fit = 1, rename = FALSE)
  ))
  # brms names mu's slope b_c_load and c's slope b_c_c_load
  local_mocked_bindings(
    variables = function(x, ...) {
      c("b_Intercept", "b_c_load", "b_c_Intercept", "b_c_c_load", "b_kappa_Intercept", "lprior", "lp__")
    },
    .package = "brms"
  )
  expect_identical(insight::find_parameters(fit), list(
    mu = c("b_Intercept", "b_c_load"),
    c = c("b_c_Intercept", "b_c_c_load"),
    kappa = "b_kappa_Intercept"
  ))
  priors <- insight::get_priors(fit)
  slopes <- priors[match(c("b_c_load", "b_c_c_load"), priors$Parameter), ]
  expect_identical(slopes$Scale, c(0.3, 2))
})

test_that("get_priors() gives the eta of an lkj prior as its location", {
  skip_on_cran()
  skip_if_not_installed("insight")
  fit <- suppressWarnings(suppressMessages(
    bmm(bmf(kappa ~ 1 + (1 | p | ID), c ~ 1 + (1 | p | ID)), small_data(), sdm("dev_rad"),
        prior = brms::set_prior("lkj(4)", class = "cor", group = "ID"),
        backend = "mock", mock_fit = 1, rename = FALSE)
  ))
  local_mocked_bindings(
    variables = function(x, ...) {
      c("b_c_Intercept", "b_kappa_Intercept", "sd_ID__c_Intercept", "sd_ID__kappa_Intercept",
        "cor_ID__c_Intercept__kappa_Intercept", "lprior", "lp__")
    },
    .package = "brms"
  )
  priors <- insight::get_priors(fit)
  expect_identical(
    as.list(priors[priors$Parameter == "cor_ID__c_Intercept__kappa_Intercept", c("Distribution", "Location", "Scale")]),
    list(Distribution = "lkj", Location = 4, Scale = NA_real_)
  )
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
    expect_identical(unique(mp$Component),
                     intersect(names(fit$bmm$model$parameters), names(fit$formula$pforms)),
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
