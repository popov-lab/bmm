dat_vwm <- data.frame(
  id = factor(rep(1:5, each = 4)),
  error = runif(20, -pi, pi)
)
dat_rt <- data.frame(
  id = factor(rep(1:5, each = 4)),
  rt = rlnorm(20, meanlog = -0.5, sdlog = 0.3),
  cond = rep(c("a", "b"), 10)
)

f_vwm <- bmf(thetat ~ 1 + (1 | p | id), kappa ~ 1 + (1 | p | id))
f_rt <- bmf(rt ~ 1 + (1 | p | id))

test_that("bmm_component() validates its arguments", {
  expect_error(
    bmm_component(f_vwm, model = mixture2p(resp_error = "error")),
    "required arguments are missing"
  )
  expect_error(
    bmm_component(f_vwm, data = dat_vwm),
    "Exactly one of 'model' or 'family'"
  )
  expect_error(
    bmm_component(f_vwm,
      model = mixture2p(resp_error = "error"),
      family = brms::lognormal(), data = dat_vwm
    ),
    "Exactly one of 'model' or 'family'"
  )
  expect_error(
    bmm_component(kappa ~ 1, model = mixture2p(resp_error = "error"), data = dat_vwm),
    "must be a bmm formula"
  )
  expect_error(
    bmm_component(f_vwm,
      model = mixture2p(resp_error = "error"), data = dat_vwm,
      prior = "normal(0, 1)"
    ),
    "must be a brmsprior object"
  )
  expect_error(
    bmm_component(f_vwm, model = brms::lognormal(), data = dat_vwm),
    "must be a bmmodel object"
  )
})

test_that("bmm_component() validates family components", {
  expect_error(
    bmm_component(bmf(acc ~ 1), family = brms::bernoulli(), data = dat_rt),
    NA
  )
  expect_error(
    bmm_component(bmf(acc ~ 1), family = binomial(), data = dat_rt),
    "requires response addition terms"
  )
  expect_error(
    bmm_component(bmf(acc ~ 1), family = brms::zero_inflated_binomial(), data = dat_rt),
    "requires response addition terms"
  )
  expect_error(
    bmm_component(f_rt,
      family = brms::mixture("gaussian", "gaussian"),
      data = dat_rt
    ),
    "Mixture and custom families are not supported"
  )
  expect_error(
    bmm_component(f_rt,
      family = brms::custom_family("foo", dpars = "mu", type = "real"),
      data = dat_rt
    ),
    "Mixture and custom families are not supported"
  )
  expect_error(
    bmm_component(f_rt, family = "lognormal", data = dat_rt),
    "must be a distributional family"
  )
  expect_error(
    bmm_component(bmf(sigma ~ 1, rt ~ 1), family = brms::lognormal(), data = dat_rt),
    "first formula of a family component"
  )
})

test_that("family components accept functions and stats families", {
  comp_fun <- bmm_component(f_rt, family = brms::lognormal, data = dat_rt)
  expect_s3_class(comp_fun$model$family, "brmsfamily")
  expect_equal(comp_fun$model$family$family, "lognormal")

  comp_stats <- bmm_component(f_rt, family = gaussian(), data = dat_rt)
  expect_s3_class(comp_stats$model$family, "brmsfamily")
  expect_equal(comp_stats$model$family$family, "gaussian")
})

test_that("the family wrapper model exposes only predicted parameters", {
  comp <- bmm_component(f_rt, family = brms::lognormal(), data = dat_rt)
  expect_s3_class(comp$model, c("bmmodel", "brms_family"), exact = TRUE)
  expect_length(comp$model$parameters, 0)
  expect_equal(comp$model$links$mu, "identity")

  comp_sigma <- bmm_component(
    bmf(rt ~ 1 + (1 | p | id), sigma ~ cond),
    family = brms::lognormal(), data = dat_rt
  )
  expect_named(comp_sigma$model$parameters, "sigma")
  expect_equal(comp_sigma$model$links$sigma, "log")
})

test_that("components combine into a multivariate specification", {
  comp_vwm <- bmm_component(f_vwm, model = mixture2p(resp_error = "error"), data = dat_vwm)
  comp_rt <- bmm_component(f_rt, family = brms::lognormal(), data = dat_rt)

  joint <- comp_vwm + comp_rt
  expect_s3_class(joint, "mvbmmformula")
  expect_length(joint, 2)

  expect_length(joint + comp_vwm, 3)
  expect_length(comp_vwm + joint, 3)
  expect_length(joint + joint, 4)

  expect_error(comp_vwm + 5, "cannot be part of a multivariate bmm model")
  expect_error(comp_vwm + (rt ~ 1), "cannot be part of a multivariate bmm model")
  expect_error(
    comp_vwm + mixture2p(resp_error = "error"),
    "cannot be part of a multivariate bmm model"
  )
})

test_that("family components run through the bmm pipeline", {
  comp <- bmm_component(
    bmf(rt ~ 1 + (1 | p | id), sigma ~ cond),
    family = brms::lognormal(), data = dat_rt
  )
  model <- check_model(comp$model, dat_rt, comp$formula)
  data <- check_data(model, dat_rt, comp$formula)
  formula <- check_formula(model, data, comp$formula)
  cfg <- configure_model(model, data, formula)

  expect_equal(cfg$formula$family$family, "lognormal")
  expect_equal(deparse1(cfg$formula$formula), "rt ~ 1 + (1 | p | id)")
  expect_named(cfg$formula$pforms, "sigma")

  wrong_data <- dat_vwm
  expect_error(
    check_data(model, wrong_data, comp$formula),
    "response variable 'rt' is not present"
  )
  expect_error(
    check_formula(model, data, bmf(rt ~ 1, sigm ~ 1)),
    "Unrecognized model parameters: 'sigm'"
  )
})

test_that("constant distributional parameters work in family components", {
  comp <- bmm_component(
    bmf(rt ~ 1 + (1 | p | id), sigma = 0.3),
    family = brms::lognormal(), data = dat_rt
  )
  model <- check_model(comp$model, dat_rt, comp$formula)
  data <- check_data(model, dat_rt, comp$formula)
  formula <- check_formula(model, data, comp$formula)
  cfg <- configure_model(model, data, formula)
  prior <- configure_prior(model, data, cfg$formula, NULL)

  expect_equal(model$fixed_parameters, list(sigma = 0.3))
  expect_named(cfg$formula$pforms, "sigma")
  constant_row <- prior[prior$dpar == "sigma", ]
  expect_equal(constant_row$prior, "constant(0.3)")
})

test_that("printing components and specifications works", {
  comp_vwm <- bmm_component(f_vwm, model = mixture2p(resp_error = "error"), data = dat_vwm)
  comp_rt <- bmm_component(f_rt, family = brms::lognormal(), data = dat_rt)

  expect_output(print(comp_vwm), "mixture2p")
  expect_output(print(comp_rt), "lognormal")
  expect_output(print(comp_vwm + comp_rt), "2 components")
})
