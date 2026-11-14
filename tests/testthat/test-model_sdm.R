test_that("no cell of sdm() messages when the model is built", {
  for (activation in c("vonmises", "laplace")) {
    for (scaling in c("density", "peak", "range")) {
      expect_no_message(sdm("y", activation = activation, scaling = scaling))
    }
  }
  expect_no_message(sdm("y", activation = "laplace"))
})

test_that("sdm() takes the scaling of the activation when none is given", {
  expect_identical(sdm("y")$other_vars, list(activation = "vonmises", scaling = "density"))
  expect_identical(
    sdm("y", activation = "laplace")$other_vars,
    list(activation = "laplace", scaling = "peak")
  )
  expect_identical(
    sdm("y", activation = "laplace", scaling = "density")$other_vars$scaling,
    "density"
  )
  expect_error(sdm("y", activation = "gauss"), "should be one of")
  expect_error(sdm("y", scaling = "height"), "should be one of")
})

test_that("every cell is named the SDM and keeps the class chain", {
  laplace <- sdm("y", activation = "laplace")
  expect_identical(laplace$name, "Signal Discrimination Model (SDM)")
  expect_identical(sdm("y", scaling = "range")$name, "Signal Discrimination Model (SDM)")
  for (model in list(sdm("y"), sdm("y", scaling = "peak"), laplace)) {
    expect_identical(class(model), c("bmmodel", "circular", "sdm", "sdm_simple"))
  }
})

test_that("printing a model says what c is unless the scaling is density", {
  expect_output(print(sdm("y", scaling = "peak")), "memory strength `c` is the activation of the target.", fixed = TRUE)
  expect_output(print(sdm("y", scaling = "range")), "`c` is the target-opposite contrast in activation.", fixed = TRUE)
  out <- capture.output(print(sdm("y")))
  expect_false(any(grepl("Scaling:", out, fixed = TRUE)))
})

test_that("each activation x scaling cell has its own main priors on c and kappa", {
  main_priors <- function(activation, scaling) {
    priors <- sdm("y", activation = activation, scaling = scaling)$default_priors
    c(c = priors$c$main, kappa = priors$kappa$main)
  }
  expect_identical(main_priors("vonmises", "density"), c(c = "student_t(5, 2, 0.75)", kappa = "student_t(5, 1.75, 0.75)"))
  expect_identical(main_priors("vonmises", "peak"), c(c = "student_t(5, 2, 0.75)", kappa = "student_t(5, 1.75, 0.75)"))
  expect_identical(main_priors("vonmises", "range"), c(c = "student_t(5, 2, 0.75)", kappa = "student_t(5, 1.75, 0.75)"))
  expect_identical(main_priors("laplace", "density"), c(c = "student_t(5, 3.25, 0.75)", kappa = "student_t(5, -0.5, 0.75)"))
  expect_identical(main_priors("laplace", "peak"), c(c = "student_t(5, 2.5, 0.75)", kappa = "student_t(5, -0.5, 0.75)"))
  expect_identical(main_priors("laplace", "range"), c(c = "student_t(5, 2.25, 0.75)", kappa = "student_t(5, -0.5, 0.75)"))
  expect_identical(sdm("y", activation = "laplace")$default_priors$c$main, "student_t(5, 2.5, 0.75)")
  expect_identical(sdm("y")$default_priors$c$effects, "normal(0, 1)")
  expect_identical(sdm("y", activation = "laplace", scaling = "density")$default_priors$kappa$sd, "exponential(1)")
})

test_that("log_lik() and posterior_predict() read the cell from the Stan data", {
  withr::local_seed(558)
  draws <- matrix(c(5, 6, 7, 5, 6, 7), 3)
  prep <- structure(list(
    dpars = list(mu = matrix(0, 3, 2), c = draws, kappa = matrix(2, 3, 2)),
    data = list(Y = c(0.1, -0.4), sdm_activation = 2L, sdm_scaling = 3L),
    ndraws = 3, nobs = 2
  ), class = "brmsprep")
  expect_equal(
    log_lik_sdm_simple(2, prep),
    dsdm(-0.4, 0, c(5, 6, 7), 2, log = TRUE, activation = "laplace", scaling = "range")
  )
  expect_length(posterior_predict_sdm_simple(1, prep), 3)
  # an sdm_simple fit has no codes in its data
  prep$data[c("sdm_activation", "sdm_scaling")] <- NULL
  expect_equal(log_lik_sdm_simple(2, prep), dsdm(-0.4, 0, c(5, 6, 7), 2, log = TRUE))
})

test_that("estimating mu with the Laplace activation warns", {
  dat <- data.frame(y = c(-0.3, 0.1, 0.5, 0.02))
  laplace <- sdm("y", activation = "laplace")
  expect_warning(
    standata(bmf(mu ~ 1, c ~ 1, kappa ~ 1), dat, laplace),
    "kink at every response", fixed = TRUE
  )
  expect_no_warning(standata(bmf(c ~ 1, kappa ~ 1), dat, laplace))
  expect_no_warning(standata(bmf(mu ~ 1, c ~ 1, kappa ~ 1), dat, sdm("y")))
})

test_that("the cells other than the default get the sdm_cont family and pass their codes", {
  dat <- data.frame(y = c(-0.3, 0.1, 0.5, 0.02), cond = c("a", "a", "b", "b"))
  formula <- bmf(c ~ cond, kappa ~ 1)
  for (i in seq_len(nrow(sdm_cells))) {
    cell <- sdm_cells[i, ]
    model <- sdm("y", activation = cell$activation, scaling = cell$scaling)
    sdata <- standata(formula, dat, model)
    code <- stancode(formula, dat, model)
    if (cell$activation == "vonmises" && cell$scaling == "density") {
      expect_null(sdata$sdm_activation)
      expect_false(grepl("sdm_cont", code, fixed = TRUE))
    } else {
      expect_identical(as.integer(sdata$sdm_activation), match(cell$activation, c("vonmises", "laplace")))
      expect_identical(as.integer(sdata$sdm_scaling), match(cell$scaling, c("density", "peak", "range")))
      expect_match(code, "target += sdm_cont_lpdf(Y | mu, c, kappa, sdm_activation, sdm_scaling);", fixed = TRUE)
      expect_match(code, "target += sdm_cont_run_ldenom(", fixed = TRUE)
      expect_identical(as.integer(sdata$G_sdm_runs), 2L)
    }
  }
  threaded <- stancode(formula, dat, sdm("y", activation = "laplace"),
    threads = brms::threading(2)
  )
  expect_match(threaded, "data int sdm_activation, data int sdm_scaling) {", fixed = TRUE)
  expect_match(threaded, "sdm_cont_run_ldenom_slice(c, kappa, sdm_activation, sdm_scaling", fixed = TRUE)
})

test_that("an sdm() fit from before activation and scaling is the default cell", {
  skip_on_cran() # loads a stored fit
  fit <- load_fixture_fit("bmmfit_example1.rds")
  expect_null(fit$bmm$model$other_vars$activation)
  fit <- restructure(fit)
  expect_identical(fit$bmm$model$other_vars, list(activation = "vonmises", scaling = "density"))

  ll <- brms::log_lik(fit, ndraws = 5)
  expect_identical(dim(ll), c(5L, nrow(fit$data)))
  expect_true(all(is.finite(ll)))

  up <- suppressMessages(update(fit, formula. = bmf(c ~ 1, kappa ~ 1), testmode = TRUE))
  expect_s3_class(up, "bmmfit")
  expect_identical(up$bmm$model$other_vars$scaling, "density")
  expect_identical(up$formula$family$name, "sdm_simple")
})

test_that("a fit from 0.4.5 to 0.4.7, whose model class is still sdmSimple, gets the default cell", {
  skip_on_cran() # loads a stored fit
  fit <- load_fixture_fit("bmmfit_example1.rds")
  class(fit$bmm$model) <- c("bmmodel", "circular", "sdmSimple")
  fit$bmm$model$other_vars <- list()
  fit$version$bmm <- package_version("0.4.6")
  fit$version$bmm_restructure <- NULL
  fit <- restructure(fit)
  expect_identical(class(fit$bmm$model), c("bmmodel", "circular", "sdm", "sdm_simple"))
  expect_identical(fit$bmm$model$other_vars, list(activation = "vonmises", scaling = "density"))
  expect_output(print(fit$bmm$model), "Parameters:")
})
