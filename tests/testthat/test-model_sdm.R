test_that("sdm() messages that the Laplace cell is not the SDM proper, and only then", {
  expect_message(sdm("y", activation = "laplace"), "is not the signal discrimination model", fixed = TRUE)
  expect_message(sdm("y", activation = "laplace", scaling = "range"), "Table 2", fixed = TRUE)
  for (scaling in c("density", "peak", "range")) {
    expect_no_message(sdm("y", scaling = scaling))
  }
  expect_no_message(sdm("y"))
})

test_that("sdm() takes the scaling of the activation when none is given", {
  expect_identical(sdm("y")$other_vars, list(activation = "vonmises", scaling = "density"))
  expect_identical(
    suppressMessages(sdm("y", activation = "laplace"))$other_vars,
    list(activation = "laplace", scaling = "peak")
  )
  expect_identical(
    suppressMessages(sdm("y", activation = "laplace", scaling = "density"))$other_vars$scaling,
    "density"
  )
  expect_error(sdm("y", activation = "gauss"), "should be one of")
  expect_error(sdm("y", scaling = "height"), "should be one of")
})

test_that("each cell names itself and keeps the class chain", {
  laplace <- suppressMessages(sdm("y", activation = "laplace"))
  expect_identical(laplace$name, "SDM with a Laplace activation function (Oberauer, 2023, Table 2); not the SDM proper")
  expect_identical(sdm("y", scaling = "range")$name, "Signal Discrimination Model (SDM) by Oberauer (2023)")
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

test_that("every scaling keeps the default prior on c", {
  for (scaling in c("density", "peak", "range")) {
    expect_identical(sdm("y", scaling = scaling)$default_priors$c$main, "student_t(5, 2, 0.75)")
  }
  expect_identical(
    suppressMessages(sdm("y", activation = "laplace"))$default_priors$c$main,
    "student_t(5, 2, 0.75)"
  )
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
  laplace <- suppressMessages(sdm("y", activation = "laplace"))
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
    model <- suppressMessages(sdm("y", activation = cell$activation, scaling = cell$scaling))
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
  threaded <- stancode(formula, dat, suppressMessages(sdm("y", activation = "laplace")),
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

# Discrete response options and option sets (#558, K2) -------------------------

grid_prep <- function(Y, dpars, extra = list(), K = 360L, activation = 1L, scaling = 1L) {
  structure(list(
    dpars = dpars,
    data = c(list(Y = Y, sdm_K = K, sdm_activation = activation, sdm_scaling = scaling), extra),
    ndraws = nrow(dpars$c), nobs = length(Y)
  ), class = "brmsprep")
}

test_that("a model with n_options gets the sdm_grid family and the density-scale constant", {
  dat <- data.frame(y = deg2rad(c(0, 10, -20, 5)), cond = c("a", "a", "b", "b"))
  formula <- bmf(c ~ cond, kappa ~ 1)
  for (i in seq_len(nrow(sdm_cells))) {
    cell <- sdm_cells[i, ]
    model <- suppressMessages(sdm("y", activation = cell$activation, scaling = cell$scaling, n_options = 360))
    code <- stancode(formula, dat, model)
    expect_match(code, "target += sdm_grid_lpmf(Y | mu, c, kappa, sdm_K, sdm_activation, sdm_scaling);",
                 fixed = TRUE)
    expect_match(code, "target += -sdm_grid_cells_lse(", fixed = TRUE)
    expect_match(code, "target += -N * log(2 * pi() / sdm_K);", fixed = TRUE)
    expect_match(code, "array[N] int Y;", fixed = TRUE)
    expect_false(grepl("sdm_simple_lpdf|sdm_cont_lpdf", code))
    sdata <- standata(formula, dat, model)
    expect_identical(as.integer(sdata$Y), c(0L, 10L, -20L, 5L))
    expect_identical(as.integer(sdata$sdm_activation), match(cell$activation, c("vonmises", "laplace")))
    expect_identical(as.integer(sdata$sdm_scaling), match(cell$scaling, c("density", "peak", "range")))
  }
  threaded <- stancode(formula, dat, sdm("y", n_options = 360), threads = brms::threading(2))
  expect_match(threaded, "ptarget += sdm_grid_lpmf(Y[start:end] | mu, c, kappa, sdm_K, sdm_activation, sdm_scaling);",
               fixed = TRUE)
  expect_match(threaded, "sdm_grid_cells_lse_slice(", fixed = TRUE)
})

# brms 2.23 passes vint columns to the threaded partial_log_lik() although the
# family does not read them; the probe that the per-trial option columns reach
# the data in both programs
test_that("the option columns reach brms as vint terms in the serial and threaded programs", {
  dat <- data.frame(y = deg2rad(c(0, 90, -90)), o1 = 0, o2 = deg2rad(c(90, 90, -90)))
  model <- sdm("y", n_options = 360, options = c("o1", "o2"))
  formula <- bmf(c ~ 1, kappa ~ 1)
  # two cells in three rows, which messages
  serial <- suppressMessages(stancode(formula, dat, model))
  expect_match(serial, "array[N] int vint1;", fixed = TRUE)
  expect_match(serial, "array[N] int vint2;", fixed = TRUE)
  threaded <- suppressMessages(stancode(formula, dat, model, threads = brms::threading(2)))
  expect_match(threaded, "data array[] int vint1, data array[] int vint2", fixed = TRUE)
  sdata <- suppressMessages(standata(formula, dat, model))
  expect_identical(as.integer(sdata$vint1), c(0L, 0L, 0L))
  expect_identical(as.integer(sdata$vint2), c(90L, 90L, -90L))
})

test_that("log_lik() of the full wheel is dsdm() with n_options, for every cell", {
  draws <- matrix(c(5, 6, 7, 5, 6, 7), 3)
  dpars <- list(mu = matrix(c(0, 0, 0, 0.1, 0.1, 0.1), 3), c = draws, kappa = matrix(2, 3, 2))
  for (a in 1:2) for (s in 1:3) {
    prep <- grid_prep(c(3L, -40L), dpars, activation = a, scaling = s)
    activation <- c("vonmises", "laplace")[a]
    scaling <- c("density", "peak", "range")[s]
    expect_equal(
      log_lik_sdm_grid(2, prep),
      dsdm(deg2rad(-40), 0.1, c(5, 6, 7), 2, n_options = 360, log = TRUE, activation = activation,
           scaling = scaling),
      tolerance = 1e-12
    )
  }
})

test_that("log_lik() of an option set is dsdm() with the trial's options", {
  draws <- matrix(c(5, 6, 7, 5, 6, 7), 3)
  dpars <- list(mu = matrix(0, 3, 2), c = draws, kappa = matrix(2, 3, 2))
  absent <- .option_absent(360L)
  extra <- list(vint1 = c(0L, 0L), vint2 = c(90L, -45L), vint3 = c(absent, 45L))
  prep <- grid_prep(c(90L, 45L), dpars, extra)
  expect_equal(
    log_lik_sdm_grid(1, prep),
    dsdm(deg2rad(90), 0, c(5, 6, 7), 2, n_options = 360, options = deg2rad(c(0, 90)), log = TRUE),
    tolerance = 1e-12
  )
  expect_equal(
    log_lik_sdm_grid(2, prep),
    dsdm(deg2rad(45), 0, c(5, 6, 7), 2, n_options = 360, options = deg2rad(c(0, -45, 45)), log = TRUE),
    tolerance = 1e-12
  )
})

test_that("posterior_predict() draws offered option offsets in steps", {
  withr::local_seed(558)
  n <- 2000
  dpars <- list(mu = matrix(0, n, 2), c = matrix(3, n, 2), kappa = matrix(2, n, 2))
  absent <- .option_absent(360L)
  prep <- grid_prep(c(0L, 120L), dpars, list(vint1 = c(0L, 0L), vint2 = c(120L, 120L),
                                             vint3 = c(-120L, absent)))
  pred <- posterior_predict_sdm_grid(1, prep)
  expect_setequal(unique(pred), c(0L, 120L, -120L))
  expected <- dsdm(0, 0, 3, 2, n_options = 360, options = deg2rad(c(0, 120, -120))) * 2 * pi / 360
  expect_lt(abs(mean(pred == 0) - expected), 0.03)
  expect_setequal(unique(posterior_predict_sdm_grid(2, prep)), c(0L, 120L))
  full <- posterior_predict_sdm_grid(1, grid_prep(c(0L, 1L), dpars, K = 8L))
  expect_true(all(full %in% -3:4))
})

test_that("a draw whose activation is narrower than a step gives a post-fit warning", {
  dat <- data.frame(y = deg2rad(c(0, 1, -2)))
  model <- sdm("y", n_options = 90)
  # w = 1 / sqrt(max(c*, 1) kappa) against the step 2 pi / 90 = 0.07: c* = 50
  # with kappa = 300 gives 0.008, c* = 2 with kappa = 3 gives 0.41
  narrow <- list(mu = matrix(0, 2, 3), c = matrix(c(50, 2), 2, 3), kappa = matrix(c(300, 3), 2, 3))
  prep <- structure(list(dpars = narrow, data = list(sdm_K = 90L, sdm_activation = 1L, sdm_scaling = 2L),
                         ndraws = 2, nobs = 3), class = "brmsprep")
  expect_warning(.sdm_check_grid(model, prep), "In 1 of 2 posterior draws")
  prep$dpars$kappa[] <- 3
  expect_no_warning(.sdm_check_grid(model, prep))
  expect_no_warning(.sdm_check_grid(sdm("y"), prep))
})

test_that("each grid cell runs the bmm() pipeline on the mock backend", {
  dat <- small_data()
  dat$opt1 <- 0
  # the response on the 1-degree grid, and the opposite option where it is the target
  dat$opt2 <- deg2rad(round(rad2deg(dat$dev_rad)))
  dat$opt2[dat$opt2 == 0] <- pi
  dat$cont <- seq(0, 1, length.out = nrow(dat))
  formula <- bmf(c ~ set_size, kappa ~ 1)
  models <- list(
    full = sdm("dev_rad", n_options = 360),
    laplace = suppressMessages(sdm("dev_rad", n_options = 360, activation = "laplace")),
    sets = sdm("dev_rad", n_options = 360, options = c("opt1", "opt2")),
    shorthand = sdm("dev_rad", n_options = 360, n_alternatives = 4)
  )
  dat_short <- dat
  dat_short$dev_rad <- round(dat$dev_rad / (pi / 2)) * pi / 2
  for (name in names(models)) {
    data <- if (name == "shorthand") dat_short else dat
    fit <- suppressMessages(bmm(formula, data, models[[name]], backend = "mock", mock_fit = 1, rename = FALSE))
    expect_s3_class(fit, "bmmfit")
    expect_identical(fit$family$name, "sdm_grid", label = name)
  }
  per_trial <- suppressMessages(bmm(bmf(c ~ cont, kappa ~ 1), dat, models$full,
                                    backend = "mock", mock_fit = 1, rename = FALSE))
  expect_identical(per_trial$family$name, "sdm_grid")
})

test_that("insight::get_data() of a grid fit returns the response column given to bmm()", {
  skip_if_not_installed("insight")
  dat <- data.frame(y = deg2rad(c(0, 10, -20, 5, NA)), cond = c("a", "a", "b", "b", "b"),
                    o1 = 0, o2 = deg2rad(c(90, 10, -20, 5, 90)))
  fit <- suppressWarnings(suppressMessages(
    bmm(bmf(c ~ cond, kappa ~ 1), dat, sdm("y", n_options = 360, options = c("o1", "o2")),
        backend = "mock", mock_fit = 1, rename = FALSE)
  ))
  got <- insight::get_data(fit)
  expect_setequal(names(got), c("y", "cond"))
  expect_equal(got$y, deg2rad(c(0, 10, -20, 5)))
  expect_identical(insight::find_response(fit), "y")
})

test_that("a column named like one check_data() builds is reported before it is overwritten", {
  dat <- data.frame(y = deg2rad(c(0, 90)), o1 = 0, o2 = deg2rad(c(90, 90)), opt_set_id = 7, OptIdx1 = 3)
  expect_warning(
    check_data(sdm("y", n_options = 360, options = c("o1", "o2")), dat, bmf(c ~ 1, kappa ~ 1)),
    "'OptIdx1', 'opt_set_id'"
  )
})
