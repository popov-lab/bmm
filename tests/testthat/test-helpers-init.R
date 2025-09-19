test_that("create_initfun returns function for sdm", {
  # prepare info for tests
  ff <- bmmformula(kappa ~ 1, c ~ 1)
  dat <- oberauer_lin_2017
  mod <- sdm(resp_error = "dev_rad")
  config_args <- configure_model(mod, data = dat, formula = ff)

  # create initfun
  init_fun <- create_initfun(mod, dat, config_args$formula)

  # run tests
  expect_equal(class(init_fun),"function")
  expect_equal(class(unlist(init_fun())),"numeric")
})


test_that("create_initfun returns 1 for mixture2p models", {
  # prepare info for tests
  dat <- oberauer_lin_2017
  model_mix2p <- mixture2p(resp_error = "dev_rad")

  ff_mix2p <- bmf(thetat ~ 1, kappa ~ 1)

  config_args_mix2p <- configure_model(model_mix2p, data = dat, formula = ff_mix2p)

  # create initfun
  init_fun <- create_initfun(model_mix2p, dat, config_args_mix2p$formula)

  # run tests
  expect_equal(class(init_fun),"numeric")
  expect_equal(init_fun,1)
})
