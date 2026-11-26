# update() is mocked to return the fixture fit as the run of the new chains,
# with the seed it was given, so nothing samples here. The sampling test is in
# tests/local/test-add_chains-fits.R

load_chains_fit <- function(name = "bmmfit_m3_ppcheck.rds") restructure(load_fixture_fit(name))

chain_seeds <- function(fit) {
  vapply(fit$fit@stan_args, function(args) as.character(args$seed), character(1))
}

# the arguments of the last update() call land in `calls`; the new chains store
# the seed passed to update() unless `stored_seed` says otherwise
mock_update <- function(new_fit, stored_seed = NULL, env = parent.frame()) {
  calls <- new.env()
  local_mocked_bindings(
    update = function(object, ...) {
      calls$args <- list(...)
      new_fit$fit@stan_args <- lapply(new_fit$fit@stan_args, function(args) {
        args$seed <- stored_seed %||% calls$args$seed
        args
      })
      new_fit
    },
    .package = "stats",
    .env = env
  )
  calls
}

test_that("add_chains() returns one fit with the chains of the fit and the new ones", {
  fit <- load_chains_fit()
  calls <- mock_update(fit)

  out <- add_chains(fit, chains = 1, seed = 99)

  expect_s3_class(out, "bmmfit")
  expect_equal(brms::nchains(out), 2)
  expect_equal(fit_info(out, "time")$chain_id, c(1, 2))
  expect_equal(chain_seeds(out), c("5", "99"))
})

test_that("add_chains() passes chains, seed and the further arguments to update()", {
  fit <- load_chains_fit()
  calls <- mock_update(fit)

  add_chains(fit, chains = 3, seed = 99, cores = 2, refresh = 0)

  expect_equal(calls$args, list(chains = 3, seed = 99, cores = 2, refresh = 0))
})

test_that("add_chains() without a seed draws one from R's random number generator", {
  fit <- load_chains_fit()
  calls <- mock_update(fit)
  withr::local_seed(555)

  out <- add_chains(fit)

  # sample.int(.Machine$integer.max, 1) after set.seed(555)
  expect_identical(calls$args$seed, 1592386033L)
  expect_equal(chain_seeds(out), c("5", "1592386033"))
})

test_that("add_chains() without a seed draws again when the draw is a seed of the fit", {
  fit <- load_chains_fit()
  fit$fit@stan_args[[1]]$seed <- 1592386033L
  calls <- mock_update(fit)
  withr::local_seed(555)

  out <- add_chains(fit)

  expect_false(identical(as.character(calls$args$seed), "1592386033"))
  expect_equal(brms::nchains(out), 2)
})

test_that("add_chains() refuses the seed of a chain the fit already has", {
  fit <- load_chains_fit()
  calls <- mock_update(fit)

  expect_error(
    add_chains(fit, seed = 5),
    "Chain 1 of 'fit' was sampled with seed 5"
  )
  expect_null(calls$args)

  combined <- add_chains(fit, seed = 99)
  expect_error(add_chains(combined, seed = 99), "Chain 2 of 'fit' was sampled with seed 99")
})

test_that("add_chains() compares seeds as numbers", {
  fit <- load_chains_fit()
  fit$fit@stan_args[[1]]$seed <- 100000L
  calls <- mock_update(fit)

  expect_error(add_chains(fit, seed = 1e5), "Chain 1 of 'fit' was sampled with seed 100000\\.")
  expect_null(calls$args)
})

test_that("add_chains() refuses new chains that ran with a seed of the fit", {
  fit <- load_chains_fit()
  # as cmdstanr, which runs a requested seed of 100000 as seed 1
  calls <- mock_update(fit, stored_seed = 5)

  expect_error(
    add_chains(fit, seed = 100000),
    "New chain 1 ran with seed 5, the seed of chain 1 of 'fit'"
  )
  expect_equal(calls$args$seed, 100000)
})

test_that("add_chains() keeps what describes the model and drops criteria of the old draws", {
  fit <- load_chains_fit()
  fit$criteria$loo <- "computed on the chains of fit"
  calls <- mock_update(fit)

  out <- add_chains(fit, seed = 99)

  expect_identical(out$criteria, list())
  expect_identical(out$bmm, fit$bmm)
  expect_identical(out$formula, fit$formula)
  expect_identical(out$prior, fit$prior)
})

test_that("add_chains() returns a fit without the file of the fit", {
  fit <- load_chains_fit()
  fit$file <- "fit.rds"
  mock_update(fit)

  expect_null(add_chains(fit, seed = 99)$file)
})

test_that("add_chains() refuses new chains from another Stan program", {
  fit <- load_chains_fit()
  new <- fit
  new$model <- sub("model {", "model {\n  // changed", new$model, fixed = TRUE)
  mock_update(new)

  expect_error(add_chains(fit, seed = 99), "different Stan program")
})

test_that("add_chains() accepts a program that differs only in the versions of its header", {
  fit <- load_chains_fit()
  new <- fit
  new$model <- sub("bmm [^\n]+\n", "bmm 9.9.9\n", new$model)
  expect_false(identical(new$model, fit$model))
  mock_update(new)

  expect_equal(brms::nchains(add_chains(fit, seed = 99)), 2)
})

test_that("add_chains() combines rstan fits", {
  fit <- load_chains_fit("bmmfit_example1.rds")
  mock_update(fit)

  out <- add_chains(fit, seed = 99)

  expect_equal(brms::nchains(out), 2)
  expect_equal(chain_seeds(out), c("6", "99"))
  expect_length(fit_info(out, "inv_metric"), 2)
})

test_that("add_chains() refuses what is not a bmm fit sampled with MCMC", {
  fit <- load_chains_fit()
  calls <- mock_update(fit)

  brmsfit <- fit
  class(brmsfit) <- "brmsfit"
  expect_error(add_chains(brmsfit), "'fit' must be a fit returned by bmm()")
  expect_error(add_chains(data.frame()), "'fit' must be a fit returned by bmm()")

  not_mcmc <- "add_chains\\(\\) needs a fit sampled with MCMC"
  expect_error(add_chains(load_fixture_fit("mock_bmmfit_mixture2p.rds")), not_mcmc)
  for (algorithm in c("meanfield", "fixed_param")) {
    variational <- fit
    variational$algorithm <- algorithm
    expect_error(add_chains(variational), not_mcmc, info = algorithm)
  }
  mock_backend <- fit
  mock_backend$backend <- "mock"
  expect_error(add_chains(mock_backend), not_mcmc)
  no_chains <- fit
  no_chains$fit@sim$samples <- list()
  expect_error(add_chains(no_chains), not_mcmc)
  empty <- fit
  empty$fit <- NULL
  expect_error(add_chains(empty), not_mcmc)

  expect_null(calls$args)
})

test_that("add_chains() refuses chains and seeds that are not whole numbers", {
  fit <- load_chains_fit()
  calls <- mock_update(fit)

  for (chains in list(0, -1, 1.5, "2", c(1, 2), NA, NULL)) {
    expect_error(
      add_chains(fit, chains = chains),
      "'chains' must be a positive whole number",
      info = deparse(chains)
    )
  }
  for (seed in list(1.5, "2", c(1, 2), NA)) {
    expect_error(
      add_chains(fit, seed = seed),
      "'seed' must be a whole number or NULL",
      info = deparse(seed)
    )
  }
  expect_null(calls$args)
})

test_that("add_chains() refuses arguments that would change the model, data or draws", {
  fit <- load_chains_fit()
  calls <- mock_update(fit)

  # written out rather than read from the function, so a name dropped from its
  # list fails here
  refused <- list(
    formula. = bmf(c ~ 1), newdata = data.frame(), data = data.frame(), data2 = list(),
    model = sdm("y"), prior = brms::prior(normal(0, 1)), sample_prior = "only",
    iter = 400, warmup = 100, thin = 2, algorithm = "meanfield", backend = "rstan",
    normalize = FALSE, file = "fit.rds", file_compress = FALSE
  )
  for (name in names(refused)) {
    expect_error(
      do.call(add_chains, c(list(fit), refused[name])),
      paste0("add_chains\\(\\) does not take '", gsub(".", "\\.", name, fixed = TRUE), "'"),
      info = name
    )
  }
  # update() and brms::brm() match these to `newdata`, `formula.` and `prior`
  for (name in c("new", "form", "pri")) {
    expect_error(
      do.call(add_chains, c(list(fit), stats::setNames(list(1), name))),
      paste0("add_chains\\(\\) does not take '", name, "'"),
      info = name
    )
  }
  # update() would read it as `formula.`
  expect_error(add_chains(fit, 1, 99, bmf(c ~ 1)), "must be named")

  expect_null(calls$args)
})

test_that("add_chains() passes on arguments whose names start like no refused name", {
  fit <- load_chains_fit()
  calls <- mock_update(fit)
  passed <- list(
    cores = 1, refresh = 0, control = list(adapt_delta = 0.9), init = 0, inits = 0,
    threads = NULL, silent = 2, save_pars = NULL, stan_model_args = list()
  )

  do.call(add_chains, c(list(fit, seed = 99), passed))

  expect_identical(names(calls$args), c("chains", "seed", names(passed)))
})
