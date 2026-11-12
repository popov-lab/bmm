# add_chains() on real fits (#555); tests/testthat/test-add_chains.R mocks
# update()

for (backend in c("cmdstanr", "rstan")) {
  test_that(paste("add_chains() samples new chains next to those of a", backend, "fit"), {
    skip_if_not_installed(backend)
    if (backend == "cmdstanr") {
      skip_if(is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)), "CmdStan is not installed")
    }
    withr::local_seed(555)
    dat <- oberauer_lin_2017[oberauer_lin_2017$ID == 1, ]
    fit <- SW(SM(bmm(
      bmf(thetat ~ 1, kappa ~ 1), dat, mixture2p("dev_rad"),
      backend = backend, chains = 1, iter = 300, warmup = 100, seed = 1, refresh = 0
    )))

    out <- SW(SM(add_chains(fit, refresh = 0)))

    expect_s3_class(out, "bmmfit")
    expect_equal(summary(out)$chains, 2)
    expect_equal(nrow(fit_info(out, "time")), 2)
    draws <- posterior::as_draws_array(out, variable = "b_kappa_Intercept")
    expect_equal(dim(draws)[1:2], c(200, 2))
    expect_false(isTRUE(all.equal(c(draws[, 1, 1]), c(draws[, 2, 1]))))
    expect_equal(c(draws[, 1, 1]), c(posterior::as_draws_array(fit, variable = "b_kappa_Intercept")))
  })
}
