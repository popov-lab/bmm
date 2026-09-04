fake_prep <- function(ndraws, nobs, dpars, data = list()) {
  structure(
    list(ndraws = ndraws, nobs = nobs,
         dpars = lapply(dpars, function(v) {
           if (length(v) == 1L) v else matrix(v, ndraws, nobs)
         }),
         data = data),
    class = "brmsprep"
  )
}

fake_bmmfit <- function(model) {
  structure(list(bmm = list(model = model)), class = "bmmfit")
}

registered_models <- list(
  ddm(rt = "rt", response = "resp"),
  cswald(rt = "rt", response = "resp", version = "simple"),
  cswald(rt = "rt", response = "resp", version = "crisk"),
  ezdm(mean_rt = "mrt", var_rt = "vrt", n_upper = "nu", n_trials = "nt"),
  ezdm(mean_rt = c("mu", "ml"), var_rt = c("vu", "vl"), n_upper = "nu",
       n_trials = "nt", version = "4par")
)

test_that("declared observables name real standata slots", {
  withr::local_seed(1)
  cases <- list(
    list(ddm(rt = "rt", response = "response"),
         rddm(20, drift = 1, bound = 1.5, ndt = 0.3)),
    list(cswald(rt = "rt", response = "response"),
         rcswald(20, drift = 2, bound = 1, ndt = 0.3)),
    list(ezdm(mean_rt = "mean_rt", var_rt = "var_rt", n_upper = "n_upper",
              n_trials = "n_trials"),
         rezdm(5, n_trials = 40, drift = 0.3, bound = 1.2, ndt = 0.3)),
    list(ezdm(mean_rt = c("mean_rt_upper", "mean_rt_lower"),
              var_rt = c("var_rt_upper", "var_rt_lower"),
              n_upper = "n_upper", n_trials = "n_trials", version = "4par"),
         rezdm(5, n_trials = 40, drift = 0.3, bound = 1.2, ndt = 0.3,
               version = "4par"))
  )
  for (case in cases) {
    slots <- names(suppressMessages(standata(bmf(drift ~ 1), case[[2]], case[[1]])))
    expect_true(all(pp_observables(case[[1]])$observed %in% slots))
  }
})

test_that("the default check of every registered spec is the brms Y observable", {
  for (model in registered_models) {
    spec <- pp_observables(model)
    vars <- pp_check_vars(fake_bmmfit(model))
    expect_identical(vars$resp_var[vars$default],
                     names(spec$observed)[spec$observed == "Y"])
  }
})

test_that("signed_rt flips the sign of lower-boundary response times", {
  signed_rt <- .pp_spec_rt_response()$checks$signed_rt$compute
  expect_identical(signed_rt(list(rt = c(0.5, 0.7, 0.4), response = c(1, 0, 1))),
                   c(0.5, -0.7, 0.4))
})

test_that("pp_simulate() errors for a model without a method", {
  model <- structure(list(name = "fakemodel"), class = c("bmmodel", "fakemodel"))
  expect_error(pp_simulate(model, fake_prep(2L, 2L, list())),
               "no pp_simulate")
})

test_that("pp_simulate.ddm() draws each cell from its own parameters", {
  withr::local_seed(1)
  prep <- fake_prep(20L, 6L, dpars = list(
    drift = rep(c(5, -5), each = 60), bound = 1.5, ndt = 0.2, zr = 0.5
  ))
  sims <- pp_simulate(ddm(rt = "rt", response = "resp"), prep)
  expect_identical(dim(sims$rt), c(20L, 6L))
  expect_gt(min(colMeans(sims$response)[1:3]), 0.95)
  expect_lt(max(colMeans(sims$response)[4:6]), 0.05)
  expect_true(all(sims$rt > 0.2))
})

test_that("cswald simple doubles the bound of the two-boundary generator", {
  simple <- withr::with_seed(1, pp_simulate(
    cswald(rt = "rt", response = "r"),
    fake_prep(3L, 4L, dpars = list(drift = 3, bound = 0.8, ndt = 0.2, s = 1))
  ))
  crisk <- withr::with_seed(1, pp_simulate(
    cswald(rt = "rt", response = "r", version = "crisk"),
    fake_prep(3L, 4L, dpars = list(drift = 3, bound = 1.6, ndt = 0.2, zr = 0.5,
                                   s = 1))
  ))
  expect_identical(simple, crisk)
})

test_that("pp_simulate() for ezdm respects per-observation trial counts", {
  model <- ezdm(mean_rt = "mrt", var_rt = "vrt", n_upper = "nu", n_trials = "nt")
  n_trials <- c(10L, 20L, 40L)
  prep <- fake_prep(5L, 3L, dpars = list(
    drift = rep(2, 15), bound = rep(1.2, 15), ndt = rep(0.3, 15), s = 1
  ), data = list(trials = n_trials))
  sims <- pp_simulate(model, prep)
  expect_identical(dim(sims$n_upper), c(5L, 3L))
  for (n in 1:3) {
    expect_true(all(sims$n_upper[, n] <= n_trials[n]))
  }
})

test_that("pp_simulate.ezdm_4par() emits NA where a boundary has < 2 responses", {
  model <- ezdm(mean_rt = c("mu", "ml"), var_rt = c("vu", "vl"), n_upper = "nu",
                n_trials = "nt", version = "4par")
  withr::local_seed(42)
  prep <- fake_prep(30L, 4L, dpars = list(
    drift = rep(0, 120), bound = rep(1, 120), ndt = rep(0.3, 120),
    zr = 0.5, s = 1
  ), data = list(vint2 = rep(3L, 4L)))
  sims <- pp_simulate(model, prep)
  expect_identical(is.na(sims$mean_rt_upper), sims$n_upper < 2)
  expect_identical(is.na(sims$mean_rt_lower), (3L - sims$n_upper) < 2)
})

test_that(".pp_build_ppc_plot() drops NA cells jointly with a warning", {
  model <- ezdm(mean_rt = c("mu", "ml"), var_rt = c("vu", "vl"), n_upper = "nu",
                n_trials = "nt", version = "4par")
  spec <- pp_observables(model)
  prep <- fake_prep(4L, 6L, dpars = list(
    drift = rep(1, 24), bound = rep(1, 24), ndt = rep(0.3, 24), zr = 0.5, s = 1
  ), data = list(Y = rep(0.5, 6), vreal1 = rep(0.5, 6), vreal2 = rep(0.02, 6),
                 vreal3 = rep(0.02, 6), vint1 = rep(10L, 6), vint2 = rep(20L, 6)))
  observed <- lapply(spec$observed, function(slot) prep$data[[slot]])
  yrep_inputs <- lapply(observed, .pp_expand_data, ndraws = prep$ndraws)
  sims <- pp_simulate(model, prep)
  sims$mean_rt_upper[1L, 2L] <- NA_real_
  yrep_inputs[names(spec$observed)] <- sims[names(spec$observed)]
  expect_warning(
    p <- .pp_build_ppc_plot(spec$checks$mean_rt_upper, observed, yrep_inputs,
                            type = NULL, group_vec = NULL, plot_dots = list()),
    "Dropped 1 of 6"
  )
  expect_s3_class(p, "ggplot")

  yrep_inputs$mean_rt_upper[1L, ] <- NA_real_
  expect_error(
    .pp_build_ppc_plot(spec$checks$mean_rt_upper, observed, yrep_inputs,
                       type = NULL, group_vec = NULL, plot_dots = list()),
    "All observations"
  )
})

test_that("pp_check_vars() lists the declared checks", {
  fit <- fake_bmmfit(ddm(rt = "rt", response = "resp"))
  out <- pp_check_vars(fit)
  expect_identical(out$resp_var, c("rt", "response", "signed_rt"))
  expect_identical(out$default, c(TRUE, FALSE, FALSE))
  expect_identical(out$default_type[out$resp_var == "response"], "bars")
})

test_that("pp_check_vars() messages and returns NULL without a declaration", {
  fit <- fake_bmmfit(sdm(resp_error = "y"))
  expect_null(pp_observables(fit$bmm$model))
  expect_message(out <- pp_check_vars(fit), "no additional observables")
  expect_null(out)
})

test_that("pp_check() rejects resp_var for models without a declaration", {
  fit <- fake_bmmfit(sdm(resp_error = "y"))
  expect_error(pp_check(fit, resp_var = "rt"), "declares no additional")
})

test_that("pp_check() rejects an unknown resp_var and lists the options", {
  fit <- fake_bmmfit(ddm(rt = "rt", response = "resp"))
  expect_error(pp_check(fit, resp_var = "accuracy"), "'rt', 'response', 'signed_rt'")
})

test_that("pp_check() rejects negative_rt without a signed_rt observable", {
  fit <- fake_bmmfit(sdm(resp_error = "y"))
  expect_error(pp_check(fit, negative_rt = TRUE), "not supported")
})

test_that("pp_check() rejects negative_rt combined with another resp_var", {
  fit <- fake_bmmfit(ddm(rt = "rt", response = "resp"))
  expect_error(pp_check(fit, resp_var = "rt", negative_rt = TRUE),
               "cannot be combined")
})
