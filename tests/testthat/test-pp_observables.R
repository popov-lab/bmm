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

two_checks <- list(
  a = .pp_observable(function(d) d$a, label = "A"),
  b = .pp_observable(function(d) d$b, label = "B")
)

test_that(".pp_reduce_na() drops observations whose observed value is NA", {
  observed <- list(a = c(1, 2, NA, 4), b = c(1, 2, 3, 4))
  yrep <- list(a = .pp_expand_data(observed$a, 3L),
               b = .pp_expand_data(observed$b, 3L))
  expect_warning(out <- .pp_reduce_na(two_checks["a"], observed, yrep),
                 "Dropped 1 of 4 observations")
  expect_identical(out$keep, c(TRUE, TRUE, FALSE, TRUE))
  expect_identical(out$values$a$y, c(1, 2, 4))
  expect_identical(dim(out$values$a$yrep), c(3L, 3L))
})

# the defect this replaces: reducing observations for an NA in any draw made
# the retained count decay as (1 - p)^ndraws
test_that(".pp_reduce_na() drops draws, not observations, for NA replicates", {
  observed <- list(a = c(1, 2, 3, 4), b = c(1, 2, 3, 4))
  yrep <- list(a = .pp_expand_data(observed$a, 3L),
               b = .pp_expand_data(observed$b, 3L))
  yrep$a[2L, 2L] <- NA_real_
  expect_warning(out <- .pp_reduce_na(two_checks["a"], observed, yrep),
                 "Dropped 1 of 3 posterior draws")
  expect_true(all(out$keep))
  expect_identical(out$values$a$y, observed$a)
  expect_identical(dim(out$values$a$yrep), c(2L, 4L))
})

test_that(".pp_reduce_na() shares one reduction across all checks", {
  observed <- list(a = c(NA, 2, 3, 4), b = c(1, 2, NA, 4))
  yrep <- list(a = .pp_expand_data(observed$a, 3L),
               b = .pp_expand_data(observed$b, 3L))
  yrep$a[1L, 2L] <- NA_real_
  out <- suppressWarnings(.pp_reduce_na(two_checks, observed, yrep))
  expect_identical(out$keep, c(FALSE, TRUE, FALSE, TRUE))
  expect_identical(out$values$a$y, c(2, 4))
  expect_identical(out$values$b$y, c(2, 4))
  expect_identical(dim(out$values$a$yrep), dim(out$values$b$yrep))
  expect_identical(dim(out$values$b$yrep), c(2L, 2L))
})

test_that(".pp_reduce_na() errors when nothing is left on either dimension", {
  observed <- list(a = rep(NA_real_, 3), b = c(1, 2, 3))
  yrep <- list(a = .pp_expand_data(observed$a, 2L),
               b = .pp_expand_data(observed$b, 2L))
  expect_error(.pp_reduce_na(two_checks["a"], observed, yrep),
               "All observations")

  observed$a <- c(1, 2, 3)
  yrep$a <- .pp_expand_data(observed$a, 2L)
  yrep$a[, 1L] <- NA_real_
  expect_error(suppressWarnings(.pp_reduce_na(two_checks["a"], observed, yrep)),
               "Every posterior draw")
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
test_that(".pp_resolve_type() resolves the check's default and validates it", {
  check <- .pp_observable(function(d) d$x, label = "X", type = "bars")
  expect_identical(.pp_resolve_type(NULL, check, NULL), "bars")
  expect_identical(.pp_resolve_type(NULL, check, "cond"), "bars_grouped")
  expect_identical(.pp_resolve_type("hist", check, NULL), "hist")
  expect_error(.pp_resolve_type("no_such_type", check, NULL), "not a supported")
  expect_error(.pp_resolve_type("loo_pit", check, NULL), "not a supported")
  expect_warning(resolved <- .pp_resolve_type("hist", NULL, NULL), "ignored")
  expect_null(resolved)
})

# the fake fit carries no draws, so prepare_predictions() would fail: reaching
# the type error proves 'type' is validated before anything is simulated
test_that("pp_check() rejects an unknown type before simulating", {
  fit <- fake_bmmfit(ddm(rt = "rt", response = "resp"))
  expect_error(pp_check(fit, resp_var = "rt", type = "no_such_type"),
               "not a supported")
})

