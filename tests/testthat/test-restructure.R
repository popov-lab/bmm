test_that("add_links works", {
  x <- .model_mixture2p()
  links <- x$links
  x$links <- NULL
  x <- add_links(x)
  expect_equal(x$links, links)
})

# bmm 0.2.1 saved its fits as plain brmsfit objects; the model and the user's
# formula live only in the environment of the custom family
legacy_fit <- function() {
  fit <- load_fixture_fit("bmmfit_legacy_sdm_v0.2.1.rds")
  class(fit) <- c("bmmfit", class(fit))
  fit
}

test_that("a fit saved by bmm 0.2.1 is rebuilt as the current model", {
  fit <- restructure(legacy_fit())
  expect_identical(class(fit$bmm$model), c("bmmodel", "circular", "sdm", "sdm_simple"))
  expect_equal(fit$bmm$model$links, .model_sdm()$links)
  expect_equal(fit$bmm$user_formula, bmf(c ~ 0 + SetSize, kappa ~ 1), ignore_formula_env = TRUE)
  expect_identical(fit$version$bmm_restructure, utils::packageVersion("bmm"))
  expect_s3_class(summary(fit), "bmmsummary")
})

test_that("an old fit whose family kept no model is refused with the advice to refit", {
  fit <- legacy_fit()
  fit$family$env <- NULL
  expect_error(restructure(fit), "Please refit")
})

# bmm 0.4.1.9000 fitted imm on a brms mixture family of five components, the
# last one guessing; its model object named the set size `setsize`, and the
# version bsc showed only in its class
legacy_imm_fit <- function() load_fixture_fit("bmmfit_legacy_imm_bsc_v0.4.1.rds")

test_that("a fit saved before bmm 0.4.5 on a mixture family is restructured as fitted", {
  fit <- restructure(legacy_imm_fit())
  model <- fit$bmm$model
  expect_identical(class(model), c("bmmodel", "vwm", "non_targets", "imm", "imm_bsc"))
  expect_identical(model$other_vars$set_size, "SetSize")
  expect_identical(model$resp_vars, list(resp_error = "dev_rad"))
  # the links of 0.4.8, which parameterised imm as 0.4.1 did
  expect_identical(model$links, list(mu1 = "tan_half", kappa = "log", c = "identity", s = "log"))
  expect_s3_class(native_parameters(fit, ndraws = 5), "data.frame")
})

test_that("prior_info() of such a fit reports bmm's defaults and drops the guessing component", {
  info <- as.data.frame(prior_info(restructure(legacy_imm_fit())))
  defaults <- info[!nzchar(info$coef), c("parameter", "link", "prior", "source")]
  row.names(defaults) <- NULL
  expect_identical(defaults, data.frame(
    parameter = c("c", "kappa", "s", "mu1"),
    link = c("identity", "log", "log", "tan_half"),
    prior = c("normal(0,1)", "normal(2,1)", "normal(0,1)", "constant(0)"),
    source = "bmm default"
  ))
})

test_that("a fit from 0.4.3 or 0.4.4, whose imm class was not renamed, is recognised by its family", {
  fit <- legacy_imm_fit()
  # the model object as the constructor of 0.4.4.9000 built it
  info <- fit$bmm$model$info
  fit$bmm$model$info <- NULL
  fit$bmm$model[names(info)] <- info
  class(fit$bmm$model) <- c("bmmmodel", "vwm", "nontargets", "IMMspatial", "IMMbsc")
  fit$version$bmm <- as.package_version("0.4.4.9000")
  model <- restructure(fit)$bmm$model
  expect_identical(class(model), c("bmmodel", "vwm", "non_targets", "imm", "imm_bsc"))
  expect_identical(model$other_vars$set_size, "SetSize")
})

test_that("update() refits a fit saved before bmm 0.4.5 under the custom family", {
  # the constant 0.4.1 set on the set-size-1 level of s is dropped (D6 of #434)
  expect_warning(
    fit <- update(legacy_imm_fit(), backend = "mock", mock_fit = 1, rename = FALSE, cores = 1),
    "set size of 1"
  )
  expect_identical(class(fit$bmm$model), c("bmmodel", "circular", "non_targets", "imm", "imm_bsc"))
  expect_identical(fit$bmm$model$other_vars$set_size, "SetSize")
  expect_s3_class(fit$formula$family, "customfamily")
})

test_that("restructure() adds posterior_epred to a custom family saved without it", {
  fit <- load_fixture_fit("bmmfit_example1.rds")
  # as a fit saved before the family carried it
  fit$formula$family$posterior_epred <- NULL
  fit$family$posterior_epred <- NULL
  fit <- restructure(fit)
  expect_identical(fit$formula$family$posterior_epred, posterior_epred_sdm_simple)
  expect_identical(fit$family$posterior_epred, posterior_epred_sdm_simple)
})

test_that("restructure() leaves native brms families without posterior_epred", {
  fit <- load_fixture_fit("bmmfit_m3_ppcheck.rds")
  fit <- restructure(fit)
  expect_null(fit$formula$family$posterior_epred)
})

# a fit saved before #731, whose non-linear formulas reset_env() pointed at the
# global environment
nl_formulas <- function(fit) {
  Filter(function(f) isTRUE(attr(f, "nl")), c(list(fit$formula$formula), fit$formula$pforms))
}

# reset_env() leaves a namespace environment alone, so the old state is set here
with_global_formula_envs <- function(fit) {
  environment(fit$formula$formula) <- globalenv()
  fit$formula$pforms <- lapply(fit$formula$pforms, function(f) {
    environment(f) <- globalenv()
    f
  })
  fit
}

nl_companion_fits <- local({
  fits <- NULL
  function() {
    if (is.null(fits)) fits <<- build_nl_companion_fits()
    fits
  }
})

build_nl_companion_fits <- function() {
  m3_gauss <- m3(c("corr", "other", "npl"), c("n_corr", "n_other", "n_npl"),
                 choice_rule = "gaussian", version = "ss")
  list(
    sdt_ranking = stored_frame_fit(model_test_cases("stored_frame")$sdt_ranking),
    sdt_cdp = stored_frame_fit(model_test_cases("stored_frame")$sdt_cdp),
    sdt_rating = stored_frame_fit(model_test_cases("stored_frame")$sdt_rating),
    m3_gaussian = bmm(bmf(c ~ 1, a ~ 1), oberauer_lewandowsky_2019_e1, m3_gauss,
                      backend = "mock", mock_fit = 1, rename = FALSE)
  )
}

for (kind in names(nl_companion_fits())) {
  test_that(paste("restructure() gives the non-linear formulas of an old", kind, "fit the namespace"), {
    skip_on_cran()
    fit <- nl_companion_fits()[[kind]]
    old <- with_global_formula_envs(fit)
    expect_true(all(vapply(nl_formulas(old), function(f) identical(environment(f), globalenv()), logical(1))))

    repaired <- restructure(old)
    expect_identical(environment(repaired$formula$formula), asNamespace("bmm"))
    expect_nl_functions_resolve(repaired)
  })
}

test_that("restructure() also repairs the formulas of the further categories", {
  skip_on_cran()
  repaired <- restructure(with_global_formula_envs(nl_companion_fits()$sdt_ranking))
  expect_identical(environment(repaired$formula$pforms$murank2), asNamespace("bmm"))
  expect_identical(environment(repaired$formula$pforms$murank5), asNamespace("bmm"))
})

test_that("restructure() repairs the formulas once and leaves a current fit as it is", {
  skip_on_cran()
  fit <- nl_companion_fits()$sdt_ranking
  repaired <- restructure(with_global_formula_envs(fit))
  again <- restructure(repaired)
  expect_identical(lapply(nl_formulas(again), environment), lapply(nl_formulas(repaired), environment))
  current <- restructure(fit)
  expect_identical(lapply(nl_formulas(current), environment), lapply(nl_formulas(fit), environment))
})

test_that("restructure() keeps the environment of a formula that calls no function of bmm", {
  fit <- load_fixture_fit("bmmfit_m3_ppcheck.rds")
  own <- new.env()
  fit$formula$formula <- reset_env(fit$formula$formula, own)
  attr(fit$formula$formula, "nl") <- TRUE
  fit <- restructure(fit)
  expect_identical(environment(fit$formula$formula), own)
})
