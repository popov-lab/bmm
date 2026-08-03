test_that("default_prior() works with brmsformula", {
  ff <- brms::bf(count ~ zAge + zBase * Trt + (1 | patient))
  prior <- default_prior(ff, data = brms::epilepsy, family = poisson())
  expect_equal(class(prior)[1], "brmsprior")
})

test_that("default_prior() works with formula", {
  ff <- count ~ zAge + zBase * Trt + (1 | patient)
  prior <- default_prior(ff, data = brms::epilepsy, family = poisson())
  expect_equal(class(prior)[1], "brmsprior")
})

test_that("default_prior() works with bmmformula", {
  ff <- bmmformula(kappa ~ 1, thetat ~ 1, thetant ~ 1)
  prior <- default_prior(ff, oberauer_lin_2017, mixture3p(
    resp_error = "dev_rad",
    nt_features = "col_nt",
    set_size = "set_size",
    regex = T
  ))
  expect_equal(class(prior)[1], "brmsprior")
})

test_that("combine prior returns a brmsprior object", {
  prior1 <- brms::prior(normal(0, 1), class = "sd", dpar = "c")
  prior2 <- brms::prior(normal(0, 1), class = "sd", dpar = "kappa")

  # combine the prior
  prior <- combine_prior(prior1, prior2)
  expect_equal(class(prior)[1], "brmsprior")
})

test_that("in combine prior, prior2 overwrites only shared components with prior1", {
  prior1 <- brms::prior(normal(0, 1), class = "sd", dpar = "c1") +
    brms::prior(normal(0, 1), class = "sd", dpar = "c2") +
    brms::prior(normal(0, 1), class = "sd", dpar = "c3")
  prior2 <- brms::prior(normal(0, 1), class = "sd", dpar = "kappa") +
    brms::prior(normal(0, 2), class = "sd", dpar = "c2")

  # combine the prior
  prior <- combine_prior(prior1, prior2)
  expect_equal(nrow(prior), 4)
  expect_equal(dplyr::filter(prior, dpar == "c1"), dplyr::filter(prior1, dpar == "c1"))
  expect_equal(dplyr::filter(prior, dpar == "c2"), dplyr::filter(prior2, dpar == "c2"))
  expect_equal(dplyr::filter(prior, dpar == "c3"), dplyr::filter(prior1, dpar == "c3"))
  expect_equal(dplyr::filter(prior, dpar == "kappa"), dplyr::filter(prior2, dpar == "kappa"))
})

test_that("fixed_pars_priors errors clearly when a fixed parameter is absent from the formula", {
  model <- list(fixed_parameters = list(sdratio = 0))
  formula <- brms::bf(y ~ 1)
  expect_error(
    fixed_pars_priors(model, formula),
    "sdratio.*not part of the model formula"
  )
})

test_that("the fixed-parameter guard fires through the full default_prior pipeline", {
  # a fixed parameter that no configure_model wires into the formula survives
  # check_model/check_formula reconciliation and must be caught before brms
  model <- sdm(resp_error = "y")
  model$fixed_parameters$sdratio <- 0
  dat <- data.frame(y = rsdm(30, mu = 0, c = 3, kappa = 4))
  expect_error(
    default_prior(bmf(c ~ 1, kappa ~ 1), dat, model),
    "sdratio.*not part of the model formula"
  )
})

test_that("fixed_pars_priors builds constant priors for dpar and nlpar fixed parameters", {
  dpar_prior <- fixed_pars_priors(
    list(fixed_parameters = list(sigma = 1)),
    brms::bf(y ~ 1, sigma ~ 1)
  )
  expect_equal(dpar_prior$prior, "constant(1)")
  expect_equal(dpar_prior$dpar, "sigma")

  nlpar_prior <- fixed_pars_priors(
    list(fixed_parameters = list(b = 0)),
    brms::bf(y ~ a + b, a ~ 1, b ~ 1, nl = TRUE)
  )
  expect_equal(nlpar_prior$prior, "constant(0)")
  expect_equal(nlpar_prior$nlpar, "b")
})

test_that("no check for sort_data with default_priors function", {
  withr::local_options("bmm.sort_data" = "check")
  res <- capture_messages(default_prior(
    bmf(kappa ~ set_size, c ~ set_size),
    oberauer_lin_2017,
    sdm("dev_rad")
  ))
  expect_false(any(grepl("sort", res)))
})

prior_row <- function(prior = "", class = "b", coef = "", group = "", resp = "",
                      dpar = "", nlpar = "", source = "user") {
  out <- brms::empty_prior()
  out[1, c("prior", "class", "coef", "group", "resp", "dpar", "nlpar", "source")] <-
    list(prior, class, coef, group, resp, dpar, nlpar, source)
  out
}

test_that("resolve_effective_prior() resolves coef -> group -> class inheritance", {
  prior <- rbind(
    prior_row("student_t(3, 0, 2.5)", class = "sd", dpar = "kappa"),
    prior_row(class = "sd", dpar = "kappa", group = "ID"),
    prior_row(class = "sd", dpar = "kappa", group = "ID", coef = "Intercept"),
    prior_row("normal(0, 1)", class = "b", nlpar = "c"),
    prior_row(class = "b", nlpar = "c", coef = "x1"),
    prior_row(class = "b", dpar = "a"),
    prior_row(class = "b", dpar = "a", coef = "x1")
  )
  expect_equal(
    resolve_effective_prior(prior),
    c(
      "student_t(3, 0, 2.5)", "student_t(3, 0, 2.5)", "student_t(3, 0, 2.5)",
      "normal(0, 1)", "normal(0, 1)", "", ""
    )
  )
})

test_that("classify_priors() distinguishes bmm defaults, brms defaults, user and flat", {
  defaults <- rbind(
    prior_row("normal(0, 1)", class = "b", nlpar = "a"),
    prior_row("student_t(5, 2, 0.75)", class = "Intercept", dpar = "k"),
    prior_row("student_t(3, 0, 2.5)", class = "sd", dpar = "k", source = "default"),
    prior_row(class = "b", dpar = "k", source = "default")
  )
  fit_prior <- rbind(
    prior_row("normal(0, 1)", class = "b", nlpar = "a"),
    prior_row("normal(2, 1)", class = "Intercept", dpar = "k"),
    prior_row("student_t(3, 0, 2.5)", class = "sd", dpar = "k"),
    prior_row(class = "b", dpar = "k", source = "default"),
    prior_row("exponential(1)", class = "sds", dpar = "a")
  )
  res <- classify_priors(fit_prior, defaults, links = list(a = "log", k = "logit"))
  expect_s3_class(res, "data.frame")
  expect_equal(
    res$source,
    c("bmm default", "user", "brms default", "flat", "user")
  )
  expect_equal(res$prior, c("normal(0, 1)", "normal(2, 1)", "student_t(3, 0, 2.5)", "", "exponential(1)"))
})

test_that("classify_priors() drops rows that inherit from a parent row in the table", {
  defaults <- rbind(
    prior_row("normal(0, 1)", class = "b", nlpar = "a"),
    prior_row("student_t(3, 0, 2.5)", class = "sd", nlpar = "a", source = "default")
  )
  fit_prior <- rbind(
    prior_row("normal(0, 1)", class = "b", nlpar = "a"),
    prior_row(class = "b", nlpar = "a", coef = "x1", source = "default"),
    prior_row(class = "b", nlpar = "a", coef = "x2", source = "default"),
    prior_row("student_t(3, 0, 2.5)", class = "sd", nlpar = "a", source = "default"),
    prior_row(class = "sd", nlpar = "a", group = "ID", source = "default"),
    prior_row(class = "sd", nlpar = "a", group = "ID", coef = "Intercept", source = "default"),
    prior_row(class = "b", dpar = "k", source = "default"),
    prior_row(class = "b", dpar = "k", coef = "x1", source = "default")
  )
  res <- classify_priors(fit_prior, defaults, links = list(a = "log"))
  expect_equal(nrow(res), 3)
  expect_equal(res$class, c("b", "sd", "b"))
  expect_equal(res$source, c("bmm default", "brms default", "flat"))
  # coef rows without any class-level parent are kept individually
  orphan <- prior_row("normal(0, 2)", class = "b", coef = "x1", nlpar = "a")
  res_orphan <- classify_priors(orphan, defaults, links = list(a = "log"))
  expect_equal(nrow(res_orphan), 1)
  expect_equal(res_orphan$source, "user")
})

test_that("classify_priors() drops vacuous flat class rows", {
  fit_prior <- rbind(
    prior_row(class = "b", nlpar = "d", source = "default"),
    prior_row("normal(1, 1)", class = "b", nlpar = "d", coef = "Intercept")
  )
  res <- classify_priors(fit_prior, brms::empty_prior(), links = list(d = "log"))
  expect_equal(nrow(res), 1)
  expect_equal(res$coef, "Intercept")
})

test_that("classify_priors() maps prior rows to model parameters and links", {
  fit_prior <- rbind(
    prior_row("normal(0, 1)", class = "b", nlpar = "a"),
    prior_row("normal(0, 1)", class = "Intercept", dpar = "k"),
    prior_row("constant(0)", class = "Intercept"),
    prior_row("student_t(3, 0, 2.5)", class = "sd"),
    prior_row("lkj(1)", class = "cor"),
    prior_row("normal(0, 1)", class = "b", dpar = "z")
  )
  res <- classify_priors(fit_prior, brms::empty_prior(), links = list(a = "log", k = "logit"))
  expect_equal(res$parameter, c("a", "k", "mu", NA, NA, "z"))
  expect_equal(res$link, c("log", "logit", "identity", NA, NA, "identity"))
})

test_that("classify_priors() is robust to schema differences across brms versions", {
  defaults <- prior_row("normal(0, 1)", class = "b", nlpar = "a")
  fit_prior <- prior_row("normal(0, 1)", class = "b", nlpar = "a")
  fit_prior$tag <- NULL
  fit_prior$lb <- ""
  fit_prior$ub <- ""
  res <- classify_priors(fit_prior, defaults, links = list(a = "log"))
  expect_equal(res$source, "bmm default")
})

test_that("report_priors() validates its input", {
  expect_error(report_priors(1), "bmmfit")
})

test_that("report_priors() rejects fits without prior information", {
  skip_on_cran()
  path <- test_path("assets/mock_bmmfit_mixture2p.rds")
  skip_if_not(file.exists(path), "Mock fixture not available (excluded by .Rbuildignore)")
  expect_error(report_priors(readRDS(path)), "no prior information")
})

test_that("report_priors() handles models whose data pipeline transforms the response (m3)", {
  skip_on_cran()
  path <- test_path("assets/bmmfit_m3_ppcheck.rds")
  skip_if_not(file.exists(path), "M3 fixture not available (excluded by .Rbuildignore)")
  fit <- readRDS(path)

  out <- report_priors(fit)
  expect_s3_class(out, "bmm_report_priors")
  # a & c: effects + intercept; d: intercept; b: fixed constant; sd: a, c, d
  expect_equal(nrow(out), 9)
  expect_setequal(out$parameter, c("a", "b", "c", "d"))
  expect_equal(sum(out$source == "bmm default"), 6)
  expect_equal(sum(out$source == "brms default"), 3)
  expect_true(all(out$class[out$source == "brms default"] == "sd"))
  expect_equal(out$prior[out$parameter == "b"], "constant(0.1)")
  expect_false(any(out$source == "flat"))
})

test_that("report_priors() classifies the fixture fit correctly", {
  skip_on_cran()
  path <- test_path("assets/bmmfit_example1.rds")
  skip_if_not(file.exists(path), "SDM fixture not available (excluded by .Rbuildignore)")
  fit <- readRDS(path)

  expect_error(report_priors(fit, format = "foo"), "'arg'")
  expect_silent(out <- report_priors(fit))
  expect_s3_class(out, "bmm_report_priors")
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 3)
  expect_setequal(out$parameter, c("c", "kappa", "mu"))
  expect_equal(out$link[match(c("c", "kappa", "mu"), out$parameter)], c("log", "log", "tan_half"))
  expect_true(all(out$source == "bmm default"))
  expect_equal(out$prior[out$parameter == "c"], "student_t(5, 2, 0.75)")
  expect_equal(out$prior[out$parameter == "kappa"], "student_t(5, 1.75, 0.75)")
  expect_equal(out$prior[out$parameter == "mu"], "constant(0)")
})

test_that("constant_native_value() translates constants to the native scale", {
  expect_equal(constant_native_value("constant(0)", "log"), 1)
  expect_equal(constant_native_value("constant(0)", "logit"), 0.5)
  # unchanged values, non-constants, and non-invertible links carry no information
  expect_true(is.na(constant_native_value("constant(0)", "identity")))
  expect_true(is.na(constant_native_value("constant(0)", "tan_half")))
  expect_true(is.na(constant_native_value("constant(0)", "softmax")))
  expect_true(is.na(constant_native_value("constant(0)", NA_character_)))
  expect_true(is.na(constant_native_value("normal(0, 1)", "log")))
})

test_that("printed reports annotate constants with their native-scale value", {
  skip_on_cran()
  path <- test_path("assets/bmmfit_example1.rds")
  skip_if_not(file.exists(path), "SDM fixture not available (excluded by .Rbuildignore)")
  out <- report_priors(readRDS(path))

  # sampling-scale footnote is always present when links are non-identity
  expect_true(any(grepl("sampling scale", capture.output(print(out)), fixed = TRUE)))

  # log-link constant: table shows the native value, text names both scales
  out$link[out$parameter == "mu"] <- "log"
  printed <- capture.output(print(out))
  expect_true(any(grepl("constant(0) [native: 1]", printed, fixed = TRUE)))
  attr(out, "format") <- "text"
  printed_text <- paste(capture.output(print(out)), collapse = " ")
  expect_match(printed_text, "log \\(sampling\\) scale")
  expect_match(printed_text, "1 on the native scale")
})

test_that("report_priors() omits the technical mu parameter for void_mu models", {
  skip_on_cran()
  path <- test_path("assets/bmmfit_example1.rds")
  skip_if_not(file.exists(path), "SDM fixture not available (excluded by .Rbuildignore)")
  fit <- readRDS(path)

  fit$bmm$model$void_mu <- TRUE
  out <- report_priors(fit)
  expect_setequal(out$parameter, c("c", "kappa"))
})

test_that("report_priors() detects user-modified priors without refitting", {
  skip_on_cran()
  path <- test_path("assets/bmmfit_example1.rds")
  skip_if_not(file.exists(path), "SDM fixture not available (excluded by .Rbuildignore)")
  fit <- readRDS(path)

  kappa_row <- fit$prior$dpar == "kappa" & fit$prior$class == "Intercept"
  fit$prior$prior[kappa_row] <- "normal(0, 1)"
  out <- report_priors(fit)
  expect_equal(out$source[out$parameter == "kappa"], "user")
  expect_true(all(out$source[out$parameter != "kappa"] == "bmm default"))
})

test_that("print.bmm_report_priors() renders table and text formats", {
  skip_on_cran()
  path <- test_path("assets/bmmfit_example1.rds")
  skip_if_not(file.exists(path), "SDM fixture not available (excluded by .Rbuildignore)")
  fit <- readRDS(path)

  out <- report_priors(fit)
  printed <- capture.output(print(out))
  expect_true(any(grepl("bmm default", printed, fixed = TRUE)))
  expect_true(any(grepl("tan_half", printed, fixed = TRUE)))
  expect_false(any(grepl("coef", printed, fixed = TRUE)))

  text_out <- report_priors(fit, format = "text")
  printed_text <- capture.output(print(text_out))
  expect_true(any(grepl("log link", printed_text, fixed = TRUE)))
  expect_true(any(grepl("fixed to 0", printed_text, fixed = TRUE)))
  expect_true(any(grepl("sampling \\(link\\) scale", printed_text)))
})
