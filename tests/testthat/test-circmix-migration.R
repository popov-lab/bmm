# old fits (D1) ---------------------------------------------------------------
# Fits saved by bmm 1.4.1, when mixture2p(), mixture3p() and imm() were brms
# mixture families with non-linear parameters (built by
# tests/local/make-old-circmix-fixtures.R). Post-processing reads them as they
# were fitted; update() refits them under the custom family.

old_circmix_data <- function() {
  dat <- oberauer_lin_2017[oberauer_lin_2017$ID %in% 1:2 &
                             oberauer_lin_2017$set_size %in% c(1, 2, 4), ]
  dat$set_size <- droplevels(dat$set_size)
  dat
}

old_nt_model <- function(name) {
  switch(name,
    mixture3p = mixture3p("dev_rad", nt_features = paste0("col_nt", 1:3), set_size = "set_size"),
    imm = imm("dev_rad", nt_features = paste0("col_nt", 1:3),
              nt_distances = paste0("dist_nt", 1:3), set_size = "set_size", version = "full")
  )
}

# brms::update.brmsfit() reads the sampler settings from fit$fit, where the
# mock fixtures store 1
stanfit_stub <- function() {
  skip_if_not_installed("rstan")
  methods::new("stanfit", sim = list(
    iter = 10L, warmup = 5L, chains = 1L, thin = 1L,
    samples = list(structure(list(), args = list(control = list())))
  ))
}

# the update() result with every message and warning it emitted
refit_old <- function(name, ...) {
  fit <- load_fixture_fit(name)
  if (is.numeric(fit$fit)) fit$fit <- stanfit_stub()
  messages <- character()
  warnings <- character()
  up <- withCallingHandlers(
    update(fit, ..., backend = "mock", mock_fit = stanfit_stub(), rename = FALSE, recompile = TRUE),
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    },
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(fit = up, messages = messages, warnings = warnings)
}

new_mock_fit <- function(formula, model, prior) {
  suppressWarnings(suppressMessages(bmm(
    formula, old_circmix_data(), model, prior = prior,
    backend = "mock", mock_fit = stanfit_stub(), rename = FALSE
  )))
}

lprior_lines <- function(fit) {
  grep("lprior += ", strsplit(brms::stancode(fit), "\n")[[1]], fixed = TRUE, value = TRUE)
}

prior_rows <- function(fit) {
  out <- as.data.frame(fit$prior)[, c("prior", "class", "coef", "group", "dpar", "nlpar")]
  rownames(out) <- NULL
  out
}

has_prior_row <- function(fit, prior, class, coef = "", group = "", dpar = "") {
  p <- prior_rows(fit)
  any(p$prior == prior & p$class == class & p$coef == coef & p$group == group &
        p$dpar == dpar & p$nlpar == "")
}

set_size1_constants <- function(fit) {
  p <- prior_rows(fit)
  sum(grepl("^constant\\(", p$prior) & p$coef == "set_size1")
}

prior_info_rows <- function(fit) {
  out <- as.data.frame(prior_info(fit))
  attributes(out) <- attributes(out)[c("names", "row.names", "class")]
  out
}

expect_restructure_keeps_fit <- function(name) {
  old <- load_fixture_fit(name)
  new <- restructure(old)
  expect_identical(new$bmm$model, old$bmm$model)
  expect_identical(new$formula, old$formula)
  expect_identical(new$prior, old$prior)
}

test_that("restructure() leaves a 1.4 mixture2p fit as it was fitted", {
  expect_restructure_keeps_fit("bmmfit_v141_mixture2p.rds")
  expect_identical(class(load_fixture_fit("bmmfit_v141_mixture2p.rds")$bmm$model),
                   c("bmmodel", "circular", "mixture2p"))
})

test_that("restructure() leaves a 1.4 mixture3p fit as it was fitted", {
  expect_restructure_keeps_fit("mock_v141_mixture3p.rds")
})

test_that("restructure() leaves a 1.4 imm fit as it was fitted", {
  expect_restructure_keeps_fit("mock_v141_imm.rds")
})

test_that("prior_info() reports a 1.4 mixture2p fit as fitted", {
  old <- load_fixture_fit("bmmfit_v141_mixture2p.rds")
  # check_data() of a 1.4 model, which has no version, once failed here and
  # prior_info() warned that the data could not be checked again
  expect_no_warning(rows <- prior_info_rows(old))
  expect_identical(rows, data.frame(
    parameter = c("kappa", "kappa", "thetat", "mu1", "kappa"),
    link = c("log", "log", "logit", "tan_half", "log"),
    class = c("b", "b", "b", "Intercept", "sd"),
    coef = c("", "Intercept", "Intercept", "", ""),
    group = c("", "", "", "", ""),
    prior = c("normal(0, 0.3)", "normal(1, 0.5)", "logistic(0, 1)", "constant(0)", "exponential(3)"),
    source = c("user", "user", "bmm default", "bmm default", "user")
  ))
})

test_that("prior_info() reports a 1.4 mixture3p fit as fitted", {
  old <- load_fixture_fit("mock_v141_mixture3p.rds")
  expect_no_warning(rows <- prior_info_rows(old))
  expect_identical(rows, data.frame(
    parameter = c("kappa", "thetant", "thetant", "thetat", "mu1", "thetant", "thetant"),
    link = c("log", "softmax", "softmax", "softmax", "tan_half", "softmax", "softmax"),
    class = c("b", "b", "b", "b", "Intercept", "sd", "sd"),
    coef = c("", "", "set_size1", "", "", "", "set_size1"),
    group = c("", "", "", "", "", "", "ID"),
    prior = c("normal(2, 0.5)", "logistic(0, 1)", "constant(-100)", "normal(0, 1)",
              "student_t(1, 0, 0.3)", "exponential(1)", "constant(1e-8)"),
    source = c("user", "bmm default", "user", "user", "user", "bmm default", "user")
  ))
})

test_that("prior_info() reports a 1.4 imm fit as fitted", {
  old <- load_fixture_fit("mock_v141_imm.rds")
  expect_no_warning(rows <- prior_info_rows(old))
  expect_identical(rows, data.frame(
    parameter = c("a", "a", "c", "kappa", "s", "mu1", "a", "a"),
    link = c("log", "log", "log", "log", "log", "tan_half", "log", "log"),
    class = c("b", "b", "b", "b", "b", "Intercept", "sd", "sd"),
    coef = c("", "set_size1", "", "Intercept", "Intercept", "", "", "set_size1"),
    group = c("", "", "", "", "", "", "", "ID"),
    prior = c("normal(0, 1)", "constant(0)", "normal(0, 0.7)", "normal(2, 1)", "normal(0, 1)",
              "constant(0)", "exponential(2)", "constant(1e-8)"),
    source = c("bmm default", "user", "user", "bmm default", "bmm default", "bmm default",
               "user", "user")
  ))
})

test_that("emmeans() routes a parameter of a 1.4 fit to its nlpar, without a warning", {
  for (name in c("bmmfit_v141_mixture2p.rds", "mock_v141_mixture3p.rds", "mock_v141_imm.rds")) {
    old <- load_fixture_fit(name)
    expect_no_warning(by_dpar <- .bmmfit_resolve_par(old, dpar = "kappa", nlpar = NULL))
    expect_identical(by_dpar, list(dpar = NULL, nlpar = "kappa"), label = name)
    expect_no_warning(by_nlpar <- .bmmfit_resolve_par(old, dpar = NULL, nlpar = "kappa"))
    expect_identical(by_nlpar, list(dpar = NULL, nlpar = "kappa"), label = name)
  }
})

test_that("emmeans() on a 1.4 mixture2p fit gives the same means by dpar and by nlpar", {
  skip_on_cran() # 0.40 s in CRAN mode (ADR 0004, R10)
  skip_if_not_installed("emmeans")
  old <- load_fixture_fit("bmmfit_v141_mixture2p.rds")
  expect_no_warning(by_dpar <- as.data.frame(emmeans::emmeans(old, ~set_size, dpar = "kappa")))
  expect_no_warning(by_nlpar <- as.data.frame(emmeans::emmeans(old, ~set_size, nlpar = "kappa")))
  expect_equal(by_dpar$emmean, c(1.74496275, 1.15032743, 0.28118060), tolerance = 1e-6)
  expect_equal(by_dpar$lower.HPD, c(0.64330591, 0.37250224, -0.62541227), tolerance = 1e-6)
  expect_identical(by_nlpar, by_dpar)
})

test_that("native_parameters() reads a 1.4 mixture2p fit under its fitted names", {
  skip_on_cran() # 0.41 s in CRAN mode (ADR 0004, R10)
  old <- load_fixture_fit("bmmfit_v141_mixture2p.rds")
  expect_no_warning(np <- native_parameters(old, re_formula = NA, summary = TRUE))
  expect_identical(np$parameter, rep(c("mu1", "kappa", "thetat"), each = 3))
  expect_identical(as.character(np$set_size), rep(c("1", "2", "4"), 3))
  expect_equal(np$Estimate, c(
    0, 0, 0,
    6.413662836189575, 3.550133240767290, 1.537117488563647,
    0.896159810841888, 0.896159810841888, 0.896159810841888
  ), tolerance = 1e-10)
})

test_that("update() refits a 1.4 mixture2p fit under the custom family", {
  skip_on_cran()
  res <- refit_old("bmmfit_v141_mixture2p.rds")
  up <- res$fit
  expect_length(grep("custom family", res$messages), 1)
  expect_match(grep("custom family", res$messages, value = TRUE), "dpar =", fixed = TRUE)
  expect_length(grep("deprecated", res$warnings), 0)
  expect_identical(class(up$bmm$model), c("bmmodel", "circular", "mixture2p", "mixture2p_simple"))
  expect_identical(up$bmm$model$dpar_since, "1.5.0")
  expect_s3_class(up$formula$family, "customfamily")
  expect_identical(unique(prior_rows(up)$nlpar), "")
  expect_true(has_prior_row(up, "normal(1, 0.5)", "Intercept", dpar = "kappa"))
  expect_true(has_prior_row(up, "normal(0, 0.3)", "b", dpar = "kappa"))
  expect_true(has_prior_row(up, "exponential(3)", "sd", dpar = "kappa"))
  expect_length(lprior_lines(up), 4)
  expect_identical(lprior_lines(up), lprior_lines(new_mock_fit(
    bmf(kappa ~ 1 + set_size + (1 | ID), thetat ~ 1), mixture2p("dev_rad"),
    brms::set_prior("normal(1, 0.5)", class = "Intercept", dpar = "kappa") +
      brms::set_prior("normal(0, 0.3)", class = "b", dpar = "kappa") +
      brms::set_prior("exponential(3)", class = "sd", dpar = "kappa")
  )))
})

test_that("update() refits a 1.4 mixture3p fit under the custom family", {
  skip_on_cran()
  res <- refit_old("mock_v141_mixture3p.rds")
  up <- res$fit
  expect_length(grep("custom family", res$messages), 1)
  expect_match(grep("custom family", res$messages, value = TRUE), "dpar =", fixed = TRUE)
  expect_match(grep("custom family", res$messages, value = TRUE),
               "set-size-1 level of 'thetant' are dropped", fixed = TRUE)
  expect_length(grep("deprecated", res$warnings), 0)
  expect_identical(class(up$bmm$model),
                   c("bmmodel", "circular", "non_targets", "mixture3p", "mixture3p_simple"))
  expect_identical(up$bmm$model$dpar_since, "1.5.0")
  expect_s3_class(up$formula$family, "customfamily")
  expect_identical(names(up$bmm$user_formula), c("mu", "kappa", "thetat", "thetant"))
  expect_identical(unique(prior_rows(up)$nlpar), "")
  # mu1 is the main dpar; a class-b row on kappa (intercept and slope) splits
  # in two; thetat has no intercept, so its class-b row stays one
  expect_true(has_prior_row(up, "student_t(1, 0, 0.3)", "Intercept"))
  expect_true(has_prior_row(up, "normal(2, 0.5)", "Intercept", dpar = "kappa"))
  expect_true(has_prior_row(up, "normal(2, 0.5)", "b", dpar = "kappa"))
  expect_true(has_prior_row(up, "normal(0, 1)", "b", dpar = "thetat"))
  expect_false(has_prior_row(up, "normal(0, 1)", "Intercept", dpar = "thetat"))
  expect_identical(set_size1_constants(up), 0L)
  expect_length(lprior_lines(up), 6)
  expect_identical(lprior_lines(up), lprior_lines(new_mock_fit(
    bmf(mu ~ 1, kappa ~ 1 + session, thetat ~ 0 + set_size,
        thetant ~ 0 + set_size + (0 + set_size || ID)),
    old_nt_model("mixture3p"),
    brms::set_prior("student_t(1, 0, 0.3)", class = "Intercept") +
      brms::set_prior("normal(2, 0.5)", class = "Intercept", dpar = "kappa") +
      brms::set_prior("normal(2, 0.5)", class = "b", dpar = "kappa") +
      brms::set_prior("normal(0, 1)", class = "b", dpar = "thetat")
  )))
})

test_that("update() refits a 1.4 imm fit under the custom family", {
  skip_on_cran()
  res <- refit_old("mock_v141_imm.rds")
  up <- res$fit
  expect_length(grep("custom family", res$messages), 1)
  expect_match(grep("custom family", res$messages, value = TRUE), "dpar =", fixed = TRUE)
  expect_match(grep("custom family", res$messages, value = TRUE),
               "set-size-1 level of 'a' are dropped", fixed = TRUE)
  expect_length(grep("deprecated", res$warnings), 0)
  expect_identical(class(up$bmm$model), c("bmmodel", "circular", "non_targets", "imm", "imm_full"))
  expect_identical(up$bmm$model$dpar_since, "1.5.0")
  expect_s3_class(up$formula$family, "customfamily")
  expect_identical(unique(prior_rows(up)$nlpar), "")
  # c has only an intercept, so its class-b row becomes an Intercept row
  expect_true(has_prior_row(up, "normal(0, 0.7)", "Intercept", dpar = "c"))
  expect_false(has_prior_row(up, "normal(0, 0.7)", "b", dpar = "c"))
  expect_true(has_prior_row(up, "exponential(2)", "sd", dpar = "a"))
  expect_identical(set_size1_constants(up), 0L)
  expect_length(lprior_lines(up), 5)
  expect_identical(lprior_lines(up), lprior_lines(new_mock_fit(
    bmf(kappa ~ 1, c ~ 1, a ~ 0 + set_size + (0 + set_size || ID), s ~ 1),
    old_nt_model("imm"),
    brms::set_prior("normal(0, 0.7)", class = "Intercept", dpar = "c") +
      brms::set_prior("exponential(2)", class = "sd", dpar = "a")
  )))
})

test_that("update() of a 1.4 fit warns about an nlpar prior the user passes", {
  skip_on_cran()
  res <- refit_old("bmmfit_v141_mixture2p.rds",
                   prior = brms::set_prior("normal(0, 0.5)", class = "b", nlpar = "thetat"))
  expect_length(grep("deprecated", res$warnings), 1)
  expect_match(grep("deprecated", res$warnings, value = TRUE), 'dpar = "thetat"', fixed = TRUE)
  expect_length(grep("custom family", res$messages), 1)
  expect_true(has_prior_row(res$fit, "normal(0, 0.5)", "Intercept", dpar = "thetat"))
  expect_true("  lprior += normal_lpdf(Intercept_thetat | 0, 0.5);" %in% lprior_lines(res$fit))
})

# 1.4 syntax keeps working (D2) ------------------------------------------------

# the three models with the formulas of their other parameters, on two
# participants; mixture3p and imm need all seven non-target columns
migration_case <- function(name, ...) {
  nt <- paste0("col_nt", 1:7)
  switch(name,
    mixture2p = list(model = mixture2p("dev_rad", ...), rest = list(thetat ~ 1)),
    mixture3p = list(
      model = mixture3p("dev_rad", nt_features = nt, set_size = "set_size", ...),
      rest = list(thetat ~ 1, thetant ~ 1)
    ),
    imm = list(
      model = imm("dev_rad", nt_features = nt, nt_distances = paste0("dist_nt", 1:7),
                  set_size = "set_size", ...),
      rest = list(c ~ 1, a ~ 1, s ~ 1)
    )
  )
}

migration_formula <- function(case, ...) {
  do.call(bmf, c(list(...), case$rest))
}

code_lprior_lines <- function(code) {
  trimws(grep("lprior \\+=", strsplit(code, "\n")[[1]], value = TRUE))
}

deprecations <- function(warnings) {
  grep("is deprecated since bmm 1.5.0", warnings, value = TRUE)
}

# the old prior gives exactly one deprecation and the Stan prior of the new one
expect_old_prior_translated <- function(name, formula_args, old, new) {
  case <- migration_case(name)
  formula <- do.call(migration_formula, c(list(case), formula_args))
  old_out <- collect_warnings(stancode(formula, small_data(), case$model, prior = old))
  new_out <- collect_warnings(stancode(formula, small_data(), case$model, prior = new))
  expect_length(deprecations(old_out$warnings), 1)
  expect_length(deprecations(new_out$warnings), 0)
  expect_identical(code_lprior_lines(old_out$value), code_lprior_lines(new_out$value), label = name)
  invisible(deprecations(old_out$warnings))
}

for (name in c("mixture2p", "mixture3p", "imm")) {
  test_that(glue("{name}: an nlpar intercept prior becomes a dpar Intercept prior"), {
    skip_on_cran()
    expect_old_prior_translated(
      name, list(kappa ~ 1),
      brms::set_prior("normal(1, 0.5)", class = "b", coef = "Intercept", nlpar = "kappa"),
      brms::set_prior("normal(1, 0.5)", class = "Intercept", dpar = "kappa")
    )
  })

  test_that(glue("{name}: an nlpar class prior covers the intercept and the slopes, as in 1.4"), {
    skip_on_cran()
    expect_old_prior_translated(
      name, list(kappa ~ 1 + set_size),
      brms::set_prior("normal(0, 0.3)", class = "b", nlpar = "kappa"),
      brms::set_prior("normal(0, 0.3)", class = "Intercept", dpar = "kappa") +
        brms::set_prior("normal(0, 0.3)", class = "b", dpar = "kappa")
    )
    # an intercept-only formula has no slopes to give the prior to
    expect_old_prior_translated(
      name, list(kappa ~ 1),
      brms::set_prior("normal(0, 0.3)", class = "b", nlpar = "kappa"),
      brms::set_prior("normal(0, 0.3)", class = "Intercept", dpar = "kappa")
    )
    # and a cell-means formula no intercept
    expect_old_prior_translated(
      name, list(kappa ~ 0 + set_size),
      brms::set_prior("normal(0, 0.3)", class = "b", nlpar = "kappa"),
      brms::set_prior("normal(0, 0.3)", class = "b", dpar = "kappa")
    )
  })

  test_that(glue("{name}: an explicit intercept row wins over the class row, as in 1.4"), {
    skip_on_cran()
    expect_old_prior_translated(
      name, list(kappa ~ 1 + set_size),
      brms::set_prior("normal(0, 0.3)", class = "b", nlpar = "kappa") +
        brms::set_prior("normal(1, 0.5)", class = "b", coef = "Intercept", nlpar = "kappa"),
      brms::set_prior("normal(0, 0.3)", class = "b", dpar = "kappa") +
        brms::set_prior("normal(1, 0.5)", class = "Intercept", dpar = "kappa")
    )
  })

  test_that(glue("{name}: slope and sd priors keep their coef and group"), {
    skip_on_cran()
    expect_old_prior_translated(
      name, list(kappa ~ 1 + set_size + (1 | ID)),
      brms::set_prior("normal(0, 0.2)", class = "b", coef = "set_size2", nlpar = "kappa") +
        brms::set_prior("exponential(3)", class = "sd", nlpar = "kappa") +
        brms::set_prior("exponential(5)", class = "sd", group = "ID", nlpar = "kappa"),
      brms::set_prior("normal(0, 0.2)", class = "b", coef = "set_size2", dpar = "kappa") +
        brms::set_prior("exponential(3)", class = "sd", dpar = "kappa") +
        brms::set_prior("exponential(5)", class = "sd", group = "ID", dpar = "kappa")
    )
  })

  test_that(glue("{name}: priors on mu1 go to mu"), {
    skip_on_cran()
    expect_old_prior_translated(
      name, list(kappa ~ 1, mu ~ 1 + (1 | ID)),
      brms::set_prior("normal(0, 0.1)", class = "Intercept", dpar = "mu1") +
        brms::set_prior("exponential(6)", class = "sd", nlpar = "mu1"),
      brms::set_prior("normal(0, 0.1)", class = "Intercept", dpar = "mu") +
        brms::set_prior("exponential(6)", class = "sd", dpar = "mu")
    )
  })

  test_that(glue("{name}: mu1 in the formula gives one deprecation and the code of mu"), {
    skip_on_cran()
    case <- migration_case(name)
    old <- collect_warnings(stancode(migration_formula(case, kappa ~ 1, mu1 ~ 1 + set_size),
                                     small_data(), case$model))
    new <- collect_warnings(stancode(migration_formula(case, kappa ~ 1, mu ~ 1 + set_size),
                                     small_data(), case$model))
    expect_identical(
      deprecations(old$warnings),
      paste(
        "`bmf(mu1 ~ 1 + set_size)` is deprecated since bmm 1.5.0 and will be removed in bmm",
        "1.6.0 at the earliest. Use `bmf(mu ~ 1 + set_size)` instead."
      )
    )
    expect_length(deprecations(new$warnings), 0)
    expect_identical(old$value, new$value)
  })

  test_that(glue("{name}: links for mu1 go to mu, with one deprecation and no spelling hint"), {
    out <- collect_warnings(migration_case(name, links = list(mu1 = "tan_half"))$model)
    expect_identical(
      out$warnings,
      paste(
        '`links = list(mu1 = "tan_half")` is deprecated since bmm 1.5.0 and will be removed in',
        'bmm 1.6.0 at the earliest. Use `links = list(mu = "tan_half")` instead.'
      )
    )
    expect_identical(out$value$links, migration_case(name)$model$links)
  })
}

test_that("the deprecation names every translated row and says when the intercept moved", {
  skip_on_cran()
  warnings <- expect_old_prior_translated(
    "mixture2p", list(kappa ~ 1 + set_size),
    brms::set_prior("normal(0, 0.3)", class = "b", nlpar = "kappa") +
      brms::set_prior("normal(0, 1)", class = "b", coef = "Intercept", nlpar = "thetat"),
    brms::set_prior("normal(0, 0.3)", class = "Intercept", dpar = "kappa") +
      brms::set_prior("normal(0, 0.3)", class = "b", dpar = "kappa") +
      brms::set_prior("normal(0, 1)", class = "Intercept", dpar = "thetat")
  )
  expect_identical(
    warnings,
    paste0(
      '`set_prior("normal(0, 0.3)", class = "b", nlpar = "kappa")` is deprecated since bmm 1.5.0 ',
      "and will be removed in bmm 1.6.0 at the earliest. Use ",
      '`set_prior("normal(0, 0.3)", class = "Intercept", dpar = "kappa") + ',
      'set_prior("normal(0, 0.3)", class = "b", dpar = "kappa")` instead.\n',
      '`set_prior("normal(0, 1)", class = "b", coef = "Intercept", nlpar = "thetat")` is deprecated ',
      "since bmm 1.5.0 and will be removed in bmm 1.6.0 at the earliest. Use ",
      '`set_prior("normal(0, 1)", class = "Intercept", dpar = "thetat")` instead.\n',
      "In bmm 1.4 a class \"b\" prior also covered the intercept, so it is now given twice, ",
      "for 'kappa'.\n",
      "brms centres the predictors of distributional parameters, so the intercept prior of ",
      "'kappa' now applies at the mean of its predictors (for a factor, the average over its ",
      "levels weighted by their number of trials), not at their reference level. The slopes ",
      "keep their meaning. If the intercept prior was meant for the reference level, give it ",
      "for the average over the predictors instead."
    )
  )
})

test_that("an intercept prior keeps its meaning without predictors and says nothing about centring", {
  skip_on_cran()
  warnings <- expect_old_prior_translated(
    "mixture2p", list(kappa ~ 1 + (1 | ID)),
    brms::set_prior("normal(1, 0.5)", class = "b", coef = "Intercept", nlpar = "kappa"),
    brms::set_prior("normal(1, 0.5)", class = "Intercept", dpar = "kappa")
  )
  expect_no_match(warnings, "centres")
})

test_that("a non-linear parameter of the formula keeps its nlpar prior", {
  skip_on_cran()
  case <- migration_case("mixture2p")
  formula <- bmf(kappa ~ exp(lk), lk ~ 1 + set_size, thetat ~ 1)
  prior <- brms::set_prior("normal(1, 0.5)", class = "b", nlpar = "lk")
  out <- collect_warnings(stancode(formula, small_data(), case$model, prior = prior))
  expect_length(deprecations(out$warnings), 0)
  expect_true("lprior += normal_lpdf(b_lk | 1, 0.5);" %in% code_lprior_lines(out$value))
})

test_that("an nlpar prior on no parameter of the model is an error that names the parameters", {
  case <- migration_case("mixture2p")
  expect_error(
    configure_prior(
      case$model, small_data(),
      configure_model(case$model, check_data(case$model, small_data(), bmf(kappa ~ 1, thetat ~ 1)),
                      check_formula(case$model, small_data(), bmf(kappa ~ 1, thetat ~ 1)))$formula,
      user_prior = brms::set_prior("normal(0, 1)", class = "b", nlpar = "kappa1")
    ),
    paste0(
      '`set_prior\\("normal\\(0, 1\\)", class = "b", nlpar = "kappa1"\\)` names no parameter of ',
      "mixture2p\\(\\)"
    )
  )
})

test_that("a model without dpar_since is not translated", {
  skip_on_cran()
  case <- migration_case("mixture2p")
  case$model$dpar_since <- NULL
  expect_error(
    stancode(migration_formula(case, kappa ~ 1), small_data(), case$model,
             prior = brms::set_prior("normal(1, 0.5)", class = "b", coef = "Intercept", nlpar = "kappa")),
    "do not correspond to any model parameter"
  )
})

test_that("update() translates an old-syntax prior given in the call, with one deprecation", {
  skip_on_cran()
  stub <- stanfit_stub()
  case <- migration_case("mixture2p")
  formula <- migration_formula(case, kappa ~ 1 + set_size)
  fit <- suppressMessages(bmm(formula, small_data(), case$model, backend = "mock", mock_fit = stub,
                              rename = FALSE))
  out <- collect_warnings(suppressMessages(update(
    fit, prior = brms::set_prior("normal(1, 0.5)", class = "b", coef = "Intercept", nlpar = "kappa"),
    backend = "mock", mock_fit = stub, rename = FALSE, recompile = TRUE
  )))
  expect_length(deprecations(out$warnings), 1)
  expect_true("lprior += normal_lpdf(Intercept_kappa | 1, 0.5);" %in% code_lprior_lines(brms::stancode(out$value)))
})

test_that("emmeans() routing takes nlpar = and mu1 as the dpar, warning only when asked", {
  case <- migration_case("mixture2p")
  fit <- suppressMessages(bmm(migration_formula(case, kappa ~ 1 + set_size, mu ~ 1),
                              small_data(), case$model, backend = "mock", mock_fit = 1, rename = FALSE))

  out <- collect_warnings(.bmmfit_resolve_par(fit, dpar = NULL, nlpar = "kappa", warn = TRUE))
  expect_identical(out$value, list(dpar = "kappa", nlpar = NULL))
  expect_identical(
    out$warnings,
    paste(
      '`emmeans(..., nlpar = "kappa")` is deprecated since bmm 1.5.0 and will be removed in bmm',
      '1.6.0 at the earliest. Use `emmeans(..., dpar = "kappa")` instead.'
    )
  )
  out <- collect_warnings(.bmmfit_resolve_par(fit, dpar = "mu1", nlpar = NULL, warn = TRUE))
  expect_identical(out$value, list(dpar = "mu", nlpar = NULL))
  expect_match(out$warnings, 'Use `emmeans(..., dpar = "mu")` instead.', fixed = TRUE)

  # emm_basis() routes the same way without a second warning
  expect_silent(expect_identical(.bmmfit_resolve_par(fit, dpar = NULL, nlpar = "kappa"),
                                 list(dpar = "kappa", nlpar = NULL)))
  expect_silent(expect_identical(.bmmfit_resolve_par(fit, dpar = "kappa", nlpar = NULL),
                                 list(dpar = "kappa", nlpar = NULL)))
})

test_that("emmeans() routing leaves a real nlpar and a model without dpar_since alone", {
  case <- migration_case("mixture2p")
  fit <- suppressMessages(bmm(bmf(kappa ~ exp(lk), lk ~ 1 + set_size, thetat ~ 1),
                              small_data(), case$model, backend = "mock", mock_fit = 1, rename = FALSE))
  expect_silent(expect_identical(.bmmfit_resolve_par(fit, dpar = NULL, nlpar = "lk", warn = TRUE),
                                 list(dpar = NULL, nlpar = "lk")))
  fit$bmm$model$dpar_since <- NULL
  expect_silent(expect_identical(.bmmfit_resolve_par(fit, dpar = NULL, nlpar = "thetat", warn = TRUE),
                                 list(dpar = NULL, nlpar = "thetat")))
})

test_that("a mu1 prior keeps its class: mu1 was a dpar in 1.4 already", {
  skip_on_cran()
  warnings <- expect_old_prior_translated(
    "mixture2p", list(kappa ~ 1, mu ~ 1 + set_size),
    brms::set_prior("normal(0, 0.1)", class = "b", dpar = "mu1"),
    brms::set_prior("normal(0, 0.1)", class = "b", dpar = "mu")
  )
  expect_identical(
    warnings,
    paste(
      '`set_prior("normal(0, 0.1)", class = "b", dpar = "mu1")` is deprecated since bmm 1.5.0 and',
      'will be removed in bmm 1.6.0 at the earliest. Use `set_prior("normal(0, 0.1)", class = "b",',
      'dpar = "mu")` instead.'
    )
  )
})

test_that("the replacement code leaves out the intercept half that an explicit intercept prior replaces", {
  skip_on_cran()
  warnings <- expect_old_prior_translated(
    "mixture2p", list(kappa ~ 1 + set_size),
    brms::set_prior("normal(0, 0.3)", class = "b", nlpar = "kappa") +
      brms::set_prior("normal(1, 0.5)", class = "b", coef = "Intercept", nlpar = "kappa"),
    brms::set_prior("normal(0, 0.3)", class = "b", dpar = "kappa") +
      brms::set_prior("normal(1, 0.5)", class = "Intercept", dpar = "kappa")
  )
  expect_match(warnings, 'Use `set_prior("normal(0, 0.3)", class = "b", dpar = "kappa")` instead.',
               fixed = TRUE)
  expect_no_match(warnings, "given twice")
})

test_that("an nlpar prior on a non-linear parameter of the new formula names the parameters", {
  case <- migration_case("mixture2p")
  expect_error(
    stancode(bmf(kappa ~ exp(lk), lk ~ 1, thetat ~ 1), small_data(), case$model,
             prior = brms::set_prior("normal(0, 1)", class = "b", nlpar = "kappa")),
    "names no parameter of mixture2p"
  )
})

test_that("an old fit is told by its family, even with a model rebuilt by restructure()", {
  old <- load_fixture_fit("mock_v141_imm.rds")
  expect_true(is_mixture_family_fit(old))
  # restructure() of fits before 0.4.5 rebuilds the model with the current constructor
  old$bmm$model$dpar_since <- "1.5.0"
  expect_true(is_mixture_family_fit(old))
})

test_that("the imm version comes from its version class, wherever it sits", {
  old <- load_fixture_fit("mock_v141_imm.rds")
  # fits from bmm 0.4.5 to 0.4.7 put non_targets after the version class
  class(old$bmm$model) <- c("bmmodel", "circular", "imm", "imm_full", "non_targets")
  expect_identical(class(upgrade_mixture_family_model(old)$bmm$model),
                   c("bmmodel", "circular", "non_targets", "imm", "imm_full"))
})
