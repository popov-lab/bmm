fixture_fits <- function() {
  paths <- list.files(test_path("assets"), pattern = "^bmmfit_.*\\.rds$", full.names = TRUE)
  stats::setNames(lapply(paths, readRDS), basename(paths))
}

all_model_versions <- function() {
  models <- model_names()
  out <- lapply(models, function(m) {
    lapply(model_versions(m), function(v) {
      if (is.na(v)) get_model(m)() else get_model(m)(version = v)
    })
  })
  unlist(out, recursive = FALSE)
}

test_that("every constructor and version meets the citation contract", {
  for (model in all_model_versions()) {
    label <- class(model)[length(class(model))]
    refs <- model_citation(model)
    expect_identical(class(refs), "character", label = label)
    expect_gt(length(refs), 0, label = label)
    expect_false(any(grepl("\n", refs, fixed = TRUE)), label = label)
    expect_false(any(grepl("^\\s*- ", refs)), label = label)
    expect_true(all(grepl("(\\.|https://doi\\.org/\\S+)$", refs)), label = label)
    expect_false(any(grepl("doi\\.org/10/", refs)), label = label)
  }
})

test_that("model_citation() splits references the constructors used to join", {
  expect_length(model_citation(.model_ezdm()), 2)
  expect_length(model_citation(.model_sdt_mafc()), 2)
  expect_length(model_citation(.model_sdt_rating(version = "dpsdt")), 2)
})

test_that("sdt_rating versions cite their own sources", {
  refs <- lapply(c("standard", "dpsdt", "metad"), function(v) {
    model_citation(.model_sdt_rating(version = v))
  })
  expect_length(unique(refs), 3)
})

test_that("a custom m3 marks its user-defined structure as uncited", {
  expect_false(is.null(attr(model_citation(.model_m3(version = "custom")), "uncited")))
  expect_null(attr(model_citation(.model_m3(version = "ss")), "uncited"))
  expect_identical(
    as.vector(model_citation(.model_m3(version = "custom"))),
    model_citation(.model_m3(version = "ss"))
  )
})

test_that("model_citation() of a fit reads the current constructor first", {
  skip_on_cran()
  path <- test_path("assets/bmmfit_example1.rds")
  skip_if_not(file.exists(path), "SDM fixture not available (excluded by .Rbuildignore)")
  fit <- readRDS(path)
  expect_identical(model_citation(fit), model_citation(.model_sdm()))
  expect_false(identical(model_citation(fit), as.character(fit$bmm$model$citation)))
})

test_that("model_citation() falls back to the copy stored on an old fit", {
  skip_on_cran()
  path <- test_path("assets/bmmfit_example1.rds")
  skip_if_not(file.exists(path), "SDM fixture not available (excluded by .Rbuildignore)")
  fit <- readRDS(path)
  blank_sdm <- function(...) {
    model <- .model_sdm(...)
    model$citation <- character()
    model
  }
  local_mocked_bindings(get_model = function(model) blank_sdm)
  expect_identical(model_citation(fit), as.character(fit$bmm$model$citation))

  fit$bmm$model$citation <- ""
  expect_identical(model_citation(fit), character(0))
})

test_that("model_citation() falls back when the model has no current constructor", {
  model <- structure(list(citation = "Doe, J. (2020). A model. Journal, 1, 1-2."), class = c("bmmodel", "retired_model"))
  expect_identical(model_citation(model), "Doe, J. (2020). A model. Journal, 1, 1-2.")
})

test_that("model_citation() falls back for a version the installed constructor dropped", {
  model <- .model_sdt_rating(version = "dpsdt")
  model$version <- "retired_version"
  expect_identical(model_citation(model), as.character(model$citation))
})

test_that("model_citation() cites the default version of a model stored without a version", {
  for (version in list(NULL, "NA", "", NA_character_)) {
    model <- .model_sdm()
    model$version <- version
    model$citation <- "stale"
    expect_identical(model_citation(model), model_citation(.model_sdm()), label = deparse(version))
  }
})

test_that("model_citation() splits the single string older fits stored into references", {
  skip_on_cran()
  paths <- c(
    cswald = test_path("assets/bmmfit_cswald_ppcheck.rds"),
    ezdm = test_path("assets/bmmfit_ezdm3_ppcheck.rds")
  )
  skip_if_not(all(file.exists(paths)), "Fixtures not available (excluded by .Rbuildignore)")
  local_mocked_bindings(current_constructor = function(model) NULL)
  refs <- lapply(paths, function(path) model_citation(readRDS(path)))
  expect_identical(lengths(refs), c(cswald = 1L, ezdm = 2L))
  refs <- unlist(refs)
  expect_false(any(grepl("\n", refs, fixed = TRUE)))
  expect_false(any(grepl("^\\s*- |  ", refs)))
  expect_true(all(startsWith(refs, c("Miller, R.", "Wagenmakers, E.-J.", "Chávez De la Peña"))))
})

test_that("model_citation() refuses objects that are not models or fits", {
  expect_error(model_citation("sdm"), "bmmodel")
  expect_error(model_citation(list(citation = "x")), "bmmodel")
})

test_that("model_docs() writes one bullet per reference", {
  docs <- model_docs(.model_ezdm(), components = "citation")
  expect_identical(lengths(regmatches(docs, gregexpr("\n   - ", docs))), 2L)
  expect_identical(
    model_docs(.model_sdm(), components = "citation"),
    paste0("* **Citation:** \n\n   - ", model_citation(.model_sdm()), "\n\n")
  )
})

test_that("use_model_template() scaffolds an empty citation vector", {
  for (versions in list(NULL, c("simple", "full"))) {
    env <- new.env()
    eval(
      parse(text = capture.output(use_model_template("tmpl_cite", versions, testing = TRUE))),
      envir = env
    )
    expect_identical(model_citation(env$.model_tmpl_cite()), character(0))
  }
})

test_that("summarise_formula() returns one line per parameter, constants included", {
  skip_on_cran()
  path <- test_path("assets/bmmfit_m3_ppcheck.rds")
  skip_if_not(file.exists(path), "M3 fixture not available (excluded by .Rbuildignore)")
  fit <- readRDS(path)
  lines <- summarise_formula(fit$bmm$user_formula, model = fit$bmm$model)
  expect_type(lines, "character")
  expect_false(any(grepl("\n", lines, fixed = TRUE)))
  constants <- fit$bmm$model$fixed_parameters
  constants <- constants[names(constants) %in% names(fit$bmm$model$parameters)]
  expect_true(all(paste0(names(constants), " = ", unlist(constants)) %in% lines))
  expect_true(all(names(fit$bmm$user_formula) %in% sub(" .*$", "", lines)))
})

test_that("print(summary()) shows the formula text summarise_formula() returns", {
  skip_on_cran()
  fits <- fixture_fits()
  skip_if(length(fits) == 0, "Fixtures not available (excluded by .Rbuildignore)")
  for (name in names(fits)) {
    fit <- fits[[name]]
    printed <- trimws(capture.output(print(suppressWarnings(summary(fit)), color = FALSE)), "right")
    lines <- summarise_formula(fit$bmm$user_formula, model = fit$bmm$model)
    expect_identical(sub("^Formula: ", "", printed[grep("^Formula: ", printed)]), lines[1], label = name)
    expect_true(all(paste0(strrep(" ", 9), lines[-1]) %in% printed), label = name)
  }
})

nuts_count <- function(fit, par, at_least = 1) {
  np <- brms::nuts_params(fit, pars = par)
  sum(np$Value >= at_least)
}

test_that("fit_settings() reads the sampler settings and versions stored on the fit", {
  skip_on_cran()
  fits <- fixture_fits()
  skip_if(length(fits) == 0, "Fixtures not available (excluded by .Rbuildignore)")
  for (name in names(fits)) {
    fit <- fits[[name]]
    settings <- fit_settings(fit)
    sim <- fit$fit@sim
    expect_identical(settings$backend, fit$backend, label = name)
    expect_identical(settings$algorithm, fit$algorithm, label = name)
    expect_equal(settings$chains, sim$chains, label = name)
    expect_equal(settings$iter, sim$iter, label = name)
    expect_equal(settings$warmup, sim$warmup, label = name)
    expect_equal(settings$thin, sim$thin, label = name)
    expect_equal(settings$ndraws_stored, posterior::ndraws(brms::as_draws_array(fit)), label = name)
    expect_identical(settings$date, fit$fit@date, label = name)
    expect_identical(
      settings$versions[c("bmm", "brms", "backend")],
      c(
        bmm = as.character(fit$version$bmm),
        brms = as.character(fit$version$brms),
        backend = as.character(fit$version[[fit$backend]])
      ),
      label = name
    )
    if (fit$backend == "cmdstanr") {
      expect_identical(settings$versions[["stan"]], as.character(fit$version$cmdstan), label = name)
    }
  }
})

test_that("fit_settings() counts the draws that summary() reports", {
  skip_on_cran()
  path <- test_path("assets/bmmfit_m3_ppcheck.rds")
  skip_if_not(file.exists(path), "M3 fixture not available (excluded by .Rbuildignore)")
  # n_save of this fixture says 80 draws, but it stores 4000
  expect_equal(fit_settings(readRDS(path))$ndraws_stored, 4000)
})

test_that("fit_settings() does not count saved warmup draws", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  path <- test_path("assets/bmmfit_example1.rds")
  skip_if_not(file.exists(path), "SDM fixture not available (excluded by .Rbuildignore)")
  fit <- readRDS(path)
  # the fixture saved no warmup; prepend 20 warmup iterations as rstan stores them
  with_warmup <- function(x) c(x[seq_len(20)], x)
  fit$fit@sim$samples <- lapply(fit$fit@sim$samples, function(chain) {
    params <- lapply(attr(chain, "sampler_params"), with_warmup)
    chain[] <- lapply(chain, with_warmup)
    attr(chain, "sampler_params") <- params
    chain
  })
  fit$fit@sim$warmup2 <- rep(20L, fit$fit@sim$chains)
  fit$fit@sim$n_save <- fit$fit@sim$n_save + 20
  expect_equal(posterior::ndraws(brms::as_draws_array(fit)), 50)
  expect_equal(fit_settings(fit)$ndraws_stored, 50)
})

test_that("fit_settings() reads the Stan version of an rstan fit from its compiled model", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  path <- test_path("assets/bmmfit_example1.rds")
  skip_if_not(file.exists(path), "SDM fixture not available (excluded by .Rbuildignore)")
  fit <- readRDS(path)
  # StanHeaders 2.32.7 built this fit, with stanc3 2.32.2
  expect_identical(fit_settings(fit)$versions[["stan"]], "2.32.2")
  fit$fit@stanmodel@model_cpp <- list()
  expect_identical(fit_settings(fit)$versions[["stan"]], as.character(fit$version$stanHeaders))
  # bmm(..., empty = TRUE) leaves no stanfit
  fit$fit <- NULL
  expect_identical(fit_settings(fit)$versions[["stan"]], as.character(fit$version$stanHeaders))
})

test_that("fit_settings() and convergence_summary() treat a stanfit without draws as a mock", {
  skip_on_cran()
  skip_if_not_installed("rstan")
  path <- test_path("assets/bmmfit_example1.rds")
  skip_if_not(file.exists(path), "SDM fixture not available (excluded by .Rbuildignore)")
  fit <- readRDS(path)
  fit$fit <- methods::new("stanfit")
  settings <- fit_settings(fit)
  expect_true(is.na(settings$ndraws_stored))
  expect_true(is.na(settings$date))
  expect_error(convergence_summary(fit), "draws")
})

test_that("fit_settings() works on mock fits, with NA where the fit has no draws", {
  skip_on_cran()
  path <- test_path("assets/mock_bmmfit_mixture2p.rds")
  skip_if_not(file.exists(path), "Mock fixture not available (excluded by .Rbuildignore)")
  settings <- fit_settings(readRDS(path))
  expect_true(all(is.na(unlist(settings[c("chains", "iter", "warmup", "thin", "ndraws_stored", "date")]))))
  expect_true(all(is.na(settings$versions)))

  dat <- data.frame(y = c(-1, 0, 1))
  mock <- bmm(bmf(c ~ 1, kappa ~ 1), dat, sdm("y"), backend = "mock", mock_fit = 1, rename = FALSE)
  settings <- fit_settings(mock)
  expect_identical(settings$backend, "mock")
  expect_identical(settings$versions[["bmm"]], as.character(utils::packageVersion("bmm")))
  expect_true(is.na(settings$ndraws_stored))
})

test_that("convergence_summary() equals R-hat and ESS computed per variable from the draws", {
  skip_on_cran()
  fits <- fixture_fits()
  skip_if(length(fits) == 0, "Fixtures not available (excluded by .Rbuildignore)")
  for (name in names(fits)) {
    fit <- fits[[name]]
    # tiny fixture fits: posterior warns that it caps the ESS estimates
    withr::local_options(list(warn = -1))
    conv <- convergence_summary(fit)
    draws <- brms::as_draws_array(fit)
    candidates <- grep("^prior_", setdiff(posterior::variables(draws), c("lp__", "lprior")), value = TRUE, invert = TRUE)
    by_var <- lapply(candidates, function(v) posterior::extract_variable_matrix(draws, v))
    rhat <- vapply(by_var, posterior::rhat, numeric(1))
    keep <- !is.na(rhat)
    expect_identical(conv$parameter, candidates[keep], label = name)
    expect_equal(conv$rhat, rhat[keep], label = name)
    expect_equal(conv$ess_bulk, vapply(by_var[keep], posterior::ess_bulk, numeric(1)), label = name)
    expect_equal(conv$ess_tail, vapply(by_var[keep], posterior::ess_tail, numeric(1)), label = name)
    expect_identical(attr(conv, "algorithm"), "sampling", label = name)
    expect_equal(attr(conv, "divergent"), nuts_count(fit, "divergent__"), label = name)
    expect_equal(
      attr(conv, "max_treedepth_hits"),
      nuts_count(fit, "treedepth__", brms::control_params(fit)$max_treedepth),
      label = name
    )
  }
})

test_that("convergence_summary() leaves out the draws from the prior", {
  skip_on_cran()
  path <- test_path("assets/bmmfit_m3_ppcheck.rds")
  skip_if_not(file.exists(path), "M3 fixture not available (excluded by .Rbuildignore)")
  withr::local_options(list(warn = -1))
  fit <- readRDS(path)
  expect_true(any(grepl("^prior_", posterior::variables(brms::as_draws_array(fit)))))
  expect_false(any(grepl("^prior_", convergence_summary(fit)$parameter)))
})

test_that("convergence_summary() marks the parameters that summary() shows", {
  skip_on_cran()
  fits <- fixture_fits()
  skip_if(length(fits) == 0, "Fixtures not available (excluded by .Rbuildignore)")
  withr::local_options(list(warn = -1))
  for (name in names(fits)) {
    fit <- fits[[name]]
    conv <- convergence_summary(fit)
    brms_summary <- brms:::summary.brmsfit(fit)
    tables <- c(list(brms_summary$fixed, brms_summary$spec_pars, brms_summary$cor_pars), brms_summary$random)
    rhat <- unlist(lapply(tables, function(x) x$Rhat))
    expect_equal(sort(conv$rhat[conv$in_summary]), sort(rhat[!is.na(rhat)]), ignore_attr = TRUE, label = name)
  }
  # distributional parameters without a formula are sampled under their own name;
  # a class name only counts as a whole prefix
  expect_identical(
    grepl(
      summary_variables_regex(fits[["bmmfit_example1.rds"]]),
      c("kappa", "c", "Intercept_c", "r_id__c[1,Intercept]", "bQ", "sdx", "cx_Intercept")
    ),
    c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE)
  )
})

test_that("convergence_summary() counts divergent and maximum-depth transitions", {
  skip_on_cran()
  path <- test_path("assets/bmmfit_example1.rds")
  skip_if_not(file.exists(path), "SDM fixture not available (excluded by .Rbuildignore)")
  withr::local_options(list(warn = -1))
  fit <- readRDS(path)
  params <- attr(fit$fit@sim$samples[[1]], "sampler_params")
  last <- length(params$divergent__) - 0:2
  params$divergent__[last] <- 1
  params$treedepth__[last] <- 10
  attr(fit$fit@sim$samples[[1]], "sampler_params") <- params
  conv <- convergence_summary(fit)
  expect_equal(attr(conv, "divergent"), 3)
  expect_equal(attr(conv, "max_treedepth_hits"), 3)
})

test_that("convergence_summary() has no rows for algorithms other than sampling", {
  skip_on_cran()
  path <- test_path("assets/bmmfit_example1.rds")
  skip_if_not(file.exists(path), "SDM fixture not available (excluded by .Rbuildignore)")
  fit <- readRDS(path)
  fit$algorithm <- "meanfield"
  conv <- convergence_summary(fit)
  expect_identical(nrow(conv), 0L)
  expect_identical(names(conv), c("parameter", "rhat", "ess_bulk", "ess_tail", "in_summary"))
  expect_identical(attr(conv, "algorithm"), "meanfield")
})

test_that("convergence_summary() refuses a mock fit", {
  dat <- data.frame(y = c(-1, 0, 1))
  mock <- bmm(bmf(c ~ 1, kappa ~ 1), dat, sdm("y"), backend = "mock", mock_fit = 1, rename = FALSE)
  expect_error(convergence_summary(mock), "draws")
})
