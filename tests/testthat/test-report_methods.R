# Every number and name in the drafted text is recomputed here from the fit's
# own slots (fit$fit@sim, fit$version, the posterior draws, brms::nuts_params),
# not through fit_settings() or convergence_summary(), and every template word
# is written out, so text from any other source fails the comparison.

report_fixture_files <- function() {
  vapply(model_test_cases("fixture"), `[[`, character(1), "file")
}

quiet_report <- function(...) {
  suppressWarnings(report_methods(...))
}

fmt_rhat <- function(x) formatC(x, format = "f", digits = 3)
fmt_ess <- function(x) format(round(x), big.mark = ",", trim = TRUE)
marker <- function(i) paste0("[", paste(i, collapse = ", "), "]")

independent_convergence <- function(fit) {
  conv <- suppressWarnings(posterior::summarise_draws(
    brms::as_draws_array(fit), "rhat", "ess_bulk", "ess_tail"
  ))
  conv <- conv[!conv$variable %in% c("lp__", "lprior") &
    !startsWith(conv$variable, "prior_") & !is.na(conv$rhat), ]
  list(
    rhat = max(conv$rhat), bulk = min(conv$ess_bulk), tail = min(conv$ess_tail),
    divergent = sum(brms::nuts_params(fit, pars = "divergent__")$Value)
  )
}

expected_convergence_text <- function(fit) {
  conv <- independent_convergence(fit)
  divergent <- switch(as.character(min(conv$divergent, 2)),
    "0" = "There were no divergent transitions after warmup.",
    "1" = "There was 1 divergent transition after warmup.",
    paste0("There were ", conv$divergent, " divergent transitions after warmup.")
  )
  paste0(
    "Across all estimated parameters, the largest rank-normalized split R-hat was ",
    fmt_rhat(conv$rhat), ", the smallest bulk effective sample size (ESS) was ",
    fmt_ess(conv$bulk), ", and the smallest tail ESS was ", fmt_ess(conv$tail),
    ". ", divergent
  )
}

stan_version_of <- function(fit) {
  if (identical(fit$backend, "cmdstanr")) {
    return(as.character(fit$version$cmdstan))
  }
  cpp <- fit$fit@stanmodel@model_cpp$model_cppcode
  regmatches(cpp, regexec("stanc3 v([0-9.]+)", cpp))[[1]][2]
}

expected_model_text <- function(fit) {
  model <- fit$bmm$model
  sim <- fit$fit@sim
  n_model <- max(1, length(model_citation(fit)))
  n_brms <- length(citation("brms"))
  i_bmm <- n_model + 1
  i_brms <- i_bmm + seq_len(n_brms)
  i_backend <- max(i_brms) + 1
  name <- sub("\\.$", "", sub("^The ", "", model$name))
  version <- if (all(model$version %in% c(NA, "NA", ""))) "" else paste0(", version \"", model$version, "\",")
  r_version <- if (is.null(fit$version$R)) "" else paste0(" (version ", fit$version$R, ")")
  chains <- sim$chains
  draws <- posterior::ndraws(brms::as_draws_array(fit))
  thinning <- if (sim$thin == 1) "without thinning" else paste0("with a thinning interval of ", sim$thin)
  paste0(
    "We fitted the ", name, " ", marker(seq_len(n_model)), version,
    " with the R package bmm (version ", fit$version$bmm, ") ", marker(i_bmm), ". ",
    "bmm translated the model into a brms model (brms version ", fit$version$brms, ") ",
    marker(i_brms), ", which was estimated in Stan (version ", stan_version_of(fit),
    ") through the R package ", fit$backend, " (version ", fit$version[[fit$backend]], ") ",
    marker(i_backend), ", in R", r_version, " ", marker(i_backend + 1), ". ",
    "We ran ", chains, if (chains == 1) " Markov chain" else " Markov chains",
    " of the No-U-Turn sampler with ", fmt_ess(sim$iter), " iterations", if (chains > 1) " each",
    ", of which the first ", fmt_ess(sim$warmup), " were warmup, ", thinning, ", which left ",
    fmt_ess(draws), " post-warmup draws. ",
    "The bmm formula of the model follows this paragraph; a parameter written with = is fixed to that constant."
  )
}

expected_references <- function(fit) {
  refs <- as.vector(model_citation(fit))
  if (!length(refs)) refs <- "[CITATION NEEDED]"
  c(
    refs,
    format(citation("bmm"), style = "text"),
    format(citation("brms"), style = "text"),
    format(citation(fit$backend), style = "text"),
    format(citation(), style = "text")
  )
}

summary_warnings <- function(fit) {
  warnings <- character()
  withCallingHandlers(
    utils::capture.output(summary(fit)),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  warnings
}

# a blank constructor reference, as for a model whose reference is still missing
blank_citation <- function(fit) {
  constructor <- get_model(intersect(rev(class(fit$bmm$model)), model_names())[1])
  local_mocked_bindings(
    get_model = function(model) {
      function(...) {
        model <- constructor(...)
        model$citation <- character()
        model
      }
    },
    .env = parent.frame()
  )
  fit$bmm$model$citation <- ""
  fit
}

test_that("report_methods() and report_convergence() refuse what is not a fit with draws", {
  expect_error(report_methods(list()), "bmmfit")
  expect_error(report_convergence(data.frame()), "bmmfit")
  not_sampled <- structure(list(fit = NULL), class = c("bmmfit", "brmsfit"))
  expect_error(report_methods(not_sampled), "no posterior draws")
  expect_error(report_convergence(not_sampled), "no posterior draws")
})

test_that("report_methods() refuses a mock fit and unknown sections or formats", {
  skip_on_cran()
  mock <- load_fixture_fit("mock_bmmfit_mixture2p.rds")
  expect_error(report_methods(mock), "no posterior draws")
  expect_error(report_convergence(mock), "no posterior draws")
  fit <- load_fixture_fit("bmmfit_example1.rds")
  expect_error(report_methods(fit, sections = "results"), "should be one of")
  expect_error(report_methods(fit, sections = NULL), "at least one of 'model', 'priors' or 'convergence'")
  expect_error(report_methods(fit, sections = character(0)), "at least one of 'model', 'priors' or 'convergence'")
  expect_error(report_methods(fit, sections = NA), "at least one")
  expect_error(report_convergence(fit, format = "html"), "should be one of")
})

test_that("the paragraph of every fixture holds only the fit's values and bmm's template", {
  skip_on_cran()
  for (file in report_fixture_files()) {
    fit <- load_fixture_fit(file)
    report <- quiet_report(fit)
    prior_text <- prior_report_text(suppressWarnings(prior_info(fit, format = "text")))
    expect_identical(
      report$paragraph,
      paste(expected_model_text(fit), prior_text, expected_convergence_text(fit)),
      label = file
    )
    expect_identical(report$sections, c("model", "priors", "convergence"), label = file)
  }
})

test_that("the template itself is written out for the sdm fixture", {
  skip_on_cran()
  fit <- load_fixture_fit("bmmfit_example1.rds")
  report <- quiet_report(fit, sections = "model")
  n_brms <- length(citation("brms"))
  expect_identical(report$paragraph, paste0(
    "We fitted the Signal Discrimination Model (SDM) by Oberauer (2023) [1], ",
    "version \"simple\", with the R package bmm (version 1.3.2.9000) [2]. ",
    "bmm translated the model into a brms model (brms version 2.23.0) ",
    marker(2 + seq_len(n_brms)), ", which was estimated in Stan (version 2.39.0) ",
    "through the R package rstan (version 2.32.7) ", marker(3 + n_brms), ", in R ",
    marker(4 + n_brms), ". We ran 1 Markov chain of the No-U-Turn sampler with 400 ",
    "iterations, of which the first 200 were warmup, without thinning, which left 200 ",
    "post-warmup draws. The bmm formula of the model follows this paragraph; a parameter ",
    "written with = is fixed to that constant."
  ))
  expect_identical(report$formula, c("mu = 0", "c ~ 0 + set_size", "kappa ~ 1"))
  expect_identical(quiet_report(fit, sections = "priors")$paragraph, paste(
    "Priors were specified on the sampling (link) scale of each parameter.",
    "For parameters with a non-identity link, priors therefore apply to the transformed parameter.",
    "For the memory strength parameter c (log link), all population-level coefficients received a",
    "student_t(5, 2, 0.75) prior (bmm default). The location parameter mu was fixed to 0.",
    "For the precision parameter kappa (log link), the intercept received a student_t(5, 1.75, 0.75)",
    "prior (bmm default)."
  ))
})

test_that("a model name that ends in the period of its reference does not end the sentence", {
  skip_on_cran()
  fit <- load_fixture_fit("bmmfit_example1.rds")
  fit$bmm$model$name <- "Two-parameter mixture model by Zhang and Luck (2008)."
  paragraph <- quiet_report(fit, sections = "model")$paragraph
  expect_true(startsWith(
    paragraph,
    "We fitted the Two-parameter mixture model by Zhang and Luck (2008) [1], version"
  ))
  expect_no_match(paragraph, "(2008). [", fixed = TRUE)
})

test_that("the formula lists the user's formula with its constants", {
  skip_on_cran()
  for (file in report_fixture_files()) {
    fit <- load_fixture_fit(file)
    expect_identical(
      quiet_report(fit)$formula,
      summarise_formula(fit$bmm$user_formula, fit$bmm$model),
      label = file
    )
  }
  m3 <- load_fixture_fit("bmmfit_m3_ppcheck.rds")
  expect_true("b = 0.1" %in% quiet_report(m3)$formula)
})

test_that("the references are the model's and the installed citation() entries, verbatim", {
  skip_on_cran()
  for (file in report_fixture_files()) {
    fit <- load_fixture_fit(file)
    expect_identical(quiet_report(fit)$references, expected_references(fit), label = file)
  }
})

test_that("a model without a reference gets a placeholder and a note", {
  skip_on_cran()
  fit <- blank_citation(load_fixture_fit("bmmfit_example1.rds"))
  expect_identical(model_citation(fit), character(0))
  report <- quiet_report(fit)
  expect_identical(report$references[1], "[CITATION NEEDED]")
  expect_identical(report$references[-1], expected_references(fit)[-1])
  expect_match(report$paragraph, "by Oberauer (2023) [1], version", fixed = TRUE)
  expect_true(any(grepl("[CITATION NEEDED]", report$notes, fixed = TRUE)))
})

test_that("the notes remind to cite bmm and name what bmm cannot cite", {
  skip_on_cran()
  fit <- load_fixture_fit("bmmfit_m3_ppcheck.rds")
  notes <- quiet_report(fit)$notes
  expect_true(any(grepl("cite bmm [2]", notes, fixed = TRUE)))
  expect_true(any(grepl(attr(model_citation(fit), "uncited"), notes, fixed = TRUE)))
  sdm_notes <- quiet_report(load_fixture_fit("bmmfit_example1.rds"))$notes
  expect_false(any(grepl("user-defined", sdm_notes)))
})

test_that("sections subset and order the paragraph", {
  skip_on_cran()
  fit <- load_fixture_fit("bmmfit_example1.rds")
  model <- expected_model_text(fit)
  priors <- prior_report_text(prior_info(fit, format = "text"))
  convergence <- expected_convergence_text(fit)

  reordered <- quiet_report(fit, sections = c("convergence", "model"))
  expect_identical(reordered$paragraph, paste(convergence, model))
  expect_identical(reordered$sections, c("convergence", "model"))

  only_priors <- quiet_report(fit, sections = "priors")
  expect_identical(only_priors$paragraph, priors)
  expect_identical(only_priors$formula, character(0))
  expect_identical(only_priors$references, expected_references(fit))

  expect_identical(
    quiet_report(fit, sections = c("priors", "model", "priors"))$paragraph,
    paste(priors, model)
  )
})

test_that("a note appears exactly when brms's summary() warns about R-hat", {
  skip_on_cran()
  for (file in report_fixture_files()) {
    fit <- load_fixture_fit(file)
    brms_warns <- any(grepl("Rhats are > 1.05", summary_warnings(fit)))
    notes <- quiet_report(fit)$notes
    expect_identical(any(grepl("R-hat above 1.05", notes)), brms_warns, label = file)
  }
  # zr is the parameter above 1.05 in the stored ezdm 4par fit
  ezdm4 <- load_fixture_fit("bmmfit_ezdm4_ppcheck.rds")
  expect_true(any(grepl("b_zr_Intercept", quiet_report(ezdm4)$notes, fixed = TRUE)))
})

test_that("divergent transitions are counted and named where brms's summary() warns", {
  skip_on_cran()
  fit <- load_fixture_fit("bmmfit_example1.rds")
  params <- attr(fit$fit@sim$samples[[1]], "sampler_params")
  params$divergent__[length(params$divergent__) - 0:2] <- 1
  attr(fit$fit@sim$samples[[1]], "sampler_params") <- params
  expect_identical(sum(brms::nuts_params(fit, pars = "divergent__")$Value), 3)
  expect_true(any(grepl("3 divergent transitions", summary_warnings(fit))))

  report <- quiet_report(fit)
  expect_match(report$paragraph, "There were 3 divergent transitions after warmup.", fixed = TRUE)
  expect_true(any(grepl("3 divergent transitions", report$notes)))

  params$divergent__[] <- 0
  params$divergent__[length(params$divergent__)] <- 1
  attr(fit$fit@sim$samples[[1]], "sampler_params") <- params
  expect_match(quiet_report(fit)$paragraph, "There was 1 divergent transition after warmup.", fixed = TRUE)
})

test_that("divergent transitions are counted after warmup only", {
  skip_on_cran()
  fit <- load_fixture_fit("bmmfit_example1.rds")
  # the fixtures store no warmup draws, so only the call's argument can pin this
  local_mocked_bindings(
    nuts_params = function(x, pars, inc_warmup = FALSE, ...) {
      stopifnot(isFALSE(inc_warmup))
      data.frame(Value = c(0, 0, 1))
    },
    .package = "brms"
  )
  expect_match(
    quiet_report(fit, sections = "convergence")$paragraph,
    "There was 1 divergent transition after warmup.",
    fixed = TRUE
  )
})

test_that("a non-sampling algorithm gets no convergence sentence, and one note", {
  skip_on_cran()
  fit <- load_fixture_fit("bmmfit_example1.rds")
  fit$algorithm <- "meanfield"
  report <- quiet_report(fit)
  expect_no_match(report$paragraph, "R-hat|Markov chain|divergent")
  expect_match(
    report$paragraph,
    "The model was run with Stan's meanfield algorithm, which returned 200 draws.",
    fixed = TRUE
  )
  expect_identical(sum(grepl("meanfield", report$notes)), 1L)
  expect_identical(
    as.character(report_convergence(fit, format = "text")),
    "The meanfield algorithm runs no Markov chains, so there are no convergence diagnostics to report."
  )
  expect_identical(nrow(report_convergence(fit)), 0L)
  expect_identical(sum(grepl("meanfield", quiet_report(fit, sections = "model")$notes)), 0L)
  convergence_only <- quiet_report(fit, sections = "convergence")
  expect_identical(convergence_only$paragraph, "")
  expect_true(startsWith(format(convergence_only), "[1] "))
  withr::local_options(bmm.color_summary = FALSE)
  printed <- utils::capture.output(print(convergence_only))
  expect_false("Methods paragraph (draft):" %in% printed)
  expect_true("References:" %in% printed)
})

test_that("thinning and the sampler are reported as the fit stores them", {
  skip_on_cran()
  fit <- load_fixture_fit("bmmfit_example1.rds")
  fit$fit@sim$thin <- 2
  expect_match(quiet_report(fit, "model")$paragraph, "warmup, with a thinning interval of 2, which", fixed = TRUE)
  fit$fit@stan_args[[1]]$algorithm <- "HMC"
  expect_match(quiet_report(fit, "model")$paragraph, "We ran 1 Markov chain with 400 iterations", fixed = TRUE)
})

test_that("the R version is reported when the fit stores it, and a note asks for it otherwise", {
  skip_on_cran()
  fit <- load_fixture_fit("bmmfit_m3_ppcheck.rds")
  expect_null(fit$version$R)
  expect_true(any(grepl("R version", quiet_report(fit)$notes)))
  fit$version$R <- package_version("4.5.2")
  report <- quiet_report(fit)
  expect_match(report$paragraph, "in R (version 4.5.2) [", fixed = TRUE)
  expect_false(any(grepl("does not store the R version", report$notes)))
  expect_true(any(grepl(paste0("R ", getRversion(), " (fit: R 4.5.2)"), report$notes, fixed = TRUE)))
})

test_that("a note says when an installed package differs from the one the fit used", {
  skip_on_cran()
  fit <- load_fixture_fit("bmmfit_ezdm3_ppcheck.rds")
  # a newer version than any installed one: an older one makes brms restructure the fit
  fit$version$brms <- package_version("99.0.0")
  notes <- quiet_report(fit)$notes
  expect_true(any(grepl(
    paste0("brms ", utils::packageVersion("brms"), " (fit: brms 99.0.0)"), notes,
    fixed = TRUE
  )))
})

test_that("print() wraps the draft and format() returns the text to paste", {
  skip_on_cran()
  fit <- load_fixture_fit("bmmfit_m3_ppcheck.rds")
  report <- quiet_report(fit)
  withr::local_options(bmm.color_summary = FALSE)
  printed <- utils::capture.output(print(report))
  expect_true(all(nchar(printed) <= 80))
  expect_true(all(c("Formula:", "References:", "Notes (for you, not for the paper):") %in% printed))
  utils::capture.output(shown <- withVisible(print(report)))
  expect_false(shown$visible)
  expect_identical(shown$value, report)

  pasted <- format(report)
  expect_length(pasted, 1)
  expect_true(startsWith(pasted, report$paragraph))
  expect_match(pasted, paste0("[1] ", report$references[1]), fixed = TRUE)
  expect_match(pasted, "corr ~ b + a + c", fixed = TRUE)
  expect_no_match(pasted, report$notes[1], fixed = TRUE)

  long <- load_fixture_fit("bmmfit_example1.rds")
  long$bmm$user_formula <- bmf(
    c ~ 0 + set_size + some_extremely_long_predictor_name_1 + another_very_long_predictor_name_2 +
      (1 + set_size | ID),
    kappa ~ 1
  )
  printed <- utils::capture.output(print(quiet_report(long, sections = "model")))
  expect_true(all(nchar(printed) <= 80))
  expect_true(any(grepl("  c ~ 0 + set_size + some_extremely_long_predictor_name_1 +", printed, fixed = TRUE)))
  expect_true(any(grepl("      another_very_long_predictor_name_2 + (1 + set_size | ID)", printed, fixed = TRUE)))
})

test_that("report_convergence() returns the numbers of every fixture", {
  skip_on_cran()
  for (file in report_fixture_files()) {
    fit <- load_fixture_fit(file)
    conv <- suppressWarnings(report_convergence(fit))
    expected <- independent_convergence(fit)
    expect_s3_class(conv, "bmm_report_convergence")
    expect_identical(max(conv$rhat), expected$rhat, label = file)
    expect_identical(min(conv$ess_bulk), expected$bulk, label = file)
    expect_identical(min(conv$ess_tail), expected$tail, label = file)
    expect_identical(attr(conv, "divergent"), expected$divergent, label = file)
    expect_identical(
      as.character(suppressWarnings(report_convergence(fit, format = "text"))),
      expected_convergence_text(fit),
      label = file
    )
  }
})

test_that("report_convergence() prints the table with its footer", {
  skip_on_cran()
  fit <- load_fixture_fit("bmmfit_example1.rds")
  withr::local_options(bmm.color_summary = FALSE)
  conv <- report_convergence(fit)
  printed <- utils::capture.output(print(conv))
  expect_true(any(grepl("^Algorithm: sampling$", printed)))
  expect_true(any(grepl("^Divergent transitions after warmup: 0$", printed)))
  expect_true(any(grepl("b_kappa_Intercept", printed)))
  expect_identical(class(conv[1:2, ]), "data.frame")
  expect_identical(conv[, "rhat"], conv$rhat)
  expect_identical(report_convergence(fit, format = "text")[1], expected_convergence_text(fit))
  text <- utils::capture.output(print(report_convergence(fit, format = "text")))
  expect_identical(paste(text, collapse = " "), expected_convergence_text(fit))
})

test_that("several chains are counted, and a single chain gets a note on R-hat", {
  skip_on_cran()
  fit <- load_fixture_fit("bmmfit_example1.rds")
  expect_true(any(grepl("single chain", quiet_report(fit)$notes)))
  combined <- brms::combine_models(fit, fit)
  report <- quiet_report(combined)
  expect_identical(
    sub(" Priors were specified.*$", "", report$paragraph),
    expected_model_text(combined)
  )
  expect_match(
    report$paragraph,
    "We ran 2 Markov chains of the No-U-Turn sampler with 400 iterations each",
    fixed = TRUE
  )
  expect_false(any(grepl("single chain", report$notes)))
})

test_that("a backend that is not installed gets a placeholder reference and a note", {
  skip_on_cran()
  fit <- load_fixture_fit("bmmfit_ezdm3_ppcheck.rds")
  fit$backend <- "notinstalledbackend"
  report <- quiet_report(fit, sections = "model")
  i_backend <- length(model_citation(fit)) + 1 + length(citation("brms")) + 1
  expect_identical(report$references[i_backend], "[CITATION NEEDED]")
  expect_match(report$paragraph, paste0("notinstalledbackend ", marker(i_backend)), fixed = TRUE)
  expect_true(any(grepl(
    paste0(
      "notinstalledbackend could not be read because it is not installed: replace [CITATION NEEDED] ",
      marker(i_backend)
    ),
    report$notes, fixed = TRUE
  )))
})

test_that("fits without a model version, and fits from old versions of bmm, are reported", {
  skip_on_cran()
  fit <- load_fixture_fit("bmmfit_example1.rds")
  fit$bmm$model$version <- NULL
  expect_no_match(quiet_report(fit, sections = "model")$paragraph, "version \"", fixed = FALSE)
  legacy <- load_fixture_fit("bmmfit_legacy_sdm_v0.2.1.rds")
  class(legacy) <- c("bmmfit", "brmsfit")
  expect_match(
    quiet_report(legacy, sections = c("model", "convergence"))$paragraph,
    "^We fitted the Signal Discrimination Model"
  )
})
