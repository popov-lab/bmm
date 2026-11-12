# Conventions of the test suite (ADR 0004). These tests read the other test
# files and pin the count of a few patterns per file, each as a literal, so a
# pattern the rules forbid or gate cannot creep in unnoticed: the failure
# names the file and the rule. Files that open pull requests add are
# registered ahead with the counts on their branch and are ignored while
# absent. Since #422 the R2, R7 and R8 counts are zero everywhere; R1 lists
# the gated compile sites.

suite_files <- function() {
  setdiff(list.files(test_path("."), pattern = "^test-.*\\.R$"), "test-conventions.R")
}

pattern_counts <- function(pattern) {
  files <- suite_files()
  counts <- vapply(files, function(f) {
    sum(grepl(pattern, readLines(test_path(f), warn = FALSE), perl = TRUE))
  }, integer(1))
  counts[order(names(counts), method = "radix")]
}

# `expected` names every file allowed a non-zero count: a present file must
# match its count, a file that is not merged yet is ignored, every other file
# must have none; integer() means no file may have one
expect_pattern_counts <- function(pattern, expected, rule) {
  if (!length(expected)) expected <- stats::setNames(integer(), character())
  counts <- pattern_counts(pattern)
  present <- intersect(names(expected), names(counts))
  # waldo reports a mismatch by position, so name the files in the message
  changed <- present[counts[present] != expected[present]]
  expect_identical(counts[present], expected[present], info = paste0(
    rule, "; found vs registered: ",
    paste0(changed, " ", counts[changed], " vs ", expected[changed], collapse = ", ")
  ))
  others <- counts[setdiff(names(counts), names(expected))]
  expect_identical(names(others)[others > 0], character(0), info = rule)
}

test_that("Stan programs are compiled only at the registered sites (R1)", {
  # each test that reaches a site, directly or through a builder function,
  # opens with skip_unless_stan_tests(). A new site is a new entry here and
  # the gate in the test. test-model_rdm.R (#502) is registered ahead.
  expect_pattern_counts(
    "cmdstan_model\\(",
    c(
      "test-helpers-circmix.R" = 3L,
      "test-model_cswald.R" = 3L,
      "test-model_ezdm.R" = 3L,
      "test-model_imm.R" = 2L,
      "test-model_m3.R" = 2L,
      "test-model_mixture2p.R" = 2L,
      "test-model_mixture3p.R" = 2L,
      "test-model_mpt.R" = 1L,
      "test-model_psychometric_function.R" = 1L,
      "test-model_rdm.R" = 1L,
      "test-model_sdm-stan.R" = 1L,
      "test-model_sdt_cdp.R" = 3L,
      "test-model_sdt_mafc.R" = 1L,
      "test-model_sdt_ranking.R" = 1L,
      "test-model_sdt_rating.R" = 2L,
      "test-model_sdt_yn.R" = 2L,
      "test-model_tcc.R" = 2L,
      "test-sdt-rating-versions.R" = 1L
    ),
    rule = "R1: cmdstan_model() sites per file"
  )
})

test_that("real-backend fits in tests/testthat are the registered ones (R1)", {
  # test-bmm_setup.R passes the backend to a mocked bmm_setup() in all but two
  # tests ("a real machine check without a smoke test finds no failure" and
  # "a real smoke test compiles and samples"), which open with
  # skip_unless_stan_tests(); test-update.R passes it to carried_control(),
  # which fits nothing. Nothing else samples here: fits live in tests/local.
  # test-model_rdm.R (#502) is registered ahead.
  expect_pattern_counts(
    "backend\\s*=\\s*[\"']cmdstanr[\"']",
    c(
      "test-bmm_setup.R" = 21L,
      "test-model_rdm.R" = 4L,
      "test-update.R" = 2L
    ),
    rule = "R1: backend = \"cmdstanr\" (any spacing or quotes) per file"
  )
})

test_that("compiled Stan programs are sampled only at the registered sites (R1)", {
  # a compiled model can sample without bmm() and its backend argument, so
  # $sample() is counted on its own. test-model_rdm.R (#502) is registered
  # ahead.
  expect_pattern_counts(
    "\\$sample\\(",
    c(
      "test-helpers-circmix.R" = 1L,
      "test-model_cswald.R" = 5L,
      "test-model_ezdm.R" = 1L,
      "test-model_imm.R" = 2L,
      "test-model_m3.R" = 1L,
      "test-model_mixture2p.R" = 2L,
      "test-model_mixture3p.R" = 2L,
      "test-model_psychometric_function.R" = 1L,
      "test-model_rdm.R" = 1L,
      "test-model_sdm-stan.R" = 1L,
      "test-model_sdt_cdp.R" = 1L,
      "test-model_sdt_mafc.R" = 1L,
      "test-model_sdt_ranking.R" = 1L,
      "test-model_sdt_rating.R" = 2L,
      "test-model_sdt_yn.R" = 1L,
      "test-model_tcc.R" = 2L,
      "test-sdt-rating-versions.R" = 1L
    ),
    rule = "R1: $sample() per file"
  )
})

test_that("skip_on_ci() is not used as a tier gate (R2)", {
  expect_pattern_counts("skip_on_ci\\(\\)", integer(), rule = "R2: skip_on_ci() per file")
})

test_that("whole files skip on CRAN only at the registered cross-model contract files (R10)", {
  expect_pattern_counts(
    "^skip_on_cran\\(\\)",
    c("test-conditional_effects.R" = 1L, "test-default-priors.R" = 1L, "test-helpers-init.R" = 1L),
    rule = "R10: file-level skip_on_cran() per file"
  )
})

test_that("stancode() tests assert content, not that the code is non-empty (R7)", {
  expect_pattern_counts("expect_true\\(nchar\\(", integer(),
                        rule = "R7: expect_true(nchar(code) > 0) per file")
})

test_that("tests seed through withr::local_seed(), not set.seed() (R8)", {
  expect_pattern_counts("^\\s*set\\.seed\\(", integer(), rule = "R8: set.seed() per file")
})

test_that("fixtures are small, short, and regenerable (R5)", {
  skip_unless_fixtures()
  fixtures <- list.files(test_path("assets"), pattern = "\\.rds$")
  size_kb <- vapply(fixtures, function(f) file.size(test_path("assets", f)) / 1024, numeric(1))
  # the draws stored in the file: the arrays are what the file carries, and
  # brms::ndraws() reads n_save, which a trimmed fit can leave wrong
  draws <- vapply(fixtures, function(f) {
    fit <- readRDS(test_path("assets", f))
    tryCatch(brms::ndraws(brms::as_draws_array(fit)), error = function(e) NA_real_)
  }, numeric(1))
  # the mock fits carry no draws
  expect_identical(names(size_kb)[size_kb > 300], character(0))
  expect_identical(names(draws)[!is.na(draws) & draws > 200], character(0))
  generator <- test_path("..", "local", "make-fixtures.R")
  skip_if_not(file.exists(generator), "fixture generator not in the built package")
  generated <- vapply(fixtures, function(f) {
    any(grepl(f, readLines(generator, warn = FALSE), fixed = TRUE))
  }, logical(1))
  expect_identical(names(generated)[!generated], character(0))
})

test_that("every model has a model-case provider with the cases every model needs (R4)", {
  # a provider is model_cases_<model>() in helper-model-<model>.R; the shared
  # tests read its cases through model_test_cases() (helper-model-cases.R)
  providers <- sub("^model_cases_", "", ls(environment(model_test_cases), pattern = "^model_cases_"))
  expect_identical(setdiff(model_names(), providers), character(0), info = paste(
    "R4: add tests/testthat/helper-model-<model>.R defining model_cases_<model>()"
  ))
  for (kind in c("stored_frame", "sd_default", "init", "links")) {
    covered <- unlist(lapply(model_test_cases(kind), function(case) {
      intersect(class(case$model), model_names())
    }))
    # a model whose check_data() consumes or creates columns also needs a
    # revert_check_data() method for its stored-frame case to pass
    expect_identical(setdiff(model_names(), covered), character(0),
                     info = paste0("R4: no '", kind, "' case for these models"))
  }
})

test_that("every model has a simulate_data() method and a simulate case per version (R4)", {
  # coverage is read per version class, because a method that dispatches on
  # model$version can miss a version silently. A model whose generator cannot
  # produce data that check_data() accepts is listed here with its reason;
  # none on develop
  exempt <- character(0)
  versions <- unlist(lapply(model_names(), function(name) {
    version <- eval(formals(get(name, mode = "function"))$version)
    if (is.null(version)) name else paste0(name, "_", version)
  }))
  covered <- unlist(lapply(model_test_cases("simulate"), function(case) class(case$model)))
  has_method <- vapply(model_names(), function(name) {
    !is.null(utils::getS3method("simulate_data", name, optional = TRUE))
  }, logical(1))
  missing <- setdiff(c(versions, model_names()[!has_method]), c(covered, exempt))
  expect_identical(missing, character(0), info = paste0(
    "R4: no 'simulate' case for these models/versions: ", collapse_comma(missing),
    "; add a `simulate` case to model_cases_<model>() and a simulate_data.<model> method; ",
    "see ?simulate_data"
  ))
})

test_that("every model constructor with ... reports arguments it does not use (#552)", {
  # a constructor is an export that builds its model with a .model_*()
  # internal; read from the namespace, so a new model is caught without a list
  ns <- asNamespace("bmm")
  exports <- Filter(function(f) is.function(f) && "..." %in% names(formals(f)),
                    mget(getNamespaceExports("bmm"), envir = ns, inherits = TRUE))
  calls <- lapply(exports, function(f) all.names(body(f)))
  constructors <- names(calls)[vapply(calls, function(x) any(grepl("^\\.model_", x)), logical(1))]
  unchecked <- constructors[!vapply(calls[constructors], function(x) {
    ".check_constructor_dots" %in% x
  }, logical(1))]
  expect_gt(length(constructors), 10)
  expect_identical(sort(unchecked), character(0),
                   info = "call .check_constructor_dots() after call <- match.call()")

  calls_to <- function(expr, pattern) {
    if (!is.call(expr)) return(list())
    own <- if (is.name(expr[[1]]) && grepl(pattern, as.character(expr[[1]]))) list(expr) else list()
    inner <- lapply(seq_along(expr)[-1], function(i) {
      if (identical(expr[[i]], quote(expr = ))) list() else calls_to(expr[[i]], pattern)
    })
    c(own, unlist(inner, recursive = FALSE))
  }
  strings <- function(expr) {
    if (is.character(expr)) return(expr)
    is_c <- is.call(expr) && identical(expr[[1]], quote(c))
    if (is_c && all(vapply(as.list(expr)[-1], is.character, logical(1)))) unlist(as.list(expr)[-1]) else NULL
  }
  # a constructor that reaches its internal through another helper or by
  # string is invisible to a static read of its body
  for (name in setdiff(constructors, unchecked)) {
    check <- as.list(match.call(ns$.check_constructor_dots,
                                calls_to(body(exports[[name]]), "^\\.check_constructor_dots$")[[1]]))
    internal_call <- calls_to(body(exports[[name]]), "^\\.model_")[[1]]
    internal <- get(as.character(internal_call[[1]]), envir = ns)
    passed <- names(internal_call)[nzchar(names(internal_call))]
    expected <- setdiff(names(formals(internal)), c(passed, "call", "..."))
    expect_identical(sort(strings(check$known) %||% character()), sort(expected), label = name,
                     info = paste0("known = must list the ", as.character(internal_call[[1]]),
                                   "() formals that ", name, "() leaves to `...`"))
    for (arg in c("known", "deprecated")) {
      expect_true(is.null(check[[arg]]) || is.character(strings(check[[arg]])),
                  label = paste0(name, ": ", arg, " is a character literal"))
    }
  }
})

test_that("distribution validators share one name and run first (#564)", {
  ns <- asNamespace("bmm")
  fns <- ls(ns, all.names = TRUE)
  exports <- getNamespaceExports("bmm")

  expect_identical(
    grep("^\\.?validate_.+_parameters$|^\\.check_.+_(pars|params|arguments|parameters)$", fns, value = TRUE),
    character(0),
    info = "name a family validator .check_<family>_args()"
  )

  # a family is whatever has a validator; the version helpers of mixture2p
  # (.check_mixture2p_version_args) have no d/p/q/r members and drop out.
  # The test is opt-in by name: a family that never adopts .check_<family>_args
  # is not inspected, and the .validate_* helpers (.validate_sdt_rates) are out
  # of scope.
  families <- sub("^\\.check_(.+)_args$", "\\1", grep("^\\.check_.+_args$", fns, value = TRUE))
  validator_of <- character()
  inspected <- character()
  for (f in families) {
    found <- grep(paste0("^[dpqr]", f, "$"), exports, value = TRUE)
    if (length(found)) inspected <- c(inspected, f)
    validator_of <- c(validator_of, stats::setNames(rep(paste0(".check_", f, "_args"), length(found)), found))
  }
  expect_gt(length(validator_of), 10)
  # a validator renamed away from this scheme drops its family from the loop
  # below without any other failure; a count floor catches that without a list
  expect_gte(length(inspected), 6L,
             label = paste0("families inspected (", paste(inspected, collapse = ", "), ")"))

  call_heads <- function(expr) {
    if (!is.call(expr)) return(character())
    fn <- expr[[1]]
    callee <- if (is.name(fn)) {
      as.character(fn)
    } else if (is.call(fn) && is.name(fn[[1]]) && as.character(fn[[1]]) %in% c("::", ":::")) {
      # base::log() must not slip past the allow-list as an anonymous head
      as.character(fn[[3]])
    } else {
      call_heads(fn)
    }
    c(callee, unlist(lapply(as.list(expr)[-1], function(e) {
      if (identical(e, quote(expr = ))) character() else call_heads(e)
    })))
  }
  preamble <- c("{", "<-", "=", "if", "||", "isTRUE", "isFALSE", "match.arg",
                ".warn_positional_log", "names", "stop_missing_args")
  statements <- function(body) {
    if (is.call(body) && identical(body[[1]], as.name("{"))) as.list(body)[-1] else list(body)
  }
  is_stub <- function(body) {
    stmts <- statements(body)
    length(stmts) == 1 && identical(stmts[[1]], quote(.NotYetImplemented()))
  }
  # covr (CI test-coverage) rewrites each statement as
  # if (TRUE) { covr:::count(<id>); <statement> }, which hides every top-level
  # call; a user-written if never opens with a covr count and is kept
  # compared by parts: a literal covr:::count in test code makes R CMD check
  # report an undeclared dependency on covr
  is_covr_count <- function(expr) {
    h <- if (is.call(expr)) expr[[1]]
    is.call(h) && identical(h[[1]], as.name(":::")) &&
      identical(as.character(h[[2]]), "covr") && identical(as.character(h[[3]]), "count")
  }
  is_covr_wrapper <- function(expr) {
    is.call(expr) && identical(expr[[1]], as.name("if")) && length(expr) == 3 &&
      isTRUE(expr[[2]]) && is.call(expr[[3]]) && identical(expr[[3]][[1]], as.name("{")) &&
      length(expr[[3]]) == 3 && is_covr_count(expr[[3]][[2]])
  }
  uncovr <- function(expr) {
    if (!is.call(expr)) return(expr)
    if (is_covr_wrapper(expr)) return(uncovr(expr[[3]][[3]]))
    as.call(lapply(as.list(expr), function(e) {
      if (identical(e, quote(expr = ))) e else uncovr(e)
    }))
  }

  for (name in names(validator_of)) {
    body <- uncovr(body(get(name, envir = ns)))
    if (is_stub(body)) next
    stmts <- statements(body)
    at <- which(vapply(stmts, function(s) {
      is.call(s) && identical(s[[1]], as.name(validator_of[[name]]))
    }, logical(1)))
    expect_gte(length(at), 1L, label = paste0(name, " calls ", validator_of[[name]], "() at top level"))
    if (length(at)) {
      before <- as.character(unlist(lapply(stmts[seq_len(at[[1]] - 1)], call_heads)))
      expect_identical(setdiff(before, preamble), character(0), label = name,
                       info = "the validator runs before anything but match.arg(), stop_missing_args() and the positional-log shim")
    }
  }
})
