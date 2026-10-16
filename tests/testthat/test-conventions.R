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
  # the gate in the test. test-helpers-circmix.R and
  # test-model_mixture2p.R (#434, counts at the tip of the stack) and
  # test-model_rdm.R (#502) are registered ahead; the circmix sites in
  # test-model_imm.R and test-model_mixture3p.R cannot be: those files exist
  # here with none and a present file must match its count, so the circmix
  # merge registers them.
  expect_pattern_counts(
    "cmdstan_model\\(",
    c(
      "test-helpers-circmix.R" = 2L,
      "test-model_cswald.R" = 3L,
      "test-model_ezdm.R" = 3L,
      "test-model_m3.R" = 2L,
      "test-model_mixture2p.R" = 1L,
      "test-model_rdm.R" = 1L,
      "test-model_sdt_cdp.R" = 3L,
      "test-model_sdt_ranking.R" = 1L,
      "test-model_sdt_rating.R" = 2L,
      "test-model_sdt_yn.R" = 1L,
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
      "test-bmm_setup.R" = 18L,
      "test-model_rdm.R" = 4L,
      "test-update.R" = 2L
    ),
    rule = "R1: backend = \"cmdstanr\" (any spacing or quotes) per file"
  )
})

test_that("compiled Stan programs are sampled only at the registered sites (R1)", {
  # a compiled model can sample without bmm() and its backend argument, so
  # $sample() is counted on its own. test-helpers-circmix.R and
  # test-model_mixture2p.R (#434) and test-model_rdm.R (#502) are registered
  # ahead; the circmix merge registers its sites in test-model_imm.R and
  # test-model_mixture3p.R, as for the compile sites.
  expect_pattern_counts(
    "\\$sample\\(",
    c(
      "test-helpers-circmix.R" = 1L,
      "test-model_cswald.R" = 3L,
      "test-model_ezdm.R" = 1L,
      "test-model_m3.R" = 1L,
      "test-model_mixture2p.R" = 1L,
      "test-model_rdm.R" = 1L,
      "test-model_sdt_cdp.R" = 1L,
      "test-model_sdt_ranking.R" = 1L,
      "test-model_sdt_rating.R" = 2L,
      "test-model_sdt_yn.R" = 1L,
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
