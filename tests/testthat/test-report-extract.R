fixture_fits <- function() {
  paths <- list.files(test_path("assets"), pattern = "^bmmfit_.*\\.rds$", full.names = TRUE)
  stats::setNames(lapply(paths, readRDS), basename(paths))
}

all_model_versions <- function() {
  models <- supported_models(print_call = FALSE)
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
    expect_gt(length(refs), 0)
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
