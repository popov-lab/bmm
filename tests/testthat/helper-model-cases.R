# Shared tests that cover every model read the cases of a model from here
# instead of naming it, so that a new model registers its cases in its own
# tests/testthat/helper-model-<model>.R and adds no line to a shared test file.
# That file defines a function model_cases_<model>() returning a list with any
# of these elements, each a named list of cases:
# - stored_frame: list(model, formula, data), for the stored-frame checks of
#   test-update.R; every model needs one
# - check_data: list(model, data), named after the model, for the check in
#   test-helpers-data.R that check_data() returns a data.frame, when the data
#   shared by the other models cannot serve the model
# - observables: list(model, data, formula), for test-pp_observables.R: the
#   model's pp_observables() must name standata slots and be elementwise
# - sd_default: list(model, formula, data, sd), for test-default-priors.R: sd
#   names the blanket sd prior default_prior() gives each parameter with
#   random effects in formula; every model needs one
# - init: list(model, formula, data), for test-helpers-init.R: the init
#   function bmm() builds names every Stan parameter with its declared shape;
#   every model needs one
# - links: list(model, settable), for test-links.R: settable is what
#   settable_links(model) returns, written out; every model needs one
# - fixture: list(file), for the tests that read a stored fit
#   (test-pp_check.R, test-report-extract.R): a file in tests/testthat/assets
#   written by tests/local/make-fixtures.R
# Case names start with the model's name and are unique across providers.
model_test_cases <- function(kind) {
  env <- environment(model_test_cases)
  providers <- ls(env, pattern = "^model_cases_")
  unlist(lapply(providers, function(p) env[[p]]()[[kind]]), recursive = FALSE)
}

# a mock fit of a stored-frame case, whose model frame is what update() and
# check_stored_data() get back
stored_frame_fit <- function(case) {
  # the toy rt_data has a 50% error rate, which cswald "simple" warns about
  suppressWarnings(suppressMessages(
    bmm(case$formula, case$data, case$model,
      backend = "mock", mock_fit = 1, rename = FALSE
    )
  ))
}
