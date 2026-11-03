# Shared tests that cover every model read the cases of a model from here
# instead of naming it, so that a new model registers its cases in its own
# tests/testthat/helper-model-<model>.R and adds no line to a shared test file.
# That file defines a function model_cases_<model>() returning a list with any
# of these elements, each a named list of cases:
# - stored_frame: list(model, formula, data), for the stored-frame checks of
#   test-update-stored-frame.R and test-update-newdata.R; every model needs one
# - check_data: list(model, data), named after the model, for the check in
#   test-helpers-data.R that check_data() returns a data.frame, when the data
#   shared by the other models cannot serve the model
#   (a name that is no model's, such as ezdm_4par, is read only by the check
#   that built_data_columns() declares what check_data() adds)
# - rt_data: list(model, formula), for the RT and response checks of
#   test-helpers-data.R, for the models whose trial data are a response time
#   and a binary decision; the data are built in the test
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
# - specimen: a built model, named after the model, for the checks that read
#   the parameters and links a model declares (test-helpers-model.R,
#   test-helpers-parameters.R), when the constructor called without arguments
#   declares none, as mpt() does without trees; read through model_specimen()
# Case names start with the model's name and are unique across providers.
model_test_cases <- function(kind) {
  env <- environment(model_test_cases)
  providers <- ls(env, pattern = "^model_cases_")
  unlist(lapply(providers, function(p) env[[p]]()[[kind]]), recursive = FALSE)
}

# the model whose declared parameters and links a check reads: the specimen
# of the model's own provider (the other providers are not built), or what the
# constructor returns without arguments
model_specimen <- function(model, version = NA) {
  provider <- environment(model_test_cases)[[paste0("model_cases_", model)]]
  specimen <- provider()$specimen[[model]]
  if (!is.null(specimen)) {
    # one specimen stands for one version; for a model with versions it would
    # be checked once per version while the versions themselves go unchecked
    if (!is.na(version)) {
      stop("model_specimen(): the specimen of '", model, "' covers no version; version '", version, "' was asked for")
    }
    return(specimen)
  }
  if (is.na(version)) get_model(model)() else get_model(model)(version = version)
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
