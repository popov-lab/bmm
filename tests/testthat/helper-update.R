# Shared by the update() test files, which were one file until its time
# set the four-worker run's wall time (ADR 0004, 2026-11-03).

update_mock <- function(object, ...) {
  suppressMessages(update(object, ..., backend = "mock", mock_fit = 1, rename = FALSE))
}

# the stored-frame cases come from the model-case providers
# (helper-model-<model>.R): one mock-fittable (model, formula, data) per model,
# each chosen so that the stored model frame is the hard case, such as an m3
# category with zero options in some rows, or a non-target model whose formula
# leaves brms no reason to keep the set_size column
# built once per test process: every provider builds all its cases on each call
stored_frame_cases <- local({
  cases <- NULL
  function() {
    if (is.null(cases)) cases <<- model_test_cases("stored_frame")
    cases
  }
})
