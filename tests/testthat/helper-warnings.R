# Every warning a call raises rather than the first matching one: some tests
# require exactly one (a second one, such as brms dropping rows, can be the bug
# itself), others look for an expected warning among several.
collect_warnings <- function(expr) {
  warnings <- character()
  value <- withCallingHandlers(expr, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = warnings)
}
