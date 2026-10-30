# the cases the shared tests run for sdm; see helper-model-cases.R. The
# default cell is "sdm"; every other activation x scaling cell (#558) is
# "sdm_<activation>_<scaling>" and runs the sdm_cont family
model_cases_sdm <- function() {
  re_formula <- bmf(kappa ~ 1 + (1 | ID), c ~ 1 + (1 | ID))
  cells <- list(sdm = sdm("dev_rad"))
  for (activation in c("vonmises", "laplace")) {
    for (scaling in c("density", "peak", "range")) {
      if (activation == "vonmises" && scaling == "density") next
      cells[[paste("sdm", activation, scaling, sep = "_")]] <-
        sdm("dev_rad", activation = activation, scaling = scaling)
    }
  }
  cases <- function(build) lapply(cells, build)
  list(
    stored_frame = cases(function(model) {
      list(model = model, formula = bmf(c ~ 1, kappa ~ 1), data = small_data())
    }),
    sd_default = cases(function(model) {
      list(
        model = model, formula = re_formula, data = small_data(),
        sd = c(kappa = "exponential(1)", c = "exponential(1)")
      )
    }),
    init = cases(function(model) list(model = model, formula = re_formula, data = small_data())),
    links = cases(function(model) list(model = model, settable = character(0))),
    fixture = list(sdm = list(file = "bmmfit_example1.rds"))
  )
}
