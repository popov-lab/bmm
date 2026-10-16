# the cases the shared tests run for sdm; see helper-model-cases.R
model_cases_sdm <- function() {
  model <- sdm("dev_rad")
  re_formula <- bmf(kappa ~ 1 + (1 | ID), c ~ 1 + (1 | ID))
  list(
    stored_frame = list(
      sdm = list(model = model, formula = bmf(c ~ 1, kappa ~ 1), data = small_data())
    ),
    sd_default = list(
      sdm = list(
        model = model, formula = re_formula, data = small_data(),
        sd = c(kappa = "exponential(1)", c = "exponential(1)")
      )
    ),
    init = list(sdm = list(model = model, formula = re_formula, data = small_data())),
    links = list(sdm = list(model = model, settable = character(0))),
    fixture = list(sdm = list(file = "bmmfit_example1.rds"))
  )
}
