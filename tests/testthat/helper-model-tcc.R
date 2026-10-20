# the cases the shared tests run for tcc; see helper-model-cases.R
model_cases_tcc <- function() {
  model <- tcc("dev_rad", similarity = exp(-(0:180) / 20))
  bins <- tcc("bin", similarity = exp(-(0:90) / 10), n_options = 180, units = "bins")
  data <- small_data()
  re_formula <- bmf(d ~ 1 + (1 | ID))
  list(
    stored_frame = list(
      # brms keeps the binned response only, so the error column is rebuilt
      tcc = list(model = model, formula = bmf(d ~ 1 + set_size), data = data[order(data$set_size), ]),
      tcc_bins = list(
        model = bins, formula = bmf(d ~ 1),
        data = data.frame(bin = c(0, 3, -7, 90, 181, 12), id = 1:6)
      )
    ),
    sd_default = list(
      tcc = list(model = model, formula = re_formula, data = data, sd = c(d = "exponential(1)"))
    ),
    init = list(tcc = list(model = model, formula = re_formula, data = data)),
    links = list(tcc = list(model = model, settable = "d")),
    observables = list(
      tcc = list(model = model, data = data, formula = bmf(d ~ 1)),
      tcc_bins = list(
        model = bins, data = data.frame(bin = c(0, 3, -7, 90, 181, 12)), formula = bmf(d ~ 1)
      )
    )
  )
}
