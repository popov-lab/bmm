# the cases the shared tests run for tcc; see helper-model-cases.R
model_cases_tcc <- function() {
  model <- tcc("dev_rad", similarity = exp(-(0:180) / 20))
  bins <- tcc("bin", similarity = exp(-(0:90) / 10), n_options = 180, units = "bins")
  # a named similarity adds a parameter and a 2-D table
  laplace <- tcc("dev_rad", similarity = "laplace")
  vonmises <- tcc("dev_rad", similarity = "vonmises")
  # Gumbel noise: a closed form, no table
  gumbel <- tcc("dev_rad", similarity = exp(-(0:180) / 20), noise = "gumbel")
  gumbel_vm <- tcc("dev_rad", similarity = "vonmises", noise = "gumbel")
  data <- small_data()
  sorted <- data[order(data$set_size), ]
  re_formula <- bmf(d ~ 1 + (1 | ID))
  list(
    stored_frame = list(
      # brms keeps the binned response only, so the error column is rebuilt
      tcc = list(model = model, formula = bmf(d ~ 1 + set_size), data = sorted),
      tcc_bins = list(
        model = bins, formula = bmf(d ~ 1),
        data = data.frame(bin = c(0, 3, -7, 90, 181, 12), id = 1:6)
      ),
      tcc_laplace = list(model = laplace, formula = bmf(d ~ 1 + set_size, lambda ~ 1), data = sorted),
      tcc_gumbel_vonmises = list(
        model = gumbel_vm, formula = bmf(d ~ 1 + set_size, kappa ~ 1), data = sorted
      )
    ),
    sd_default = list(
      tcc = list(model = model, formula = re_formula, data = data, sd = c(d = "exponential(1)")),
      tcc_vonmises = list(
        model = vonmises, formula = bmf(d ~ 1 + (1 | ID), kappa ~ 1 + (1 | ID)), data = data,
        sd = c(d = "exponential(1)", kappa = "exponential(1)")
      ),
      tcc_gumbel = list(model = gumbel, formula = re_formula, data = data, sd = c(d = "exponential(1)"))
    ),
    init = list(
      tcc = list(model = model, formula = re_formula, data = data),
      tcc_laplace = list(model = laplace, formula = bmf(d ~ 1 + (1 | ID), lambda ~ 1), data = data),
      tcc_gumbel_vonmises = list(model = gumbel_vm, formula = bmf(d ~ 1 + (1 | ID), kappa ~ 1),
                                 data = data)
    ),
    links = list(
      tcc = list(model = model, settable = "d"),
      tcc_laplace = list(model = laplace, settable = c("d", "lambda")),
      tcc_vonmises = list(model = vonmises, settable = c("d", "kappa")),
      tcc_gumbel_vonmises = list(model = gumbel_vm, settable = c("d", "kappa"))
    ),
    observables = list(
      tcc = list(model = model, data = data, formula = bmf(d ~ 1)),
      tcc_gumbel = list(model = gumbel, data = data, formula = bmf(d ~ 1)),
      tcc_bins = list(
        model = bins, data = data.frame(bin = c(0, 3, -7, 90, 181, 12)), formula = bmf(d ~ 1)
      )
    )
  )
}
