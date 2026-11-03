# the cases the shared tests run for tcc; see helper-model-cases.R
model_cases_tcc <- function() {
  model <- tcc("dev_rad", similarity = exp(-(0:180) / 20))
  bins <- tcc("bin", similarity = exp(-(0:90) / 10), n_options = 180, units = "bins")
  # Gumbel noise: a closed form, no table
  gumbel <- tcc("dev_rad", similarity = exp(-(0:180) / 20), noise = "gumbel")
  # correlated and motor noise: a Monte Carlo table from R, small to stay cheap
  correlated <- tcc("dev_rad", similarity = exp(-(0:180) / 20), d_max = 4,
                    correlation = exp(-(0:180) / 10), correlation_draws = 6000, motor_sd = 0.03)
  data <- small_data()
  re_formula <- bmf(d ~ 1 + (1 | ID))
  # forced-choice option sets (#566, D4): two- and three-option trials
  options <- tcc("err", similarity = exp(-(0:180) / 20), options = c("o1", "o2", "o3"))
  opt_data <- data.frame(
    err = c(0, 90, 0, -90, 180, 0) * pi / 180, o1 = 0,
    o2 = c(90, 90, 180, -90, 180, -90) * pi / 180, o3 = c(NA, -90, NA, 90, -90, NA) * pi / 180,
    ID = rep(1:2, 3)
  )
  list(
    stored_frame = list(
      # brms keeps the binned response only, so the error column is rebuilt
      tcc = list(model = model, formula = bmf(d ~ 1 + set_size), data = data[order(data$set_size), ]),
      tcc_bins = list(
        model = bins, formula = bmf(d ~ 1),
        data = data.frame(bin = c(0, 3, -7, 90, 181, 12), id = 1:6)
      ),
      # the response and the option columns are rebuilt from the offsets
      tcc_options = list(model = options, formula = bmf(d ~ 1), data = opt_data),
      tcc_correlated = list(model = correlated, formula = bmf(d ~ 1), data = data)
    ),
    sd_default = list(
      tcc = list(model = model, formula = re_formula, data = data, sd = c(d = "exponential(1)")),
      tcc_gumbel = list(model = gumbel, formula = re_formula, data = data, sd = c(d = "exponential(1)")),
      tcc_options = list(model = options, formula = re_formula, data = opt_data, sd = c(d = "exponential(1)"))
    ),
    init = list(
      tcc = list(model = model, formula = re_formula, data = data),
      tcc_gumbel = list(model = gumbel, formula = re_formula, data = data),
      tcc_options = list(model = options, formula = re_formula, data = opt_data)
    ),
    links = list(
      tcc = list(model = model, settable = "d"),
      tcc_gumbel = list(model = gumbel, settable = "d"),
      tcc_options = list(model = options, settable = "d")
    ),
    observables = list(
      tcc = list(model = model, data = data, formula = bmf(d ~ 1)),
      tcc_gumbel = list(model = gumbel, data = data, formula = bmf(d ~ 1)),
      tcc_bins = list(
        model = bins, data = data.frame(bin = c(0, 3, -7, 90, 181, 12)), formula = bmf(d ~ 1)
      ),
      tcc_options = list(model = options, data = opt_data, formula = bmf(d ~ 1)),
      tcc_correlated = list(model = correlated, data = data, formula = bmf(d ~ 1))
    ),
    simulate = list(
      tcc = list(
        model = model, pars = list(d = 0.5), design = data.frame(ID = rep(1:2, 50)),
        formula = bmf(d ~ 1), simulate_pars = list(d = 0.5)
      ),
      tcc_bins = list(
        model = bins, pars = list(d = 0.5), design = data.frame(ID = rep(1:2, 50)),
        formula = bmf(d ~ 1), simulate_pars = list(d = 0.5)
      ),
      tcc_gumbel = list(
        model = gumbel, pars = list(d = 0.5), design = data.frame(ID = rep(1:2, 50)),
        formula = bmf(d ~ 1), simulate_pars = list(d = 0.5)
      ),
      tcc_options = list(
        model = options, pars = list(d = 0.5), design = opt_data[rep(1:6, 10), c("o1", "o2", "o3", "ID")],
        formula = bmf(d ~ 1), simulate_pars = list(d = 0.5)
      ),
      tcc_correlated = list(
        model = correlated, pars = list(d = 0.5), design = data.frame(ID = rep(1:2, 50)),
        formula = bmf(d ~ 1), simulate_pars = list(d = 0.5)
      ),
      tcc_motor = list(
        model = tcc("dev_rad", similarity = exp(-(0:180) / 20), motor_sd = 0.03),
        pars = list(d = 0.5), design = data.frame(ID = rep(1:2, 50)),
        formula = bmf(d ~ 1), simulate_pars = list(d = 0.5)
      )
    )
  )
}
