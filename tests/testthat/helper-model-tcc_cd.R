# the cases the shared tests run for tcc_cd; see helper-model-cases.R
model_cases_tcc_cd <- function() {
  model <- tcc_cd("change", "dist", similarity = exp(-(0:180) / 20))
  bins <- tcc_cd("change", "step", similarity = exp(-(0:90) / 10), units = "bins")
  data <- data.frame(
    change = c(0, 1, 1, 0, 1, 0, 0, 1), dist = c(0, pi, 0.5, -0.5, 2, 0, 1, -pi),
    set_size = rep(c(2, 4), 4), ID = rep(1:2, each = 4)
  )
  list(
    stored_frame = list(
      # brms keeps the binned distance only, so the distance column is rebuilt
      tcc_cd = list(model = model, formula = bmf(d ~ 1 + set_size, criterion ~ 1), data = data),
      tcc_cd_bins = list(
        model = bins, formula = bmf(d ~ 1, bias ~ 1),
        data = data.frame(change = c(0, 1, 1, 0, 1, 0), step = c(0, 3, -7, 90, 181, 12))
      )
    ),
    check_data = list(tcc_cd = list(model = model, data = data)),
    sd_default = list(
      tcc_cd = list(model = model, data = data,
                    formula = bmf(d ~ 1 + (1 | ID), criterion ~ 1 + (1 | ID)),
                    sd = c(d = "exponential(1)", criterion = "exponential(1)")),
      tcc_cd_bias = list(model = model, data = data,
                         formula = bmf(d ~ 1 + (1 | ID), bias ~ 1 + (1 | ID)),
                         sd = c(d = "exponential(1)", bias = "exponential(2)"))
    ),
    init = list(tcc_cd = list(model = model, data = data,
                              formula = bmf(d ~ 1 + (1 | ID), criterion ~ 1 + (1 | ID)))),
    links = list(tcc_cd = list(model = model, settable = c("d", "criterion", "bias")))
  )
}
