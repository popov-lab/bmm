# the cases the shared tests run for sdt_yn; see helper-model-cases.R
model_cases_sdt_yn <- function() {
  model <- sdt_yn(response = "n_old", stimulus = "stimulus", n_trials = "n_trials")
  data <- small_data(broeder_schuetz_2009_e3, id = "id")
  list(
    stored_frame = list(
      sdt_yn = list(model = model, formula = bmf(d ~ 1, criterion ~ 1, sdratio ~ 1), data = data)
    ),
    sd_default = list(
      # criterion carries half of d's between-subject spread, so it joins d at rate 1
      sdt_yn = list(
        model = model, data = data,
        formula = bmf(d ~ 1 + (1 | id), criterion ~ 0 + condition + (1 | id),
                      sdratio ~ 1 + (1 | id)),
        sd = c(d = "exponential(1)", criterion = "exponential(1)", sdratio = "exponential(2)")
      ),
      sdt_yn_bias = list(
        model = model, data = data,
        formula = bmf(d ~ 1 + (1 | id), bias ~ 0 + condition + (1 | id)),
        sd = c(d = "exponential(1)", bias = "exponential(2)")
      )
    ),
    init = list(
      sdt_yn = list(model = model, data = data,
                    formula = bmf(d ~ 1 + (1 | id), criterion ~ 1 + (1 | id), sdratio ~ 1))
    ),
    links = list(sdt_yn = list(model = model, settable = c("d", "criterion", "bias")))
  )
}
