# the cases the shared tests run for sdt_mafc; see helper-model-cases.R
model_cases_sdt_mafc <- function() {
  model <- sdt_mafc(response = "n_correct", n_trials = "n_trials", m = 4)
  re_data <- data.frame(
    n_correct = rep(c(30, 24), 10), n_trials = 40, id = factor(rep(1:10, each = 2))
  )
  list(
    stored_frame = list(
      sdt_mafc = list(model = model, formula = bmf(d ~ 1 + cond), data = toy_mafc_data()),
      sdt_mafc_m_column = list(
        model = sdt_mafc("n_correct", "n_trials", m = "n_afc"),
        formula = bmf(d ~ 1 + cond), data = toy_mafc_data()
      ),
      sdt_mafc_m_predictor = list(
        model = sdt_mafc("n_correct", "n_trials", m = "n_afc"),
        formula = bmf(d ~ 1 + n_afc), data = toy_mafc_data()
      )
    ),
    sd_default = list(
      sdt_mafc = list(model = model, formula = bmf(d ~ 1 + (1 | id)), data = re_data,
                      sd = c(d = "exponential(1)"))
    ),
    init = list(sdt_mafc = list(model = model, formula = bmf(d ~ 1 + (1 | id)), data = re_data)),
    links = list(sdt_mafc = list(model = model, settable = "d"))
  )
}
