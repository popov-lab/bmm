# the cases the shared tests run for sdt_yn; see helper-model-cases.R
model_cases_sdt_yn <- function() {
  model <- sdt_yn(response = "n_old", stimulus = "stimulus", n_trials = "n_trials")
  data <- small_data(broeder_schuetz_2009_e3, id = "id")
  trial_model <- sdt_yn(response = "resp", stimulus = "stimulus")
  # one row per trial: 2 participants x 2 stimulus types x 5 trials
  trial_data <- data.frame(
    id = rep(1:2, each = 10), stimulus = rep(rep(0:1, each = 5), 2),
    resp = c(0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1)
  )
  list(
    stored_frame = list(
      sdt_yn = list(model = model, formula = bmf(d ~ 1, criterion ~ 1, sdratio ~ 1), data = data),
      sdt_yn_trial_level = list(
        model = trial_model, formula = bmf(d ~ 1 + (1 | id), criterion ~ 1), data = trial_data
      )
    ),
    check_data = list(sdt_yn_trial_level = list(model = trial_model, data = trial_data)),
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
    links = list(sdt_yn = list(model = model, settable = c("d", "criterion", "bias"))),
    simulate = list(
      sdt_yn = list(
        model = sdt_yn("n_old", "stimulus", "n_trials"), pars = list(d = 1.2, criterion = 0.4),
        design = data.frame(id = rep(1:2, each = 2), stimulus = rep(0:1, 2), n_trials = 50),
        formula = bmf(d ~ 1, criterion ~ 1), simulate_pars = list(d = 1.2, criterion = 0.4, sdratio = 0)
      ),
      sdt_yn_bias = list(
        model = sdt_yn("n_old", "stimulus", "n_trials"), pars = list(d = 1.2, bias = 0.2),
        design = data.frame(stimulus = 0:1, n_trials = 50),
        formula = bmf(d ~ 1, bias ~ 1), simulate_pars = list(d = 1.2, bias = 0.2, sdratio = 0)
      ),
      sdt_yn_trial_level = list(
        model = trial_model, pars = list(d = 1.2, criterion = 0.4),
        design = data.frame(id = rep(1:2, each = 6), stimulus = rep(0:1, 6)),
        formula = bmf(d ~ 1, criterion ~ 1), simulate_pars = list(d = 1.2, criterion = 0.4, sdratio = 0)
      )
    )
  )
}
