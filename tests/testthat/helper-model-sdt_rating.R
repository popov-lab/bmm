# the cases the shared tests run for sdt_rating; see helper-model-cases.R
model_cases_sdt_rating <- function() {
  responses <- paste0("r", 1:5)
  model <- sdt_rating(responses, "stimulus")
  data <- toy_rating_data()
  trial_model <- sdt_rating("rating", "stimulus", n_ratings = 4)
  # one row per trial: 2 participants x 2 stimulus types x 6 trials
  trial_data <- data.frame(
    id = rep(1:2, each = 12), stimulus = rep(rep(0:1, each = 6), 2),
    rating = c(1, 2, 1, 3, 2, 1, 3, 4, 4, 2, 4, 3, 2, 1, 1, 1, 4, 2, 4, 3, 4, 4, 2, 4)
  )
  list(
    stored_frame = list(
      sdt_rating = list(
        model = model,
        formula = bmf(d ~ 1 + (1 | id), criterion ~ 1, spacing ~ 1, sdratio ~ 1),
        data = data
      ),
      sdt_rating_deltas = list(
        model = sdt_rating(responses, "stimulus", threshold_type = "log_distance"),
        formula = bmf(d ~ 1, criterion ~ 1, delta1 ~ 1, delta2 ~ 1, delta3 ~ 1),
        data = data
      ),
      sdt_rating_dpsdt = list(
        model = sdt_rating(responses, "stimulus", version = "dpsdt"),
        formula = bmf(d ~ 1 + (1 | id), criterion ~ 1, spacing ~ 1, Ro ~ 1),
        data = data
      ),
      # meta-d' needs an even number of categories
      sdt_rating_metad = list(
        model = sdt_rating(paste0("r", 1:4), "stimulus", version = "metad"),
        formula = bmf(d ~ 1, criterion ~ 1, spacing ~ 1, logmratio ~ 1 + (1 | id)),
        data = data[setdiff(names(data), "r5")]
      ),
      sdt_rating_trial_level = list(
        model = trial_model, formula = bmf(d ~ 1 + (1 | id), criterion ~ 1, spacing ~ 1),
        data = trial_data
      )
    ),
    sd_default = list(
      sdt_rating = list(
        model = model, data = data,
        formula = bmf(d ~ 1 + (1 | id), criterion ~ 1 + (1 | id),
                      spacing ~ 1 + (1 | id), sdratio ~ 1 + (1 | id)),
        sd = c(d = "exponential(1)", criterion = "exponential(1)",
               spacing = "exponential(2)", sdratio = "exponential(2)")
      )
    ),
    init = list(
      sdt_rating = list(model = model, data = data,
                        formula = bmf(d ~ 1 + (1 | id), criterion ~ 1, spacing ~ 1))
    ),
    # one count column per rating category
    check_data = list(
      sdt_rating = list(model = model, data = data),
      sdt_rating_trial_level = list(model = trial_model, data = trial_data)
    ),
    links = list(sdt_rating = list(model = model, settable = c("d", "criterion", "bias"))),
    simulate = list(
      sdt_rating = list(
        model = model, pars = list(d = 1, criterion = 0.3, spacing = 0.1),
        design = data.frame(stimulus = rep(0:1, 3), n_trials = 60), formula = bmf(d ~ 1, criterion ~ 1, spacing ~ 1),
        simulate_pars = list(d = 1, criterion = 0.3, spacing = 0.1, sdratio = 0)
      ),
      sdt_rating_log_distance = list(
        model = sdt_rating(responses, "stimulus", threshold_type = "log_distance"),
        pars = list(d = 1, criterion = 0.3, delta1 = 0, delta2 = 0.1, delta3 = -0.1, sdratio = 0.2),
        design = data.frame(stimulus = rep(0:1, 3), n_trials = 60),
        formula = bmf(d ~ 1, criterion ~ 1, delta1 ~ 1, delta2 ~ 1, delta3 ~ 1, sdratio ~ 1),
        simulate_pars = list(d = 1, criterion = 0.3, delta1 = 0, delta2 = 0.1, delta3 = -0.1, sdratio = 0.2)
      ),
      sdt_rating_softmax = list(
        model = sdt_rating(responses, "stimulus", threshold_type = "softmax"),
        pars = list(d = 1, bias = 0, spacing = 0, delta1 = 0.2, delta2 = -0.2),
        design = data.frame(stimulus = rep(0:1, 3), n_trials = 60),
        formula = bmf(d ~ 1, bias ~ 1, spacing ~ 1, delta1 ~ 1, delta2 ~ 1),
        simulate_pars = list(d = 1, bias = 0, spacing = 0, delta1 = 0.2, delta2 = -0.2, sdratio = 0)
      ),
      sdt_rating_dpsdt = list(
        model = sdt_rating(responses, "stimulus", version = "dpsdt"),
        pars = list(d = 1, criterion = 0.3, spacing = 0, Ro = -0.5),
        design = data.frame(stimulus = rep(0:1, 3), n_trials = 60),
        formula = bmf(d ~ 1, criterion ~ 1, spacing ~ 1, Ro ~ 1),
        simulate_pars = list(d = 1, criterion = 0.3, spacing = 0, sdratio = 0, Ro = -0.5, Rn = -100)
      ),
      sdt_rating_metad = list(
        model = sdt_rating(paste0("r", 1:4), "stimulus", version = "metad"),
        pars = list(d = 1, criterion = 0.3, spacing = 0, logmratio = -0.2),
        design = data.frame(stimulus = rep(0:1, 3), n_trials = 60),
        formula = bmf(d ~ 1, criterion ~ 1, spacing ~ 1, logmratio ~ 1),
        simulate_pars = list(d = 1, criterion = 0.3, spacing = 0, sdratio = 0, logmratio = -0.2)
      ),
      sdt_rating_trial_level = list(
        model = trial_model, pars = list(d = 1, criterion = 0.3, spacing = 0.1),
        design = data.frame(stimulus = rep(0:1, 30)), formula = bmf(d ~ 1, criterion ~ 1, spacing ~ 1),
        simulate_pars = list(d = 1, criterion = 0.3, spacing = 0.1, sdratio = 0)
      )
    )
  )
}
