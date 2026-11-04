# the cases the shared tests run for psychometric_function(); see helper-model-cases.R
model_cases_psychometric_function <- function() {
  model <- psychometric_function("y", "x")
  trial_data <- data.frame(x = rep(c(-1, 0, 1, 2), 5), y = rep(c(0L, 1L), 10),
                           id = factor(rep(1:5, each = 4)))
  re_formula <- bmf(midpoint ~ 1 + (1 | id), width ~ 1 + (1 | id), guess ~ 1 + (1 | id),
                    lapse ~ 1 + (1 | id))
  counts <- psychometric_function("n_correct", "contrast", n_trials = "n_trials",
                                  sigmoid = "weibull", guess = "chance", lapse = 0)
  count_data <- data.frame(
    contrast = rep(c(0.01, 0.02, 0.04, 0.08), 4), n_trials = 20L,
    n_correct = rep(c(9L, 12L, 16L, 19L), 4), chance = rep(c(0.5, 0.25), each = 8)
  )
  list(
    stored_frame = list(
      psychometric_function = list(
        model = model, formula = bmf(midpoint ~ 1, width ~ 1, guess ~ 1, lapse ~ 1),
        data = trial_data
      ),
      psychometric_function_counts = list(model = counts, formula = bmf(midpoint ~ 1, width ~ 1),
                                 data = count_data)
    ),
    sd_default = list(
      # midpoint's rate is scaled to the intensity range, 4 / R with R = 3 here
      psychometric_function = list(
        model = model, formula = re_formula, data = trial_data,
        sd = c(midpoint = "exponential(1.33)", width = "exponential(2)",
               guess = "exponential(1)", lapse = "exponential(1)")
      )
    ),
    init = list(
      psychometric_function = list(model = model, formula = re_formula, data = trial_data)
    ),
    # the shared columns have no intensity, so stimulus doubles as a two-level one
    check_data = list(psychometric_function = list(
      model = psychometric_function("response", "stimulus", n_trials = "n_trials"),
      data = data.frame(response = rep(80, 50), stimulus = rep(c(0L, 1L), 25),
                        n_trials = rep(100, 50))
    )),
    links = list(psychometric_function = list(model = model, settable = character(0))),
    simulate = list(
      psychometric_function = list(
        model = model, pars = list(midpoint = 0.5, width = 0, guess = -2, lapse = -3),
        design = data.frame(x = rep(c(-1, 0, 1, 2), 5), id = rep(1:2, each = 10)),
        formula = bmf(midpoint ~ 1, width ~ 1, guess ~ 1, lapse ~ 1),
        simulate_pars = list(midpoint = 0.5, width = 0, guess = -2, lapse = -3)
      ),
      psychometric_function_counts = list(
        model = counts, pars = list(midpoint = log(0.03), width = 0.5),
        design = count_data[c("contrast", "n_trials", "chance")],
        formula = bmf(midpoint ~ 1, width ~ 1),
        simulate_pars = list(midpoint = log(0.03), width = 0.5)
      )
    )
  )
}
