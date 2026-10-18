# the cases the shared tests run for ezcdm; see helper-model-cases.R
model_cases_ezcdm <- function() {
  three <- ezcdm("mean_angle", "var_angle", "mean_rt", "var_rt", "n_trials")
  four <- ezcdm("mean_angle", "var_angle", "mean_rt", "var_rt", "n_trials", version = "4par")
  toy_data <- data.frame(
    mean_angle = rep(c(0.1, -0.1), 10), var_angle = rep(c(0.3, 0.4), 10),
    mean_rt = rep(c(0.8, 0.9), 10), var_rt = rep(c(0.05, 0.06), 10),
    n_trials = 100, id = factor(rep(1:10, each = 2))
  )
  sim_data <- withr::with_seed(1, rezcdm(5, n_trials = 40, drift = 2, bound = 1.5, ndt = 0.3))
  re_formula <- bmf(drift ~ 1 + (1 | id), bound ~ 1 + (1 | id), ndt ~ 1 + (1 | id))
  re_formula_4par <- re_formula + bmf(driftangle ~ 1 + (1 | id))
  list(
    stored_frame = list(
      ezcdm = list(model = three, formula = bmf(drift ~ 1, bound ~ 1, ndt ~ 1), data = toy_data)
    ),
    sd_default = list(
      ezcdm_4par = list(
        model = four, formula = re_formula_4par, data = toy_data,
        sd = c(drift = "exponential(2)", driftangle = "exponential(4)",
               bound = "exponential(2)", ndt = "exponential(2)")
      )
    ),
    init = list(
      ezcdm = list(model = three, formula = re_formula, data = toy_data),
      ezcdm_4par = list(model = four, formula = re_formula_4par, data = toy_data)
    ),
    links = list(
      ezcdm = list(model = three, settable = c("drift", "driftangle", "bound", "ndt")),
      ezcdm_4par = list(model = four, settable = c("drift", "driftangle", "bound", "ndt"))
    ),
    fixture = list(ezcdm = list(file = "bmmfit_ezcdm_ppcheck.rds")),
    observables = list(
      ezcdm = list(model = three, data = sim_data, formula = bmf(drift ~ 1, bound ~ 1, ndt ~ 1)),
      ezcdm_4par = list(model = four, data = sim_data, formula = bmf(drift ~ 1, driftangle ~ 1))
    )
  )
}
