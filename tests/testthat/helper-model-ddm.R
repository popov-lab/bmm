# the cases the shared tests run for ddm; see helper-model-cases.R
model_cases_ddm <- function() {
  model <- ddm("rt", "response")
  re_formula <- bmf(drift ~ 1 + (1 | id), bound ~ 1 + (1 | id), ndt ~ 1 + (1 | id))
  list(
    stored_frame = list(
      ddm = list(model = model, formula = bmf(drift ~ 1, bound ~ 1, ndt ~ 1), data = toy_rt_data())
    ),
    sd_default = list(
      ddm = list(
        model = model, formula = re_formula, data = toy_rt_data(),
        sd = c(drift = "exponential(1)", bound = "exponential(2)", ndt = "exponential(2)")
      )
    ),
    init = list(ddm = list(model = model, formula = re_formula, data = toy_rt_data())),
    links = list(ddm = list(model = model, settable = c("drift", "bound", "ndt", "zr"))),
    fixture = list(ddm = list(file = "bmmfit_ddm_ppcheck.rds")),
    observables = list(
      ddm = list(
        model = model,
        data = withr::with_seed(1, rddm(20, drift = 1, bound = 1.5, ndt = 0.3)),
        formula = bmf(drift ~ 1)
      )
    )
  )
}
