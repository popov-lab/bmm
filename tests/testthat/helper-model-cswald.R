# the cases the shared tests run for cswald; see helper-model-cases.R
model_cases_cswald <- function() {
  rt_formula <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1)
  re_formula <- bmf(drift ~ 1 + (1 | id), bound ~ 1 + (1 | id), ndt ~ 1 + (1 | id))
  simple <- cswald("rt", "response", version = "simple")
  crisk <- cswald("rt", "response", version = "crisk")
  list(
    stored_frame = list(
      cswald = list(model = simple, formula = rt_formula, data = toy_rt_data())
    ),
    sd_default = list(
      cswald = list(
        model = simple, formula = re_formula, data = toy_rt_data(),
        sd = c(drift = "exponential(2)", bound = "exponential(2)", ndt = "exponential(2)")
      ),
      cswald_crisk = list(
        model = crisk, formula = re_formula + bmf(zr ~ 1 + (1 | id)), data = toy_rt_data(),
        sd = c(drift = "exponential(1)", bound = "exponential(2)", ndt = "exponential(2)",
               zr = "exponential(2)")
      )
    ),
    init = list(
      cswald = list(model = simple, formula = re_formula, data = toy_rt_data()),
      cswald_crisk = list(
        model = crisk, formula = re_formula + bmf(zr ~ 1), data = toy_rt_data()
      )
    ),
    links = list(
      cswald = list(model = simple, settable = c("drift", "bound", "ndt", "s")),
      cswald_crisk = list(model = crisk, settable = c("drift", "bound", "ndt", "zr", "s"))
    ),
    fixture = list(cswald = list(file = "bmmfit_cswald_ppcheck.rds")),
    observables = list(
      cswald = list(
        model = cswald(rt = "rt", response = "response"),
        data = withr::with_seed(1, rcswald(20, drift = 2, bound = 1, ndt = 0.3)),
        formula = bmf(drift ~ 1)
      ),
      cswald_crisk = list(
        model = cswald(rt = "rt", response = "response", version = "crisk"),
        data = withr::with_seed(1, rcswald(20, drift = 2, bound = 1, ndt = 0.3)),
        formula = bmf(drift ~ 1)
      )
    )
  )
}
