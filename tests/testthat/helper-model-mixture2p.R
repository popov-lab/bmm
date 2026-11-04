# the cases the shared tests run for mixture2p; see helper-model-cases.R
model_cases_mixture2p <- function() {
  model <- mixture2p("dev_rad")
  re_formula <- bmf(kappa ~ 1 + (1 | ID), thetat ~ 1 + (1 | ID))
  list(
    stored_frame = list(
      mixture2p = list(model = model, formula = bmf(kappa ~ 1, thetat ~ 1), data = small_data())
    ),
    sd_default = list(
      mixture2p = list(
        model = model, formula = re_formula, data = small_data(),
        sd = c(kappa = "exponential(1)", thetat = "exponential(1)")
      )
    ),
    init = list(mixture2p = list(model = model, formula = re_formula, data = small_data())),
    links = list(mixture2p = list(model = model, settable = character(0))),
    simulate = list(
      mixture2p = list(
        model = model, pars = list(kappa = 2, thetat = 0.5), design = data.frame(ID = rep(1:2, 50)),
        formula = bmf(kappa ~ 1, thetat ~ 1), simulate_pars = list(mu1 = 0, kappa = 2, thetat = 0.5)
      )
    )
  )
}
