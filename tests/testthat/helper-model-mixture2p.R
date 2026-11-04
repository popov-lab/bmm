# the cases the shared tests run for mixture2p; see helper-model-cases.R
model_cases_mixture2p <- function() {
  model <- mixture2p("dev_rad")
  # slot averaging under variable precision has every parameter a version adds
  capacity_vp <- mixture2p("dev_rad", set_size = "set_size", version = "slot_averaging",
                           variable_precision = TRUE)
  re_formula <- bmf(kappa ~ 1 + (1 | ID), thetat ~ 1 + (1 | ID))
  capacity_vp_formula <- bmf(kappa ~ 1 + (1 | ID), K ~ 1 + (1 | ID), tau ~ 1 + (1 | ID))
  list(
    # the capacity versions consume set_size, which brms has no reason to keep
    stored_frame = list(
      mixture2p = list(model = model, formula = bmf(kappa ~ 1, thetat ~ 1), data = small_data()),
      mixture2p_slot = list(
        model = mixture2p("dev_rad", set_size = "set_size", version = "slot"),
        formula = bmf(kappa ~ 1, K ~ 1), data = small_data()
      )
    ),
    sd_default = list(
      mixture2p = list(
        model = model, formula = re_formula, data = small_data(),
        sd = c(kappa = "exponential(1)", thetat = "exponential(1)")
      ),
      mixture2p_slot_averaging_vp = list(
        model = capacity_vp, formula = capacity_vp_formula, data = small_data(),
        sd = c(kappa = "exponential(1)", K = "exponential(2)", tau = "exponential(2)")
      )
    ),
    init = list(
      mixture2p = list(model = model, formula = re_formula, data = small_data()),
      mixture2p_slot_averaging_vp = list(
        model = capacity_vp, formula = capacity_vp_formula, data = small_data()
      )
    ),
    links = list(
      mixture2p = list(model = model, settable = c("mu", "kappa", "thetat")),
      mixture2p_slot_averaging_vp = list(model = capacity_vp, settable = c("mu", "kappa", "tau", "K"))
    ),
    # the capacity versions read set_size from the design; slot averaging with
    # variable precision covers tau
    simulate = list(
      mixture2p = list(
        model = model, pars = list(kappa = 2, thetat = 0.5), design = data.frame(ID = rep(1:2, 50)),
        formula = bmf(kappa ~ 1, thetat ~ 1), simulate_pars = list(mu = 0, kappa = 2, thetat = 0.5)
      ),
      mixture2p_slot = list(
        model = mixture2p("dev_rad", set_size = "set_size", version = "slot"),
        pars = list(kappa = 2, K = 1), design = data.frame(ID = rep(1:2, 50), set_size = rep(c(2, 4), each = 50)),
        formula = bmf(kappa ~ 1, K ~ 1), simulate_pars = list(mu = 0, kappa = 2, K = 1)
      ),
      mixture2p_slot_averaging_vp = list(
        model = capacity_vp, pars = list(kappa = 2, K = 1, tau = 0),
        design = data.frame(ID = rep(1:2, 50), set_size = rep(c(2, 4), each = 50)),
        formula = bmf(kappa ~ 1, K ~ 1, tau ~ 1), simulate_pars = list(mu = 0, kappa = 2, tau = 0, K = 1)
      )
    )
  )
}
