# the cases the shared tests run for mixture3p; see helper-model-cases.R
model_cases_mixture3p <- function() {
  nt <- paste0("col_nt", 1:7)
  model <- mixture3p("dev_rad", nt_features = nt, set_size = "set_size")
  # slot averaging under variable precision has every parameter a version adds
  capacity_vp <- mixture3p("dev_rad", nt_features = nt, set_size = "set_size",
                           version = "slot_averaging", variable_precision = TRUE)
  re_formula <- bmf(kappa ~ 1 + (1 | ID), thetat ~ 1 + (1 | ID), thetant ~ 1 + (1 | ID))
  capacity_vp_formula <- bmf(kappa ~ 1 + (1 | ID), K ~ 1 + (1 | ID), ptarget ~ 1 + (1 | ID),
                             tau ~ 1 + (1 | ID))
  list(
    # brms has no reason to keep set_size in the first case's frame; the second
    # names it, so brms keeps the real column instead
    stored_frame = list(
      mixture3p = list(model = model, formula = bmf(kappa ~ 1, thetat ~ 1, thetant ~ 1),
                       data = small_data()),
      mixture3p_set_size = list(model = model,
                                formula = bmf(kappa ~ 1, thetat ~ 1, thetant ~ 0 + set_size),
                                data = small_data()),
      mixture3p_slot = list(
        model = mixture3p("dev_rad", nt_features = nt, set_size = "set_size", version = "slot"),
        formula = bmf(kappa ~ 1, K ~ 1, ptarget ~ 1), data = small_data()
      )
    ),
    sd_default = list(
      mixture3p = list(
        model = model, formula = re_formula, data = small_data(),
        sd = c(kappa = "exponential(1)", thetat = "exponential(1)", thetant = "exponential(1)")
      ),
      mixture3p_slot_averaging_vp = list(
        model = capacity_vp, formula = capacity_vp_formula, data = small_data(),
        sd = c(kappa = "exponential(1)", K = "exponential(2)", ptarget = "exponential(1)",
               tau = "exponential(2)")
      )
    ),
    init = list(
      mixture3p = list(model = model, formula = re_formula, data = small_data()),
      mixture3p_slot_averaging_vp = list(
        model = capacity_vp, formula = capacity_vp_formula, data = small_data()
      )
    ),
    links = list(
      mixture3p = list(model = model, settable = c("mu", "kappa")),
      mixture3p_slot_averaging_vp = list(
        model = capacity_vp, settable = c("mu", "kappa", "tau", "K", "ptarget")
      )
    )
  )
}
