# the cases the shared tests run for sdt_cdp; see helper-model-cases.R
model_cases_sdt_cdp <- function() {
  model <- sdt_cdp(stimulus = "stimulus", n_new = 2, n_old = 2)
  re_formula <- bmf(
    dfam ~ 1 + (1 | id), drec ~ 1 + (1 | id), criterion ~ 1 + (1 | id),
    spacing ~ 1 + (1 | id), rcrit ~ 1 + (1 | id), sigmar ~ 1 + (1 | id)
  )
  list(
    stored_frame = list(
      # a column prefix and the Know/Guess split
      sdt_cdp = list(
        model = sdt_cdp("rk", "stimulus", n_new = 2, n_old = 2),
        formula = bmf(dfam ~ 1 + (1 | id), drec ~ 1, criterion ~ 1, spacing ~ 1,
                      rcrit ~ 1, kcrit ~ 1),
        data = cdp_data(guess = TRUE)
      ),
      sdt_cdp_deltas = list(
        model = sdt_cdp(stimulus = "stimulus", n_new = 2, n_old = 2,
                        threshold_type = "log_distance"),
        formula = bmf(dfam ~ 1, drec ~ 1, criterion ~ 1, rcrit ~ 1, sigmar ~ 1,
                      delta1 ~ 1, delta2 ~ 1),
        data = cdp_data(guess = FALSE, prefix = "")
      )
    ),
    sd_default = list(
      sdt_cdp = list(
        model = model, formula = re_formula, data = cdp_data(guess = FALSE, prefix = ""),
        sd = c(dfam = "exponential(1)", drec = "exponential(1)", criterion = "exponential(2)",
               spacing = "exponential(2)", rcrit = "exponential(2)", sigmar = "exponential(2)")
      )
    ),
    init = list(
      sdt_cdp = list(model = model, formula = re_formula, data = cdp_data(guess = FALSE, prefix = ""))
    ),
    # count columns named after the declared scale (new1, know2, ...)
    check_data = list(sdt_cdp = list(
      model = sdt_cdp(stimulus = "stimulus", n_new = 1, n_old = 2),
      data = data.frame(stimulus = rep(c(0L, 1L), 25), new1 = 30L, know2 = 10L,
                        know3 = 5L, remember2 = 10L, remember3 = 5L)
    )),
    links = list(sdt_cdp = list(model = model, settable = c("dfam", "drec", "criterion", "rcrit"))),
    simulate = list(
      sdt_cdp = list(
        model = sdt_cdp("", "stimulus", n_new = 2, n_old = 2),
        pars = list(dfam = 1, drec = 0.8, criterion = 0, spacing = 0.1, rcrit = 0.2),
        design = data.frame(stimulus = rep(0:1, 3), n_trials = 60),
        formula = bmf(dfam ~ 1, drec ~ 1, criterion ~ 1, spacing ~ 1, rcrit ~ 1),
        simulate_pars = list(dfam = 1, drec = 0.8, criterion = 0, spacing = 0.1, rcrit = 0.2, sigmar = 0, rho = 0,
                             kcrit = -100),
        response = c("new1", "new2", "know3", "know4", "remember3", "remember4")
      ),
      # a free know criterion adds the guess categories
      sdt_cdp_guess = list(
        model = sdt_cdp("rk", "stimulus", n_new = 2, n_old = 2, threshold_type = "log_distance"),
        pars = list(dfam = 1, drec = 0.8, criterion = 0, delta1 = 0, delta2 = 0.1, rcrit = 0.2, kcrit = 0.5),
        design = data.frame(stimulus = rep(0:1, 3), n_trials = 60),
        formula = bmf(dfam ~ 1, drec ~ 1, criterion ~ 1, delta1 ~ 1, delta2 ~ 1, rcrit ~ 1, kcrit ~ 1),
        simulate_pars = list(dfam = 1, drec = 0.8, criterion = 0, delta1 = 0, delta2 = 0.1, rcrit = 0.2, sigmar = 0,
                             rho = 0, kcrit = 0.5),
        response = paste0("rk", c("new1", "new2", "guess3", "guess4", "know3", "know4", "remember3", "remember4"))
      )
    )
  )
}
