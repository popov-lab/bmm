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
    links = list(sdt_cdp = list(model = model, settable = c("dfam", "drec", "criterion", "rcrit")))
  )
}
