# the cases the shared tests run for rdm; see helper-model-cases.R
model_cases_rdm <- function() {
  rt_data <- data.frame(
    rt = rep(c(0.6, 0.8, 1.1, 0.7), 5),
    response = rep(c(1, 2, 3, 1), 5),
    id = factor(rep(1:5, each = 4))
  )
  sim <- withr::with_seed(1, rrdm(20, drift = c(3, 1.4, 1.4), gap = 1, ndt = 0.3))
  sim$label <- c("correct", "lure", "lure")[sim$response]
  custom_formula <- bmf(correct ~ 1, lure ~ 1, gap ~ 1, ndt ~ 1)
  re_data <- transform(sim, id = rep(1:2, each = 10))
  re_shared <- bmf(gap ~ 1 + (1 | id), ndt ~ 1 + (1 | id))
  re_formula <- bmf(driftc ~ 1 + (1 | id), drifte ~ 1 + (1 | id)) + re_shared
  re_custom <- bmf(correct ~ 1 + (1 | id), lure ~ 1 + (1 | id)) + re_shared
  simple <- rdm("rt", "response", n_choices = 3)
  custom <- rdm("rt", "label", version = "custom")

  list(
    # three alternatives, so the error code check_data() lumps into 2 stands
    # for two different responses
    stored_frame = list(
      rdm = list(
        model = rdm("rt", "response", n_choices = 3),
        formula = bmf(driftc ~ 1 + (1 | id), drifte ~ 1, gap ~ 1, ndt ~ 1),
        data = rt_data
      )
    ),
    sd_default = list(
      rdm = list(
        model = simple, formula = re_formula, data = re_data,
        sd = c(driftc = "exponential(2)", drifte = "exponential(2)", gap = "exponential(2)",
               ndt = "exponential(4)")
      ),
      rdm_custom = list(
        model = custom, formula = re_custom, data = re_data,
        sd = c(correct = "exponential(2)", lure = "exponential(2)", gap = "exponential(2)",
               ndt = "exponential(4)")
      )
    ),
    init = list(
      rdm = list(model = simple, formula = re_formula, data = re_data),
      rdm_custom = list(model = custom, formula = re_custom, data = re_data)
    ),
    # the custom version learns its accumulators from the formula, so it has no
    # closed set of names to refuse
    links = list(
      rdm = list(model = simple, settable = c("driftc", "drifte", "gap", "ndt", "s", "sp")),
      rdm_custom = list(model = custom, settable = NULL)
    ),
    fixture = list(rdm = list(file = "bmmfit_rdm_ppcheck.rds")),
    check_data = list(
      rdm = list(
        model = rdm("rt", "response", n_choices = 2),
        data = data.frame(rt = rep(c(0.6, 0.8), 25), response = rep(1:2, 25))
      )
    ),
    observables = list(
      rdm = list(
        model = rdm(rt = "rt", response = "response", n_choices = 3),
        data = sim, formula = bmf(driftc ~ 1, drifte ~ 1, gap ~ 1, ndt ~ 1)
      ),
      rdm_custom = list(
        model = check_model(rdm(rt = "rt", response = "label", version = "custom"),
                            formula = custom_formula),
        data = sim, formula = custom_formula
      )
    )
  )
}
