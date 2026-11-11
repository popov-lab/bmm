# the cases the shared tests run for sdm_cd; see helper-model-cases.R

# two participants of oberauer_lin_2017 as change-detection trials: half the
# trials probe the target (a probe of 0), the others an angle away from it
sdm_cd_data <- function() {
  data <- small_data()
  data$probe <- rep(c(0, 1.5, 0, -2), length.out = nrow(data))
  data$change <- rep(c(0L, 1L, 1L, 0L), length.out = nrow(data))
  data
}

model_cases_sdm_cd <- function() {
  model <- sdm_cd("change", "probe")
  re_formula <- bmf(c ~ 1 + (1 | ID), kappa ~ 1 + (1 | ID), criterion ~ 1 + (1 | ID))
  list(
    # the response consumes the probe and the target, which brms has no
    # reason to keep
    stored_frame = list(
      sdm_cd = list(
        model = model, formula = bmf(c ~ 1, kappa ~ 1, criterion ~ 1), data = sdm_cd_data()
      )
    ),
    check_data = list(sdm_cd = list(model = model, data = sdm_cd_data())),
    sd_default = list(
      sdm_cd = list(
        model = model, formula = re_formula, data = sdm_cd_data(),
        sd = c(c = "exponential(1)", kappa = "exponential(1)", criterion = "exponential(2)")
      )
    ),
    init = list(sdm_cd = list(model = model, formula = re_formula, data = sdm_cd_data())),
    links = list(sdm_cd = list(model = model, settable = character(0))),
    simulate = list(
      sdm_cd = list(
        model = model, pars = list(c = 1.5, kappa = 1, criterion = 0.3),
        design = data.frame(
          ID = rep(1:2, 50), probe = rep(c(0, 0, 2.5, -2), 25)
        ),
        formula = bmf(c ~ 1, kappa ~ 1, criterion ~ 1),
        simulate_pars = list(mu = 0, c = 1.5, kappa = 1, criterion = 0.3)
      )
    )
  )
}
