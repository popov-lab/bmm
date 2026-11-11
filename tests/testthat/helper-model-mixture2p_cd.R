# the cases the shared tests run for mixture2p_cd; see helper-model-cases.R

# two participants of oberauer_lin_2017 as change-detection trials: half the
# trials probe the target (a probe of 0), the others an angle away from it
mixture2p_cd_data <- function() {
  data <- small_data()
  data$probe <- rep(c(0, 1.5, 0, -2), length.out = nrow(data))
  data$change <- rep(c(0L, 1L, 1L, 0L), length.out = nrow(data))
  data
}

model_cases_mixture2p_cd <- function() {
  model <- mixture2p_cd("change", "probe")
  capacity <- function(version) {
    mixture2p_cd("change", "probe", set_size = "set_size", version = version)
  }
  re_formula <- bmf(kappa ~ 1 + (1 | ID), thetat ~ 1 + (1 | ID), criterion ~ 1 + (1 | ID))
  capacity_formula <- bmf(kappa ~ 1 + (1 | ID), K ~ 1 + (1 | ID), criterion ~ 1 + (1 | ID))
  # one row per trial: the probe relative to the target, and a set size for
  # the capacity versions
  design <- data.frame(
    ID = rep(1:2, 50),
    probe = rep(c(0, 0, 2.5, -2), 25), set_size = rep(c(2, 4), each = 50)
  )
  list(
    # the capacity versions consume set_size, which brms has no reason to keep
    stored_frame = list(
      mixture2p_cd = list(
        model = model, formula = bmf(kappa ~ 1, thetat ~ 1, criterion ~ 1),
        data = mixture2p_cd_data()
      ),
      mixture2p_cd_slot_averaging = list(
        model = capacity("slot_averaging"), formula = bmf(kappa ~ 1, K ~ 1),
        data = mixture2p_cd_data()
      )
    ),
    check_data = list(
      mixture2p_cd = list(model = model, data = mixture2p_cd_data())
    ),
    sd_default = list(
      mixture2p_cd = list(
        model = model, formula = re_formula, data = mixture2p_cd_data(),
        sd = c(kappa = "exponential(1)", thetat = "exponential(1)", criterion = "exponential(2)")
      ),
      mixture2p_cd_slot = list(
        model = capacity("slot"), formula = capacity_formula, data = mixture2p_cd_data(),
        sd = c(kappa = "exponential(1)", K = "exponential(2)", criterion = "exponential(2)")
      )
    ),
    init = list(
      mixture2p_cd = list(model = model, formula = re_formula, data = mixture2p_cd_data()),
      mixture2p_cd_slot_averaging = list(
        model = capacity("slot_averaging"), formula = capacity_formula,
        data = mixture2p_cd_data()
      )
    ),
    links = list(
      mixture2p_cd = list(model = model, settable = c("mu", "kappa", "thetat")),
      mixture2p_cd_slot = list(
        model = capacity("slot"), settable = c("mu", "kappa", "K")
      )
    ),
    simulate = list(
      mixture2p_cd = list(
        model = model, pars = list(kappa = 2, thetat = 0.5), design = design,
        formula = bmf(kappa ~ 1, thetat ~ 1),
        simulate_pars = list(mu = 0, kappa = 2, thetat = 0.5, criterion = 0)
      ),
      mixture2p_cd_slot = list(
        model = capacity("slot"), pars = list(kappa = 2, K = 1, criterion = 0.3),
        design = design, formula = bmf(kappa ~ 1, K ~ 1, criterion ~ 1),
        simulate_pars = list(mu = 0, kappa = 2, K = 1, criterion = 0.3)
      ),
      mixture2p_cd_slot_averaging = list(
        model = capacity("slot_averaging"), pars = list(kappa = 2, K = 1),
        design = design, formula = bmf(kappa ~ 1, K ~ 1),
        simulate_pars = list(mu = 0, kappa = 2, K = 1, criterion = 0)
      )
    )
  )
}
