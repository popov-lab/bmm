# the cases the shared tests run for mixture3p_cd; see helper-model-cases.R

# two participants of oberauer_lin_2017 (set sizes 1 to 8) as change-detection
# trials: the probes cycle through the target, the first non-target where a
# trial has one (an intrusion probe) and a new colour
mixture3p_cd_data <- function() {
  data <- small_data()
  type <- rep_len(1:3, nrow(data))
  intrusion <- type == 2 & as.numeric(as.character(data$set_size)) > 1
  offset <- ifelse(intrusion, data$col_nt1, c(0, 0, 1.5)[type])
  data$probe <- offset
  data$change <- rep_len(c(0L, 1L, 1L, 0L, 1L), nrow(data))
  data
}

model_cases_mixture3p_cd <- function() {
  nt <- paste0("col_nt", 1:7)
  cd3p <- function(version = "simple") {
    mixture3p_cd("change", "probe", nt_features = nt, set_size = "set_size",
                 version = version)
  }
  model <- cd3p()
  re_formula <- bmf(kappa ~ 1 + (1 | ID), thetat ~ 1 + (1 | ID), thetant ~ 1 + (1 | ID),
                    criterion ~ 1 + (1 | ID))
  capacity_formula <- bmf(kappa ~ 1 + (1 | ID), K ~ 1 + (1 | ID), ptarget ~ 1 + (1 | ID),
                          criterion ~ 1 + (1 | ID))
  # one row per trial: displays of set size 1 to 3, and the probe relative to
  # the target
  design <- nt_display_design()
  design$probe <- rep_len(c(0, 0, 2.5, -2), nrow(design))
  simulate_case <- function(version, pars, formula) {
    list(
      model = mixture3p_cd("change", "probe", nt_features = c("nt1", "nt2"),
                           set_size = "set_size", version = version),
      pars = pars, design = design, formula = formula,
      simulate_pars = c(list(mu = 0), pars, if (is.null(pars$criterion)) list(criterion = 0))
    )
  }
  list(
    # brms has no reason to keep set_size, the probe or the target in the
    # first case's frame; the second names the set size in a formula
    stored_frame = list(
      mixture3p_cd = list(model = model, formula = bmf(kappa ~ 1, thetat ~ 1, thetant ~ 1),
                          data = mixture3p_cd_data()),
      mixture3p_cd_slot_averaging = list(
        model = cd3p("slot_averaging"), formula = bmf(kappa ~ 1, K ~ 1, ptarget ~ 1 + set_size),
        data = mixture3p_cd_data()
      )
    ),
    check_data = list(
      mixture3p_cd = list(model = model, data = mixture3p_cd_data())
    ),
    sd_default = list(
      mixture3p_cd = list(
        model = model, formula = re_formula, data = mixture3p_cd_data(),
        sd = c(kappa = "exponential(1)", thetat = "exponential(1)", thetant = "exponential(1)",
               criterion = "exponential(2)")
      ),
      mixture3p_cd_slot = list(
        model = cd3p("slot"), formula = capacity_formula, data = mixture3p_cd_data(),
        sd = c(kappa = "exponential(1)", K = "exponential(2)", ptarget = "exponential(1)",
               criterion = "exponential(2)")
      )
    ),
    init = list(
      mixture3p_cd = list(model = model, formula = re_formula, data = mixture3p_cd_data()),
      mixture3p_cd_slot_averaging = list(
        model = cd3p("slot_averaging"), formula = capacity_formula, data = mixture3p_cd_data()
      )
    ),
    links = list(
      mixture3p_cd = list(model = model, settable = c("mu", "kappa")),
      mixture3p_cd_slot = list(
        model = cd3p("slot"), settable = c("mu", "kappa", "K", "ptarget")
      )
    ),
    simulate = list(
      mixture3p_cd = simulate_case(
        "simple", list(kappa = 2, thetat = 0.5, thetant = -0.5),
        bmf(kappa ~ 1, thetat ~ 1, thetant ~ 1)
      ),
      mixture3p_cd_slot = simulate_case(
        "slot", list(kappa = 2, K = 0.5, ptarget = 1, criterion = 0.3),
        bmf(kappa ~ 1, K ~ 1, ptarget ~ 1, criterion ~ 1)
      ),
      mixture3p_cd_slot_averaging = simulate_case(
        "slot_averaging", list(kappa = 2, K = 0.5, ptarget = 1),
        bmf(kappa ~ 1, K ~ 1, ptarget ~ 1)
      )
    )
  )
}
