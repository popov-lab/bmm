# the cases the shared tests run for imm_cd; see helper-model-cases.R

# two participants of oberauer_lin_2017 (set sizes 1 to 8) as change-detection
# trials: the probes cycle through the target, the first non-target where a
# trial has one (an intrusion probe) and a new colour
imm_cd_data <- function() {
  data <- small_data()
  type <- rep_len(1:3, nrow(data))
  intrusion <- type == 2 & as.numeric(as.character(data$set_size)) > 1
  offset <- ifelse(intrusion, data$col_nt1, c(0, 0, 1.5)[type])
  data$probe <- offset
  data$change <- rep_len(c(0L, 1L, 1L, 0L, 1L), nrow(data))
  data
}

model_cases_imm_cd <- function() {
  cdimm <- function(version = "full") {
    imm_cd("change", "probe", nt_features = paste0("col_nt", 1:7),
           nt_distances = paste0("dist_nt", 1:7), set_size = "set_size", version = version)
  }
  model <- cdimm()
  re_formula <- bmf(kappa ~ 1 + (1 | ID), c ~ 1 + (1 | ID), a ~ 1 + (1 | ID), s ~ 1 + (1 | ID),
                    criterion ~ 1 + (1 | ID))
  abc_formula <- bmf(kappa ~ 1 + (1 | ID), c ~ 1 + (1 | ID), a ~ 1 + (1 | ID))
  # one row per trial: displays of set size 1 to 3, and the probe relative to
  # the target
  design <- nt_display_design()
  design$probe <- rep_len(c(0, 0, 2.5, -2), nrow(design))
  # b is fixed to 0 on its log link, so it reaches simulate_pars as 0
  simulate_case <- function(version, pars, formula) {
    distances <- if (version != "abc") c("dist1", "dist2")
    list(
      model = imm_cd("change", "probe", nt_features = c("nt1", "nt2"),
                     nt_distances = distances, set_size = "set_size", version = version),
      pars = pars, design = if (version == "abc") design[setdiff(names(design), c("dist1", "dist2"))] else design,
      formula = formula,
      simulate_pars = c(list(mu = 0), pars, list(b = 0),
                        if (is.null(pars$criterion)) list(criterion = 0))
    )
  }
  list(
    # brms has no reason to keep set_size, the probe or the target in the
    # first case's frame; the second names the set size in a formula
    stored_frame = list(
      imm_cd = list(model = model, formula = bmf(kappa ~ 1, c ~ 1, a ~ 1, s ~ 1),
                    data = imm_cd_data()),
      imm_cd_abc = list(model = cdimm("abc"), formula = bmf(kappa ~ 1, c ~ 1 + set_size, a ~ 1),
                        data = imm_cd_data())
    ),
    check_data = list(
      imm_cd = list(model = model, data = imm_cd_data())
    ),
    sd_default = list(
      imm_cd = list(
        model = model, formula = re_formula, data = imm_cd_data(),
        sd = c(kappa = "exponential(1)", c = "exponential(1)", a = "exponential(1)",
               s = "exponential(1)", criterion = "exponential(2)")
      ),
      imm_cd_abc = list(
        model = cdimm("abc"), formula = abc_formula, data = imm_cd_data(),
        sd = c(kappa = "exponential(1)", c = "exponential(1)", a = "exponential(1)")
      )
    ),
    init = list(
      imm_cd = list(model = model, formula = re_formula, data = imm_cd_data()),
      imm_cd_bsc = list(model = cdimm("bsc"), formula = bmf(kappa ~ 1 + (1 | ID), c ~ 1 + (1 | ID),
                                                            s ~ 1 + (1 | ID)),
                        data = imm_cd_data())
    ),
    links = list(
      imm_cd = list(model = model, settable = c("mu", "kappa", "a", "c", "s")),
      imm_cd_bsc = list(model = cdimm("bsc"), settable = c("mu", "kappa", "c", "s"))
    ),
    simulate = list(
      imm_cd_full = simulate_case(
        "full", list(kappa = 2, c = 0.5, a = -0.5, s = 0.3, criterion = 0.3),
        bmf(kappa ~ 1, c ~ 1, a ~ 1, s ~ 1, criterion ~ 1)
      ),
      imm_cd_bsc = simulate_case(
        "bsc", list(kappa = 2, c = 0.5, s = 0.3), bmf(kappa ~ 1, c ~ 1, s ~ 1)
      ),
      imm_cd_abc = simulate_case(
        "abc", list(kappa = 2, c = 0.5, a = -0.5), bmf(kappa ~ 1, c ~ 1, a ~ 1)
      )
    )
  )
}
