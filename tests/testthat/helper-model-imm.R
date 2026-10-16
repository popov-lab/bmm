# the cases the shared tests run for imm; see helper-model-cases.R
model_cases_imm <- function() {
  model <- imm("dev_rad",
    nt_features = paste0("col_nt", 1:7), nt_distances = paste0("dist_nt", 1:7),
    set_size = "set_size"
  )
  re_formula <- bmf(kappa ~ 1 + (1 | ID), a ~ 1 + (1 | ID), c ~ 1 + (1 | ID), s ~ 1 + (1 | ID))
  list(
    stored_frame = list(
      imm = list(model = model, formula = bmf(c ~ 1, a ~ 1, s ~ 1, kappa ~ 1), data = small_data())
    ),
    sd_default = list(
      imm = list(
        model = model, formula = re_formula, data = small_data(),
        sd = c(kappa = "exponential(1)", a = "exponential(1)", c = "exponential(1)",
               s = "exponential(1)")
      )
    ),
    init = list(imm = list(model = model, formula = re_formula, data = small_data())),
    links = list(imm = list(model = model, settable = character(0)))
  )
}
