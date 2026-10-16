# the cases the shared tests run for mpt; see helper-model-cases.R
model_cases_mpt <- function() {
  model <- mpt(mpt_2htm_trees(), tree_id = "item_type")
  re_formula <- bmf(D ~ 1 + (1 | id), g ~ 1 + (1 | id))
  list(
    stored_frame = list(
      # brms drops the tree column unless a formula names it, so it is rebuilt
      mpt = list(model = model, formula = bmf(D ~ 1 + (1 | id), g ~ 1), data = toy_mpt_data()),
      mpt_tree_predictor = list(model = model, formula = bmf(D ~ 1 + item_type, g ~ 1),
                                data = toy_mpt_data()),
      mpt_single_tree = list(
        model = mpt(mpt_tree("study", list(
          C = "cp + (1 - cp) * rp * rp",
          E = "2 * (1 - cp) * rp * (1 - rp)",
          U = "(1 - cp) * (1 - rp) * (1 - rp)"
        ))),
        formula = bmf(cp ~ 1 + (1 | id), rp ~ 1),
        data = data.frame(id = factor(1:6), C = 30:35, E = 20L, U = 10L)
      ),
      mpt_impossible = list(
        model = mpt(mpt_impossible_trees(), tree_id = "tree"),
        formula = bmf(Pm ~ 1 + (1 | id), Pb ~ 1),
        data = mpt_impossible_data()
      ),
      mpt_covariates = list(
        model = mpt(
          mpt_tree("main", list(
            correct = "D + (1 - D) * Gcorr",
            incorrect = "(1 - D) * (1 - Gcorr)"
          )),
          covariates = "Gcorr"
        ),
        formula = bmf(D ~ 1 + (1 | id)),
        data = data.frame(
          id = factor(1:10), Gcorr = rep(c(0.25, 0.5), 5),
          correct = 21:30, incorrect = 19:10
        )
      )
    ),
    sd_default = list(
      mpt = list(model = model, formula = re_formula, data = toy_mpt_data(),
                 sd = c(D = "exponential(1)", g = "exponential(1)"))
    ),
    init = list(mpt = list(model = model, formula = re_formula, data = toy_mpt_data())),
    links = list(mpt = list(model = model, settable = c("D", "g")))
  )
}
