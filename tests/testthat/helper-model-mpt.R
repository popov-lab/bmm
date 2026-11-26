# the cases the shared tests run for mpt; see helper-model-cases.R
model_cases_mpt <- function() {
  model <- mpt(mpt_2htm_trees(), tree_id = "item_type")
  trial_model <- mpt(mpt_2htm_trees(), tree_id = "item_type", response = "resp")
  re_formula <- bmf(D ~ 1 + (1 | id), g ~ 1 + (1 | id))
  list(
    # the shared check_data() test names every constructor argument, and
    # response would make mpt() read trial rows
    check_data = list(
      mpt = list(model = model, data = toy_mpt_data()),
      mpt_trial_level = list(model = trial_model, data = mpt_trial_data())
    ),
    stored_frame = list(
      # the label column is rebuilt from the one-hot response matrix
      mpt_trial_level = list(model = trial_model, formula = bmf(D ~ 1 + (1 | id), g ~ 1),
                             data = mpt_trial_data()),
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
    links = list(mpt = list(model = model, settable = c("D", "g"))),
    fixture = list(mpt = list(file = "bmmfit_mpt_ppcheck.rds")),
    # mpt() without trees has no parameters and no links
    specimen = list(mpt = model),
    simulate = list(
      # one row per trial: the design has no n_trials
      mpt_trial_level = list(
        model = mpt(mpt_2htm_trees(), tree_id = "item_type", response = "resp"),
        pars = list(D = 0.5, g = -0.2), response = "resp",
        design = data.frame(id = rep(1:2, each = 10), item_type = rep(c("old", "new"), 10)),
        formula = bmf(D ~ 1 + (1 | id), g ~ 1), simulate_pars = list(D = 0.5, g = -0.2)
      ),
      mpt = list(
        model = model, pars = list(D = 0.5, g = -0.2),
        design = data.frame(item_type = rep(c("old", "new"), 10), n_trials = rep(c(30, 50), each = 10)),
        formula = bmf(D ~ 1, g ~ 1), simulate_pars = list(D = 0.5, g = -0.2)
      ),
      # the simplex members are probabilities on both scales
      mpt_simplex = list(
        model = mpt(mpt_tree("t", list(A = "gA", B = "gB", C = "gC")), simplex = c("gA", "gB", "gC")),
        pars = list(gA = 0.2, gB = 0.3, gC = 0.5), design = data.frame(n_trials = rep(40, 10)),
        formula = bmf(gA ~ 1, gB ~ 1), simulate_pars = list(gA = 0.2, gB = 0.3, gC = 0.5)
      ),
      mpt_covariates = list(
        model = mpt(mpt_tree("main", list(correct = "D + (1 - D) * Gcorr", incorrect = "(1 - D) * (1 - Gcorr)")),
                    covariates = "Gcorr"),
        pars = list(D = 0.3), design = data.frame(Gcorr = rep(c(0.25, 0.5), 5), n_trials = 40),
        formula = bmf(D ~ 1), simulate_pars = list(D = 0.3)
      )
    )
  )
}
