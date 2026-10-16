# the cases the shared tests run for m3; see helper-model-cases.R
model_cases_m3 <- function() {
  cats <- c("corr", "other", "dist", "npl")
  formula <- bmf(
    corr ~ b + a + c, other ~ b + a, dist ~ b + d, npl ~ b,
    c ~ 1, a ~ 1, d ~ 1
  )
  links <- list(c = "log", a = "log", d = "log")
  stored <- function(num_options) {
    list(
      model = m3(resp_cats = cats, num_options = num_options, choice_rule = "simple", links = links),
      formula = formula, data = oberauer_lewandowsky_2019_e1
    )
  }
  ss <- function(rule) {
    m3(resp_cats = c("corr", "other", "npl"), num_options = c("n_corr", "n_other", "n_npl"),
       choice_rule = rule, version = "ss")
  }
  custom <- m3(resp_cats = c("corr", "other", "npl"), num_options = c("n_corr", "n_other", "n_npl"),
               choice_rule = "softmax", version = "custom", links = list(c = "log", a = "logit"))
  re_formula <- bmf(c ~ 1 + (1 | ID), a ~ 1 + (1 | ID))
  sd <- c(c = "exponential(1)", a = "exponential(1)")
  list(
    stored_frame = list(
      m3 = stored(paste0("n_", cats)),
      m3_num_options = stored(c(opt_corr = 1, opt_other = 4, opt_dist = 5, opt_npl = 5)),
      m3_num_options_by_category = stored(c(other = 4, corr = 1, npl = 5, dist = 5)),
      m3_char_num_options_by_category = stored(
        c(other = "n_other", corr = "n_corr", npl = "n_npl", dist = "n_dist")
      )
    ),
    sd_default = list(
      m3_simple = list(model = ss("simple"), formula = re_formula,
                       data = oberauer_lewandowsky_2019_e1, sd = sd),
      m3_softmax = list(model = ss("softmax"), formula = re_formula,
                        data = oberauer_lewandowsky_2019_e1, sd = sd),
      # the custom version derives its defaults from the link of each parameter
      m3_custom = list(
        model = custom,
        formula = re_formula + bmf(corr ~ b + a + c, other ~ b + a, npl ~ b),
        data = oberauer_lewandowsky_2019_e1, sd = sd
      )
    ),
    init = list(
      m3 = list(model = ss("simple"), formula = re_formula, data = oberauer_lewandowsky_2019_e1)
    ),
    links = list(
      m3 = list(model = ss("simple"), settable = c("c", "a")),
      m3_custom = list(model = custom, settable = NULL)
    ),
    fixture = list(m3 = list(file = "bmmfit_m3_ppcheck.rds"))
  )
}
