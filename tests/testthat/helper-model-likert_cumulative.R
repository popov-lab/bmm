# Likert responses drawn from the latent variable directly, as the recovery
# generators do: activation = person + item location, normal noise, and the
# K - 1 thresholds written out here rather than built by the package
likert_sim_frame <- function(n_id = 6, n_item = 4, K = 5L, seed = 662) {
  withr::local_seed(seed)
  thr <- stats::qlogis(seq_len(K - 1) / K)
  dat <- expand.grid(item = paste0("i", seq_len(n_item)), id = seq_len(n_id),
                     stringsAsFactors = FALSE)
  theta <- stats::rnorm(n_id, 0, 0.8)
  beta <- seq(-0.5, 0.5, length.out = n_item)
  latent <- theta[dat$id] + beta[match(dat$item, unique(dat$item))] + stats::rnorm(nrow(dat))
  dat$y <- 1L + vapply(latent, function(x) sum(x > thr), integer(1))
  dat$keying <- ifelse(match(dat$item, unique(dat$item)) %% 2 == 0, -1, 1)
  dat
}

# the cases the shared tests run for likert_cumulative(); see helper-model-cases.R
model_cases_likert_cumulative <- function() {
  dat <- likert_sim_frame()
  # every category observed, so check_data() raises no warning
  dat$y[1:5] <- 1:5
  ordered_dat <- dat
  ordered_dat$y <- factor(c("sd", "d", "n", "a", "sa")[dat$y],
                          levels = c("sd", "d", "n", "a", "sa"), ordered = TRUE)
  model <- likert_cumulative("y", "item", n_categories = 5)
  keyed <- likert_cumulative("y", "item", keying = "keying", threshold_type = "parsimonious")
  formula <- bmf(activation ~ 0 + item + (1 | id))
  re_formula <- bmf(activation ~ 0 + item + (1 | id), delta1 ~ 1 + (1 | id),
                    delta2 ~ 1, delta3 ~ 1)
  list(
    stored_frame = list(
      likert_cumulative = list(model = model, formula = formula, data = dat),
      likert_cumulative_ordered = list(
        model = keyed, formula = bmf(activation ~ 0 + item + (1 | id), spacing ~ 1 + (1 | id)),
        data = ordered_dat
      )
    ),
    sd_default = list(
      likert_cumulative = list(
        model = model, formula = re_formula, data = dat,
        sd = c(activation = "exponential(1)", delta1 = "exponential(2)")
      )
    ),
    init = list(
      likert_cumulative = list(model = model, formula = re_formula, data = dat)
    ),
    check_data = list(likert_cumulative = list(
      model = model, data = dat, formula = formula
    )),
    links = list(likert_cumulative = list(model = model, settable = "activation")),
    specimen = list(likert_cumulative = model),
    simulate = list(
      likert_cumulative = list(
        model = model,
        pars = list(activation = 0.3, delta1 = 0, delta2 = -0.1, delta3 = 0.2),
        # 400 rows leave the rarest category (p 0.04) unobserved with
        # probability 1e-7, which check_data() would warn about
        design = expand.grid(item = paste0("i", 1:4), id = 1:100),
        formula = bmf(activation ~ 1, delta1 ~ 1, delta2 ~ 1, delta3 ~ 1),
        simulate_pars = list(activation = 0.3, delta1 = 0, delta2 = -0.1, delta3 = 0.2)
      )
    )
  )
}
