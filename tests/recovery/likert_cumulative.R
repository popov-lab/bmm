# Recovery plan for likert_cumulative(): persons answer the same 10 five-point
# items, whose locations are known and enter as an offset, so the plan scores
# the person trait (activation) and the person-specific threshold spread
# (spacing, the extreme response style) against their true values.
#
# The generator draws the latent response variable itself, the person's
# activation plus the item location plus noise from the named distribution,
# and counts the thresholds below it, written out here for the parsimonious
# type; it does not call rlikert_cumulative() or bmm's threshold builders.
#
# The noise distribution is set by BMM_RECOVERY_DIST ("normal", the default,
# or "logistic"); each run is its own card.
likert_dist <- Sys.getenv("BMM_RECOVERY_DIST", "normal")
n_items <- 10L

bmmtools::validation_plan(
  bmm::likert_cumulative("y", "item", n_categories = 5, dist = likert_dist,
                         threshold_type = "parsimonious"),
  pars = c(activation = 0.2, spacing = 0),
  sds = c(activation = 1, spacing = 0.3),
  grid = data.frame(n_subjects = c(100L, 200L), n_trials = n_items),
  reps = 10L,
  formula = bmm::bmf(activation ~ 1 + offset(item_loc) + (1 | p | id),
                     spacing ~ 1 + (1 | p | id)),
  # the model counts item among its response columns, so the generator returns
  # it, in the order of the design's rows
  trial_design = data.frame(item_loc = seq(-1.2, 1.2, length.out = n_items)),
  # SBC hands its generator one prior draw under brms's draw names; the r_
  # draws are deviations from the population value
  sbc = list(
    n_subjects = 50L, n_trials = n_items, n_sims = 100L,
    generator = function(draws, data) {
      K <- 5L
      person <- function(par) {
        draws[[paste0("b_", par, "_Intercept")]] +
          draws[paste0("r_id__", par, "[", data$id, ",Intercept]")]
      }
      spacing <- person("spacing")
      latent <- person("activation") + data$item_loc + switch(likert_dist,
        normal = stats::rnorm(nrow(data)),
        logistic = stats::rlogis(nrow(data))
      )
      data$y <- vapply(seq_len(nrow(data)), function(i) {
        thresholds <- exp(spacing[i]) * log(seq_len(K - 1) / (K - seq_len(K - 1)))
        1L + sum(latent[i] > thresholds)
      }, integer(1))
      data
    }
  ),
  prior_check = list(n_subjects = 50L, n_trials = n_items),
  generator = function(pars, n_trials, model, trial_design) {
    K <- 5L
    thresholds <- exp(pars$spacing) * log(seq_len(K - 1) / (K - seq_len(K - 1)))
    noise <- switch(model$other_vars$dist,
      normal = stats::rnorm(n_trials),
      logistic = stats::rlogis(n_trials)
    )
    latent <- pars$activation + trial_design$item_loc + noise
    data.frame(y = 1L + vapply(latent, function(x) sum(x > thresholds), integer(1)),
               item = paste0("i", seq_len(n_trials)))
  },
  seed = 2026
)
