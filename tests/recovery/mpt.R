# Recovery plan for mpt(): the two-high-threshold model (2HTM) of old/new
# recognition on trial-level data. Two trees, because one old/new tree has one
# degree of freedom for D and g (bmm prints "rank 1 of 2"); targets and lures
# together identify both.
#
# The generator walks the tree with one uniform draw per branch point and
# shares no code with bmm's rmpt() or its likelihood: a target is detected
# with probability D (answer "old"), otherwise the answer is a guess; a lure is
# detected with probability D (answer "new"), otherwise guessed. bmm's
# likelihood is built from the branch expressions D + (1 - D) * g and so on,
# so a mistake in the tree-to-likelihood step does not cancel.
mpt_2htm_respond <- function(is_target, D, g) {
  n <- length(is_target)
  detected <- stats::runif(n) < D
  guessed_old <- stats::runif(n) < g
  # detection reveals the truth; without it only the guess decides
  old <- ifelse(detected, is_target, guessed_old)
  ifelse(old, "old", "new")
}

# one subject's trials; the design (half targets, half lures, interleaved)
# comes in as trial_design, so n_trials is a multiple of two
mpt_2htm_generator <- function(pars, n_trials, model, trial_design) {
  resp <- mpt_2htm_respond(trial_design$item_type == "target", pars$D, pars$g)
  out <- data.frame(resp = resp)
  # bmmtools' check_generated() asks for one column per response category;
  # bmm's check_data() builds them from `resp`, so they are filled here
  # consistently rather than left out
  out$old <- as.integer(resp == "old")
  out$new <- as.integer(resp == "new")
  out
}

# the calibration's data sets: the responses are drawn again from the subject
# values of one prior draw, under the fit's draw names (logit link)
mpt_2htm_sbc_generator <- function(draws, data) {
  subject <- function(par) {
    draws[[sprintf("b_%s_Intercept", par)]] +
      draws[sprintf("r_id__%s[%s,Intercept]", par, as.character(data$id))]
  }
  D <- stats::plogis(subject("D"))
  g <- stats::plogis(subject("g"))
  resp <- mpt_2htm_respond(data$item_type == "target", D, g)
  data$resp <- resp
  data$old <- as.integer(resp == "old")
  data$new <- as.integer(resp == "new")
  data
}

mpt_2htm <- bmm::mpt(
  list(
    bmm::mpt_tree("target", old ~ D + (1 - D) * g, new ~ (1 - D) * (1 - g)),
    bmm::mpt_tree("lure", old ~ (1 - D) * g, new ~ D + (1 - D) * (1 - g))
  ),
  tree_id = "item_type", response = "resp"
)

bmmtools::validation_plan(
  mpt_2htm,
  # default priors are logistic(0, 1) on both logits, centred on p = 0.5:
  # D = 0.7 (logit 0.85, 0.85 scale units or 0.47 prior SDs above the centre) and
  # g = 0.3 (logit -0.85, the same distance below), so neither sits at the
  # centre and the guessing rate is below the unbiased 0.5
  pars = c(D = stats::qlogis(0.7), g = stats::qlogis(0.3)),
  # person SDs on the logit scale: ordinary between-person spread in
  # recognition, and below the prior's exponential(1) mean of 1
  sds = c(D = 0.6, g = 0.5),
  # bmmtools' default grid: a trial is one item, so 30 trials is 15 targets and
  # 15 lures, enough for the card's read-off of the smallest identified design
  trial_design = function(n_trials) {
    data.frame(item_type = rep(c("target", "lure"), length.out = n_trials))
  },
  generator = mpt_2htm_generator,
  sbc = list(
    generator = mpt_2htm_sbc_generator, cores_per_fit = 4L,
    sampler = list(chains = 4L, iter = 1000L, warmup = 500L, control = list(adapt_delta = 0.95))
  ),
  seed = 2026
)
