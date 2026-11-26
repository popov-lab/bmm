# Recovery plan for sdt_mafc() on one row per trial (ADR 0006). The generators
# run the observer forward and share no code with bmm: on each trial the signal
# alternative's evidence is a normal draw with mean d, each of the m - 1
# distractors' evidence a standard normal draw, and the answer is correct (1)
# when the signal's evidence is the largest. bmm computes
# P(correct) = integral of dnorm(x - d) * pnorm(x)^(m - 1) by quadrature, so an
# error in dsdt_mafc(), rsdt_mafc() or the Stan likelihood does not cancel.
sdt_mafc_m <- 4L

# 1 = correct; one d per trial
sdt_mafc_respond <- function(d, m = sdt_mafc_m) {
  n <- length(d)
  signal <- stats::rnorm(n, mean = d, sd = 1)
  distractors <- matrix(stats::rnorm(n * (m - 1L)), nrow = n)
  as.integer(signal > apply(distractors, 1L, max))
}

# one subject's trials; the model has no design variable (every trial holds
# exactly one signal alternative)
sdt_mafc_generator <- function(pars, n_trials, model) {
  out <- data.frame(sdt_mafc_respond(rep(pars$d, n_trials), model$other_vars$m))
  names(out) <- model$resp_vars$response
  out
}

# the calibration's data sets: one prior draw under the fit's names, the
# responses drawn again from the subject values of that draw
sdt_mafc_sbc_generator <- function(draws, data) {
  d <- draws[["b_d_Intercept"]] +
    unname(draws[sprintf("r_id__d[%s,Intercept]", as.character(data$id))])
  data$correct <- sdt_mafc_respond(d)
  data
}

bmmtools::validation_plan(
  bmm::sdt_mafc(response = "correct", m = sdt_mafc_m),
  # m = 4 alternatives: accuracy is the only datum, and 4 alternatives leave
  # more room between chance (0.25) and ceiling than 2 do. d = 1.6 sits 0.6
  # prior SDs above the prior centre (1); P(correct) is then about 0.6.
  pars = c(d = 1.6),
  # below the prior's mean SD (exponential(1) has mean 1)
  sds = c(d = 0.4),
  # bmmtools' default grid
  grid = expand.grid(n_subjects = c(20L, 50L), n_trials = c(30L, 100L, 300L)),
  generator = sdt_mafc_generator,
  sbc = list(
    generator = sdt_mafc_sbc_generator, cores_per_fit = 4L,
    sampler = list(chains = 4L, iter = 1000L, warmup = 500L, control = list(adapt_delta = 0.95))
  ),
  seed = 2026
)
