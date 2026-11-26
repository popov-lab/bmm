# Recovery plan for sdt_yn() on one row per trial (ADR 0006). The generators
# run the observer forward and share no code with bmm: each trial draws one
# normal evidence value (mean 0 on noise trials, mean d on signal trials, SD 1,
# the equal-variance default) and the observer says "old" (1) when the evidence
# exceeds the criterion. bmm computes P(old) = 1 - pnorm(criterion - d * stimulus)
# instead, so an error in dsdt_yn(), rsdt_yn() or the Stan likelihood does not
# cancel.

# 1 = "old"; one d and criterion per trial, so a subject's values recycle
sdt_yn_respond <- function(stimulus, d, criterion) {
  evidence <- stats::rnorm(length(stimulus), mean = d * stimulus, sd = 1)
  as.integer(evidence > criterion)
}

# one subject's trials, for the recovery grid; the stimulus alternates so that
# noise and signal trials are balanced (the stimulus is a design variable and
# is part of the generated data)
sdt_yn_generator <- function(pars, n_trials, model) {
  stimulus <- rep_len(c(0L, 1L), n_trials)
  out <- data.frame(stimulus, sdt_yn_respond(stimulus, pars$d, pars$criterion))
  names(out) <- c(model$other_vars$stimulus, model$resp_vars$response)
  out
}

# the calibration's data sets: bmmtools hands over one prior draw under the
# fit's names and the design data (from sdt_yn_generator), and the responses are
# drawn again from the subject values of that draw
sdt_yn_sbc_generator <- function(draws, data) {
  subject <- function(par) {
    unname(draws[sprintf("r_id__%s[%s,Intercept]", par, as.character(data$id))])
  }
  d <- draws[["b_d_Intercept"]] + subject("d")
  criterion <- draws[["b_criterion_Intercept"]] + subject("criterion")
  data$response <- sdt_yn_respond(data$stimulus, d, criterion)
  data
}

bmmtools::validation_plan(
  bmm::sdt_yn(response = "response", stimulus = "stimulus"),
  # d = 1.8 sits 0.8 prior SDs above the prior centre (1); criterion = 1.2 is
  # 0.5 SDs above its centre (0.5, SD 1.5) and conservative, since the unbiased
  # observer sits at d / 2 = 0.9. The default criterion member of the
  # criterion / bias group is estimated; bias is not named.
  pars = c(d = 1.8, criterion = 1.2),
  # both below the prior's mean SD (exponential(1) has mean 1)
  sds = c(d = 0.4, criterion = 0.3),
  # bmmtools' default grid: a trial row is cheap, and 30 trials (15 per
  # stimulus class) still pin d and criterion once subjects share a population
  grid = expand.grid(n_subjects = c(20L, 50L), n_trials = c(30L, 100L, 300L)),
  generator = sdt_yn_generator,
  # cores per fit, one per chain; sampler budget as tcc_cd needed it
  sbc = list(
    generator = sdt_yn_sbc_generator, cores_per_fit = 4L,
    sampler = list(chains = 4L, iter = 1000L, warmup = 500L, control = list(adapt_delta = 0.95))
  ),
  seed = 2026
)
