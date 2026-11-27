# Recovery plan for sdt_rating() version "standard" on one row per trial
# (ADR 0006), with the default "parsimonious" thresholds and the default
# criterion parameter (bias is the alternative member of its group and is not
# named). The generators run the observer forward and share no code with bmm:
# each trial draws one normal evidence value (mean 0 on noise trials, mean d on
# signal trials, SD 1; sdratio keeps its fixed value of equal variance), and the
# rating is 1 plus the number of the K - 1 thresholds the evidence exceeds. The
# thresholds sit at criterion + exp(spacing) * log(k / (K - k)), k = 1 .. K - 1,
# the canonical positions of Selker et al. (2019) as sdt_rating() documents
# them, so that with an even K the middle threshold is the criterion. bmm
# computes the category probabilities from pnorm() differences, so an error in
# dsdt_rating(), rsdt_rating() or the Stan likelihood does not cancel.
sdt_rating_k <- 6L

# a rating from 1 to K; one d, criterion and spacing per trial
sdt_rating_respond <- function(stimulus, d, criterion, spacing, k = sdt_rating_k) {
  evidence <- stats::rnorm(length(stimulus), mean = d * stimulus, sd = 1)
  canonical <- log(seq_len(k - 1L) / (k - seq_len(k - 1L)))
  n <- length(stimulus)
  thresholds <- rep_len(criterion, n) + outer(rep_len(exp(spacing), n), canonical)
  1L + as.integer(rowSums(evidence > thresholds))
}

# one subject's trials; the stimulus alternates so that noise and signal trials
# are balanced (the stimulus is a design variable and is part of the generated
# data)
sdt_rating_generator <- function(pars, n_trials, model) {
  stimulus <- rep_len(c(0L, 1L), n_trials)
  rating <- sdt_rating_respond(
    stimulus, pars$d, pars$criterion, pars$spacing, model$other_vars$n_ratings
  )
  out <- data.frame(stimulus, rating)
  names(out) <- c(model$other_vars$stimulus, model$resp_vars$response)
  out
}

# the calibration's data sets: one prior draw under the fit's names and the
# design data (from sdt_rating_generator), the ratings drawn again from the
# subject values of that draw
sdt_rating_sbc_generator <- function(draws, data) {
  subject <- function(par) {
    unname(draws[sprintf("r_id__%s[%s,Intercept]", par, as.character(data$id))])
  }
  data$rating <- sdt_rating_respond(
    data$stimulus,
    d = draws[["b_d_Intercept"]] + subject("d"),
    criterion = draws[["b_criterion_Intercept"]] + subject("criterion"),
    spacing = draws[["b_spacing_Intercept"]] + subject("spacing")
  )
  data
}

bmmtools::validation_plan(
  bmm::sdt_rating(
    response = "rating", stimulus = "stimulus", n_ratings = sdt_rating_k,
    version = "standard", threshold_type = "parsimonious"
  ),
  # d = 1.7 sits 0.7 prior SDs above the prior centre (1); criterion = 0.9 is
  # 0.3 SDs above its centre (0.5, SD 1.5), the unbiased point d / 2 = 0.85;
  # spacing = 0.3 (thresholds 1.35 times the canonical positions) sits 0.6
  # prior SDs above its centre (0, SD 0.5)
  pars = c(d = 1.7, criterion = 0.9, spacing = 0.3),
  # all below the prior's mean SD (exponential(1) for d and criterion,
  # exponential(2), mean 0.5, for spacing)
  sds = c(d = 0.4, criterion = 0.3, spacing = 0.15),
  # bmmtools' default grid plus 600 trials: K = 6 categories give five
  # thresholds, and at 30 trials a subject holds about 2.5 trials per category,
  # so the largest design adds the trial count the subject-level thresholds need
  grid = expand.grid(n_subjects = c(20L, 50L), n_trials = c(30L, 100L, 300L, 600L)),
  generator = sdt_rating_generator,
  # SBC at 100 trials per subject rather than bmmtools' 50, for the same reason
  sbc = list(
    n_trials = 100L, generator = sdt_rating_sbc_generator, cores_per_fit = 4L,
    sampler = list(chains = 4L, iter = 1000L, warmup = 500L, control = list(adapt_delta = 0.95))
  ),
  seed = 2026
)
