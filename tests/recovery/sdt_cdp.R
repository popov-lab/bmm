# Recovery plan for sdt_cdp() (continuous dual-process, Wixted & Mickes 2010)
# on one row per trial (ADR 0006): the remember / know / new judgment and its
# confidence on the unified scale, with 4 "new" and 2 "old" confidence levels
# and the default "parsimonious" thresholds. The generators run the observer
# forward and share no code with bmm: each trial draws a familiarity F and a
# recollection R from independent normals with SD 1 (the defaults sigmar = 0,
# rho = 0), with means 0 on new items and dfam and drec on old items. The
# strength is S = F + R, the confidence is 1 plus the number of the K - 1
# thresholds S exceeds (K = n_new + n_old), and the judgment is "new" for
# confidence up to n_new, otherwise "remember" when R exceeds rcrit and "know"
# when it does not. The thresholds lie on the raw S axis at
# criterion + exp(spacing) * (log(k / (K - k)) - log(n_new / (K - n_new))), so
# that threshold n_new is the criterion. bmm computes the category
# probabilities from bivariate-normal masses, so an error in dsdt_cdp(),
# rsdt_cdp() or the Stan likelihood does not cancel.
# Unequal on purpose: with as many "new" as "old" levels the shift
# log(n_new / (K - n_new)) is 0, and a likelihood that left it out or took it
# from n_old would fit these data exactly as well.
sdt_cdp_n_new <- 4L
sdt_cdp_n_old <- 2L

# judgment and confidence; one value of each parameter per trial
sdt_cdp_respond <- function(stimulus, dfam, drec, criterion, spacing, rcrit,
                            n_new = sdt_cdp_n_new, n_old = sdt_cdp_n_old) {
  n <- length(stimulus)
  k <- n_new + n_old
  familiarity <- stats::rnorm(n, mean = dfam * stimulus, sd = 1)
  recollection <- stats::rnorm(n, mean = drec * stimulus, sd = 1)
  j <- seq_len(k - 1L)
  canonical <- log(j / (k - j)) - log(n_new / (k - n_new))
  thresholds <- rep_len(criterion, n) + outer(rep_len(exp(spacing), n), canonical)
  confidence <- 1L + as.integer(rowSums((familiarity + recollection) > thresholds))
  judgment <- ifelse(
    confidence <= n_new, "new", ifelse(recollection > rcrit, "remember", "know")
  )
  data.frame(judgment, confidence, stringsAsFactors = FALSE)
}

# one subject's trials; the stimulus alternates so that new and old items are
# balanced (the stimulus is a design variable and is part of the generated data)
sdt_cdp_generator <- function(pars, n_trials, model) {
  stimulus <- rep_len(c(0L, 1L), n_trials)
  out <- sdt_cdp_respond(
    stimulus, pars$dfam, pars$drec, pars$criterion, pars$spacing, pars$rcrit,
    model$other_vars$n_new, model$other_vars$n_old
  )
  out <- cbind(stimulus, out)
  names(out) <- c(
    model$other_vars$stimulus, model$resp_vars$judgment, model$resp_vars$confidence
  )
  out
}

# the calibration's data sets: one prior draw under the fit's names and the
# design data (from sdt_cdp_generator), the judgments and confidences drawn
# again from the subject values of that draw
sdt_cdp_sbc_generator <- function(draws, data) {
  value <- function(par) {
    draws[[sprintf("b_%s_Intercept", par)]] +
      unname(draws[sprintf("r_id__%s[%s,Intercept]", par, as.character(data$id))])
  }
  out <- sdt_cdp_respond(
    data$stimulus, value("dfam"), value("drec"), value("criterion"),
    value("spacing"), value("rcrit")
  )
  data$judgment <- out$judgment
  data$confidence <- out$confidence
  data
}

bmmtools::validation_plan(
  bmm::sdt_cdp(
    stimulus = "stimulus", judgment = "judgment", confidence = "confidence",
    n_new = sdt_cdp_n_new, n_old = sdt_cdp_n_old
  ),
  # dfam = 1.4 and drec = 1.7 sit 0.4 and 0.7 prior SDs above the prior centre
  # (1, SD 1); criterion = 1.3 on the strength axis is 0.9 prior SDs above its
  # centre (0, SD 1.5), near the midpoint of the lure and target strength means
  # (1.55); spacing = 0.2 sits 0.4 prior SDs above its centre (0, SD 0.5);
  # rcrit = 1.1 sits 1.1 prior SDs above its centre (0, SD 1) and between the
  # lure and target recollection means (0 and 1.7). sigmar, rho and kcrit keep
  # their fixed values (equal recollection variance, independent F and R, no
  # Know/Guess split), so the data hold no "guess" judgment.
  pars = c(dfam = 1.4, drec = 1.7, criterion = 1.3, spacing = 0.2, rcrit = 1.1),
  # all below the prior's mean SD (exponential(1) for dfam and drec,
  # exponential(2), mean 0.5, for the rest)
  sds = c(dfam = 0.3, drec = 0.3, criterion = 0.3, spacing = 0.15, rcrit = 0.3),
  # bmmtools' default grid plus 600 trials: eight response categories and five
  # free parameters per subject; the 600-trial rows are the costly ones (about
  # 30000 rows per fit at 50 subjects)
  grid = expand.grid(n_subjects = c(20L, 50L), n_trials = c(30L, 100L, 300L, 600L)),
  generator = sdt_cdp_generator,
  # SBC at 100 trials per subject rather than bmmtools' 50, as for sdt_rating
  sbc = list(
    n_trials = 100L, generator = sdt_cdp_sbc_generator, cores_per_fit = 4L,
    sampler = list(chains = 4L, iter = 1000L, warmup = 500L, control = list(adapt_delta = 0.95))
  ),
  seed = 2026
)
