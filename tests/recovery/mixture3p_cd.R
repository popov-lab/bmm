# Recovery plan for mixture3p_cd(), version "simple" with the criterion fixed
# at 0. The generator runs the observer forward and shares no code with bmm:
# it builds a display per trial (set sizes 2 to 4; the probe shows the target,
# one of the non-targets or a new feature, a third each, because the swap
# weight is identified only by intrusion probes), draws which item is
# retrieved and the remembered feature (its own von Mises sampler), computes
# Lin and Oberauer's (2022) log-likelihood ratio of "change" against "same" at
# that feature, and answers "change" when it exceeds the criterion, bounded to
# the trial's range of evidence as the model bounds it. The observer's prior
# p_s is the probability of retrieving the target. thetat and thetant arrive
# on the link scale (softmax against guessing), the others on the natural one.
mixture3p_cd_generator <- function(pars, n_trials, model) {
  rvm <- function(n, kappa) {
    out <- numeric(0)
    while (length(out) < n) {
      x <- stats::runif(n, -pi, pi)
      out <- c(out, x[stats::runif(n) < exp(kappa * (cos(x) - 1))])
    }
    out[seq_len(n)]
  }
  wrap <- function(x) (x + pi) %% (2 * pi) - pi
  vm_2pi <- function(cos_d, kappa) exp(kappa * (cos_d - 1)) / besselI(kappa, 0, TRUE)

  nt_names <- model$other_vars$nt_features
  max_nt <- length(nt_names)
  set_size <- sample(2:(max_nt + 1), n_trials, replace = TRUE)
  target <- stats::runif(n_trials, -pi, pi)
  nt_abs <- matrix(stats::runif(n_trials * max_nt, -pi, pi), n_trials, max_nt)
  nt_abs[col(nt_abs) >= set_size] <- NA
  type <- sample(c("same", "intrusion", "new"), n_trials, replace = TRUE)
  intruder <- ceiling(stats::runif(n_trials) * (set_size - 1))
  probe <- target
  probe[type == "intrusion"] <- nt_abs[cbind(which(type == "intrusion"), intruder[type == "intrusion"])]
  probe[type == "new"] <- stats::runif(sum(type == "new"), -pi, pi)

  weights <- exp(c(pars$thetat, pars$thetant, 0))
  weights <- weights / sum(weights)
  u <- stats::runif(n_trials)
  location <- target + pars$mu
  swapped <- which(u >= weights[1] & u < weights[1] + weights[2])
  which_nt <- ceiling(stats::runif(length(swapped)) * (set_size[swapped] - 1))
  location[swapped] <- nt_abs[cbind(swapped, which_nt)]
  remembered <- u < weights[1] + weights[2]
  x <- stats::runif(n_trials, -pi, pi)
  x[remembered] <- location[remembered] + rvm(sum(remembered), pars$kappa)

  p_s <- weights[1]
  llr <- -log(1 - p_s + p_s * vm_2pi(cos(x - probe), pars$kappa))
  llr_max <- -log(1 - p_s + p_s * vm_2pi(-1, pars$kappa))
  llr_min <- -log(1 - p_s + p_s * vm_2pi(1, pars$kappa))
  criterion <- pars$criterion
  bounded <- if (criterion > 0) {
    llr_max * tanh(criterion / llr_max)
  } else if (criterion < 0) {
    -llr_min * tanh(criterion / -llr_min)
  } else {
    0
  }

  out <- data.frame(as.integer(llr > bounded), wrap(probe), wrap(target))
  names(out) <- unlist(model$resp_vars[c("response", "probe", "target")], use.names = FALSE)
  out[nt_names] <- wrap(nt_abs - target)
  out[[model$other_vars$set_size]] <- set_size
  out
}

bmmtools::validation_plan(
  bmm::mixture3p_cd(
    response = "change", probe = "probe", target = "target",
    nt_features = c("nt1", "nt2", "nt3"), set_size = "set_size"
  ),
  # link scale: kappa e^1.4 = 4.1, and the softmax of (1.2, -0.6, 0) gives
  # P(target) 0.68, P(non-target) 0.11, P(guess) 0.21; away from the centres
  # of the default priors (log kappa 2, thetat and thetant 0)
  pars = c(kappa = 1.4, thetat = 1.2, thetant = -0.6),
  sds = c(kappa = 0.3, thetat = 0.5, thetant = 0.5),
  grid = expand.grid(n_subjects = c(20, 40), n_trials = c(300, 600)),
  generator = mixture3p_cd_generator,
  seed = 2026
)
