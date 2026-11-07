# Recovery plan for imm_cd(), version "full" with the criterion fixed at 0.
# The generator runs the observer forward and shares no code with bmm: it
# builds a display per trial (set sizes 2 to 4, each non-target at its own
# spatial distance from the target; the probe shows the target, one of the
# non-targets or a new feature, a third each, because a and s are identified
# only by intrusion probes), computes each item's activation (c exp(-s d) + a,
# the target at d = 0) against the background b, draws which wins and the
# remembered feature (its own von Mises sampler), computes Lin and Oberauer's
# (2022) log-likelihood ratio of "change" against "same" at that feature, and
# answers "change" when it exceeds the criterion, bounded to the trial's range
# of evidence as the model bounds it. The observer's prior p_s is the
# probability of retrieving the target. All parameters arrive on the natural
# scale.
imm_cd_generator <- function(pars, n_trials, model) {
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
  dist_names <- model$other_vars$nt_distances
  max_nt <- length(nt_names)
  set_size <- sample(2:(max_nt + 1), n_trials, replace = TRUE)
  inactive <- col(matrix(0, n_trials, max_nt)) >= set_size
  target <- stats::runif(n_trials, -pi, pi)
  nt_abs <- matrix(stats::runif(n_trials * max_nt, -pi, pi), n_trials, max_nt)
  distance <- matrix(stats::runif(n_trials * max_nt, 0.3, 2.5), n_trials, max_nt)
  nt_abs[inactive] <- NA
  distance[inactive] <- NA
  type <- sample(c("same", "intrusion", "new"), n_trials, replace = TRUE)
  intruder <- ceiling(stats::runif(n_trials) * (set_size - 1))
  probe <- target
  is_intrusion <- which(type == "intrusion")
  probe[is_intrusion] <- nt_abs[cbind(is_intrusion, intruder[is_intrusion])]
  probe[type == "new"] <- stats::runif(sum(type == "new"), -pi, pi)

  activation <- cbind(pars$c + pars$a, pars$c * exp(-pars$s * distance) + pars$a)
  activation[, -1][inactive] <- 0
  total <- rowSums(activation) + pars$b
  u <- stats::runif(n_trials) * total
  # an inactive non-target adds nothing to the running sum and never wins
  winner <- rowSums(t(apply(activation, 1, cumsum)) <= u) + 1
  location <- target + pars$mu
  swapped <- which(winner > 1 & winner <= max_nt + 1)
  location[swapped] <- nt_abs[cbind(swapped, winner[swapped] - 1)]
  remembered <- winner <= max_nt + 1
  x <- stats::runif(n_trials, -pi, pi)
  x[remembered] <- location[remembered] + rvm(sum(remembered), pars$kappa)

  p_s <- activation[, 1] / total
  llr <- -log(1 - p_s + p_s * vm_2pi(cos(x - probe), pars$kappa))
  llr_max <- -log(1 - p_s + p_s * vm_2pi(-1, pars$kappa))
  llr_min <- -log(1 - p_s + p_s * vm_2pi(1, pars$kappa))
  criterion <- pars$criterion
  bounded <- if (criterion > 0) {
    llr_max * tanh(criterion / llr_max)
  } else if (criterion < 0) {
    -llr_min * tanh(criterion / -llr_min)
  } else {
    rep(0, n_trials)
  }

  out <- data.frame(as.integer(llr > bounded), wrap(probe), wrap(target))
  names(out) <- unlist(model$resp_vars[c("response", "probe", "target")], use.names = FALSE)
  out[nt_names] <- wrap(nt_abs - target)
  out[dist_names] <- distance
  out[[model$other_vars$set_size]] <- set_size
  out
}

bmmtools::validation_plan(
  bmm::imm_cd(
    response = "change", probe = "probe", target = "target",
    nt_features = c("nt1", "nt2", "nt3"), nt_distances = c("dist1", "dist2", "dist3"),
    set_size = "set_size"
  ),
  # link scale: kappa e^1.4 = 4.1, c e^1.1 = 3.0, a e^-0.9 = 0.41, s e^0.5 =
  # 1.65; away from the centres of the default priors (log kappa 2, log c,
  # log a and log s 0)
  pars = c(kappa = 1.4, c = 1.1, a = -0.9, s = 0.5),
  sds = c(kappa = 0.3, c = 0.4, a = 0.4, s = 0.3),
  grid = expand.grid(n_subjects = c(20, 40), n_trials = c(300, 600)),
  generator = imm_cd_generator,
  seed = 2026
)
