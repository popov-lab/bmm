# Exact posterior draws for conjugate normal designs with known variance
# components. These are the oracles for the reliability tests: every target
# below has a closed form, so the tests do not depend on brms or Stan.

# person-only design: theta_i ~ N(0, tau2), ybar_i ~ N(theta_i, sigma2 / K_i)
conjugate_person_draws <- function(N, tau2, sigma2, K, ndraws) {
  K <- rep_len(K, N)
  theta <- stats::rnorm(N, 0, sqrt(tau2))
  ybar <- stats::rnorm(N, theta, sqrt(sigma2 / K))
  lambda <- tau2 / (tau2 + sigma2 / K)
  m <- lambda * ybar
  e2 <- lambda * sigma2 / K
  draws <- matrix(stats::rnorm(N * ndraws, rep(m, each = ndraws), rep(sqrt(e2), each = ndraws)), ndraws, N)
  list(draws = draws, theta = theta, m = m, e2 = e2, K = K, lambda = lambda)
}

# reliability of the posterior means for the person-only design, from the
# closed form of the posterior rather than from simulated true scores
conjugate_person_reliability <- function(tau2, sigma2, K) {
  1 - mean(tau2 * (sigma2 / K) / (tau2 + sigma2 / K)) / tau2
}

# person x facet design: theta_po = u_p + u_po, ybar_po ~ N(theta_po, s2).
# Returns draws of the universe score u_p, of the cell scores theta_po (one
# matrix per facet level) and the closed-form G coefficient.
conjugate_facet_draws <- function(N, sp2, spo2, s2, n_o, ndraws) {
  Z <- cbind(1, diag(n_o))
  V <- solve(crossprod(Z) / s2 + diag(1 / c(sp2, rep(spo2, n_o))))
  L <- chol(V)
  universe <- matrix(0, ndraws, N)
  cells <- replicate(n_o, matrix(0, ndraws, N), simplify = FALSE)
  for (i in seq_len(N)) {
    u <- c(stats::rnorm(1, 0, sqrt(sp2)), stats::rnorm(n_o, 0, sqrt(spo2)))
    ybar <- u[1] + u[-1] + stats::rnorm(n_o, 0, sqrt(s2))
    mean_post <- V %*% (crossprod(Z, ybar) / s2)
    d <- matrix(stats::rnorm(ndraws * (n_o + 1)), ndraws) %*% L
    d <- sweep(d, 2, mean_post, `+`)
    universe[, i] <- d[, 1]
    for (o in seq_len(n_o)) cells[[o]][, i] <- d[, 1] + d[, o + 1]
  }
  list(
    universe = universe, cells = cells,
    g_true = sp2 / (sp2 + (spo2 + s2) / n_o)
  )
}
