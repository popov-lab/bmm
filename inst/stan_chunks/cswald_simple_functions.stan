// log-PDF of the censored shifted Wald model
real cswald_lpdf(real rt, real mu, real drift, real bound, real ndt, real s, real sndt, int response) {
  if (response == 1) {
    return swald_sndt_lpdf(rt | drift, bound, ndt, sndt, s);
  } else {
    return swald_sndt_lccdf(rt | drift, bound, ndt, sndt, s);
  }
}

// vectorized overload used as a loop = FALSE family when sndt is fixed at 0:
// returns the summed log-likelihood over all observations. Observations are
// split by decision; correct responses take a fully vectorized closed-form
// shifted Wald density (no CDF calls), censored responses vectorize the
// z-score preparation and keep a scalar loop only for the std_normal_lcdf
// pairs (Stan has no elementwise log-CDF; log(Phi(z)) underflows for
// z < ~-37). With sndt > 0 the convolution is lccdf-bound and the scalar
// per-observation path is faster, so the general path simply delegates to
// the scalar helpers (and configure_model keeps loop = TRUE in that case)
real cswald_lpdf(vector rt, vector mu, vector drift, vector bound,
                 vector ndt, vector s, vector sndt, array[] int dec) {
  int N = rows(rt);
  real lp = 0;
  if (min(sndt) < 0) return negative_infinity();
  if (max(sndt) >= 1e-8) {
    for (n in 1:N) {
      if (dec[n] == 1) {
        lp += swald_sndt_lpdf(rt[n] | drift[n], bound[n], ndt[n], sndt[n], s[n]);
      } else {
        lp += swald_sndt_lccdf(rt[n] | drift[n], bound[n], ndt[n], sndt[n], s[n]);
      }
    }
    return lp;
  }
  array[N] int idx1;
  array[N] int idx0;
  int n1 = 0;
  int n0 = 0;
  for (n in 1:N) {
    if (dec[n] == 1) {
      n1 += 1;
      idx1[n1] = n;
    } else {
      n0 += 1;
      idx0[n0] = n;
    }
  }
  if (n1 > 0) {
    // upper-boundary responses: vectorized shifted Wald log-density; a single
    // rt <= ndt makes the summed target -Inf, exactly as the scalar sum would
    array[n1] int ii = idx1[1:n1];
    vector[n1] t = rt[ii] - ndt[ii];
    if (min(t) <= 0) return negative_infinity();
    vector[n1] b = bound[ii];
    vector[n1] sg = s[ii];
    vector[n1] resid = b - drift[ii] .* t;
    lp += sum(log(b))
          - 0.5 * (n1 * log(2 * pi()) + 2 * sum(log(sg)) + 3 * sum(log(t)))
          - 0.5 * sum(square(resid) ./ (square(sg) .* t));
  }
  if (n0 > 0) {
    // censored responses: log survivor; elements with rt <= ndt have
    // survival 1 and contribute 0 (the process has not started yet)
    array[n0] int jj = idx0[1:n0];
    vector[n0] t0 = rt[jj] - ndt[jj];
    array[n0] int keep;
    int nk = 0;
    for (k in 1:n0) {
      if (t0[k] > 0) {
        nk += 1;
        keep[nk] = k;
      }
    }
    if (nk > 0) {
      array[nk] int kk = keep[1:nk];
      array[nk] int oo = jj[kk];
      vector[nk] t = t0[kk];
      vector[nk] b = bound[oo];
      vector[nk] d = drift[oo];
      vector[nk] sg = s[oo];
      vector[nk] denom = sg .* sqrt(t);
      vector[nk] dxt = d .* t;
      vector[nk] z1 = (dxt - b) ./ denom;
      vector[nk] z2 = -(dxt + b) ./ denom;
      vector[nk] log_c = 2 * (b .* d) ./ square(sg);
      for (k in 1:nk) {
        lp += log_diff_exp(std_normal_lcdf(-z1[k] | ),
                           log_c[k] + std_normal_lcdf(z2[k] | ));
      }
    }
  }
  return lp;
}
