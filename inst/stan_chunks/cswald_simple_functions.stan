// log-PDF of the censored shifted Wald model
real cswald_lpdf(real rt, real mu, real drift, real bound, real ndt, real s, int response) {
  if (response == 1) {
    return swald_lpdf(rt | drift, bound, ndt, s);
  } else {
    return swald_lccdf(rt | drift, bound, ndt, s);
  }
}

// vectorized overload for the loop = FALSE family: density for responses, survivor for censored trials
real cswald_lpdf(vector rt, vector mu, vector drift, vector bound,
                 vector ndt, vector s, array[] int dec) {
  int N = rows(rt);
  array[N] int idx1;
  array[N] int idx0;
  int n1 = 0;
  int n0 = 0;
  real lp = 0;

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
    array[n1] int ii = idx1[1:n1];
    lp += sum(swald_log_dens_vec(rt[ii] - ndt[ii], drift[ii], bound[ii], s[ii]));
  }

  if (n0 > 0) {
    // censored trials with rt <= ndt have survival 1, so they are dropped
    array[n0] int jj = idx0[1:n0];
    vector[n0] t0 = rt[jj] - ndt[jj];
    array[n0] int keep;
    int nk = 0;
    for (k in 1:n0) {
      if (t0[k] > 0) {
        nk += 1;
        keep[nk] = jj[k];
      }
    }
    if (nk > 0) {
      array[nk] int oo = keep[1:nk];
      lp += sum(swald_log_surv_vec(rt[oo] - ndt[oo], drift[oo], bound[oo], s[oo]));
    }
  }

  return lp;
}
