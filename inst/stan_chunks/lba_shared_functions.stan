real lba_softplus(real x) {
  return 1e-12 * log1p_exp(x / 1e-12);
}

real lba_log_positive(real x) {
  return log(lba_softplus(x));
}

real lba_race_loglik(int response, array[] int n, real log_pdf_response,
                     array[] real log_surv) {
  int n_cats = size(n);
  real log_lik = log(n[response]) + log_pdf_response;
  for (j in 1:n_cats) {
    if (j == response) {
      if (n[j] > 1) {
        log_lik += (n[j] - 1) * log_surv[j];
      }
    } else if (n[j] > 0) {
      log_lik += n[j] * log_surv[j];
    }
  }
  return log_lik;
}
