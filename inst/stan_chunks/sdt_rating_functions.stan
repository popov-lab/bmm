// parsimonious thresholds: logit-spaced around the criterion, scaled by exp(spacing)
vector sdt_thresholds_parsimonious_rating(real criterion, real spacing,
                                          int K_full) {
  int n_thresh = K_full - 1;
  vector[n_thresh] thresholds;
  for (k in 1:n_thresh) {
    real gk = log(k * 1.0 / (K_full - k));
    thresholds[k] = criterion + exp(spacing) * gk;
  }
  return thresholds;
}

// equidistant thresholds: exp(spacing) apart, centred on the criterion
vector sdt_thresholds_equidistant_rating(real criterion, real spacing,
                                         int K_full) {
  int n_thresh = K_full - 1;
  vector[n_thresh] thresholds;
  for (k in 1:n_thresh) {
    thresholds[k] = criterion + (k - K_full / 2.0) * exp(spacing);
  }
  return thresholds;
}

// K - 1 thresholds from K - 2 interval widths, centred on the criterion
vector sdt_place_thresholds_rating(real criterion, vector gaps, int K_full) {
  int n_thresh = K_full - 1;
  vector[n_thresh] thresholds;
  int lo;
  int hi;
  if (K_full % 2 == 0) {
    lo = K_full %/% 2;
    hi = lo;
    thresholds[lo] = criterion;
  } else {
    lo = (K_full - 1) %/% 2;
    hi = lo + 1;
    thresholds[lo] = criterion - gaps[lo] / 2;
    thresholds[hi] = criterion + gaps[lo] / 2;
  }
  if (hi < n_thresh) {
    for (k in (hi + 1):n_thresh) {
      thresholds[k] = thresholds[k - 1] + gaps[k - 1];
    }
  }
  if (lo > 1) {
    for (j in 1:(lo - 1)) {
      int k = lo - j;
      thresholds[k] = thresholds[k + 1] - gaps[k];
    }
  }
  return thresholds;
}

// log_distance thresholds: each interval width is exp(delta)
vector sdt_thresholds_log_distance_rating(real criterion,
                                          array[] real deltas,
                                          int K_full) {
  int n_gaps = K_full - 2;
  vector[n_gaps] gaps;
  for (j in 1:n_gaps) {
    gaps[j] = exp(deltas[j]);
  }
  return sdt_place_thresholds_rating(criterion, gaps, K_full);
}

// log_ratio thresholds (Paulewicz & Blaut, 2022): outer widths are multiples of the central widths
vector sdt_thresholds_log_ratio_rating(real criterion,
                                       array[] real deltas,
                                       int K_full) {
  int n_gaps = K_full - 2;
  vector[n_gaps] gaps;
  for (j in 1:n_gaps) {
    gaps[j] = exp(deltas[j]);
  }
  if (K_full % 2 == 0) {
    int m = K_full %/% 2;
    gaps[m - 1] *= gaps[m];
    if (m + 1 <= n_gaps) {
      for (j in (m + 1):n_gaps) gaps[j] *= gaps[m];
    }
    if (m - 2 >= 1) {
      for (j in 1:(m - 2)) gaps[j] *= gaps[m - 1];
    }
  } else {
    int g = (K_full - 1) %/% 2;
    gaps[g + 1] *= gaps[g];
    gaps[g - 1] *= gaps[g];
    if (g + 2 <= n_gaps) {
      for (j in (g + 2):n_gaps) gaps[j] *= gaps[g + 1];
    }
    if (g - 2 >= 1) {
      for (j in 1:(g - 2)) gaps[j] *= gaps[g - 1];
    }
  }
  return sdt_place_thresholds_rating(criterion, gaps, K_full);
}

// softmax thresholds: interval widths are shares of a total width set by spacing
vector sdt_thresholds_softmax_rating(real criterion, real spacing,
                                     array[] real deltas,
                                     int K_full) {
  int n_intervals = K_full - 2;
  vector[n_intervals] logits;
  if (n_intervals > 1) {
    for (j in 1:(n_intervals - 1)) {
      logits[j] = deltas[j];
    }
  }
  logits[n_intervals] = 0;
  return sdt_place_thresholds_rating(
    criterion, softmax(logits) * n_intervals * exp(spacing), K_full
  );
}

// thresh_type: 1 = parsimonious, 2 = equidistant, 3 = log_distance, 4 = log_ratio, 5 = softmax
vector sdt_make_thresholds_rating(real criterion, real spacing,
                                  array[] real deltas,
                                  int K_full, int thresh_type) {
  if (thresh_type == 1)
    return sdt_thresholds_parsimonious_rating(criterion, spacing, K_full);
  if (thresh_type == 2)
    return sdt_thresholds_equidistant_rating(criterion, spacing, K_full);
  if (thresh_type == 3)
    return sdt_thresholds_log_distance_rating(criterion, deltas, K_full);
  if (thresh_type == 4)
    return sdt_thresholds_log_ratio_rating(criterion, deltas, K_full);
  return sdt_thresholds_softmax_rating(criterion, spacing, deltas, K_full);
}

// log probability of rating category cat; d is d_a and sdratio the log of the SD ratio
real sdt_rating_logmu_cat(int cat, vector thresholds,
                          real d, real sdratio, real stimulus,
                          int dist_type) {
  int K_full = num_elements(thresholds) + 1;
  real sigma = exp(sdratio);
  real shift = d * sdt_rms_scale(sigma) * stimulus;
  real scale = stimulus > 0.5 ? sigma : 1.0;

  if (cat == 1) {
    return sdt_log_cumprob((thresholds[1] - shift) / scale, dist_type);
  }
  if (cat == K_full) {
    return sdt_log_one_minus_cumprob(
      (thresholds[K_full - 1] - shift) / scale, dist_type
    );
  }

  real eta_lo = (thresholds[cat - 1] - shift) / scale;
  real eta_hi = (thresholds[cat] - shift) / scale;
  if (eta_lo > 0) {
    return log_diff_exp(sdt_log_one_minus_cumprob(eta_lo, dist_type),
                        sdt_log_one_minus_cumprob(eta_hi, dist_type));
  }
  return log_diff_exp(sdt_log_cumprob(eta_hi, dist_type),
                      sdt_log_cumprob(eta_lo, dist_type));
}
