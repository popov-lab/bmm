// K - 1 ordered thresholds on a latent axis from a criterion, a log spacing and log interval deltas;
// shared by the ordered-threshold models (sdt_rating, likert_cumulative); needs nothing

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

vector sdt_thresholds_equidistant_rating(real criterion, real spacing,
                                         int K_full) {
  int n_thresh = K_full - 1;
  vector[n_thresh] thresholds;
  for (k in 1:n_thresh) {
    thresholds[k] = criterion + (k - K_full / 2.0) * exp(spacing);
  }
  return thresholds;
}

// K - 1 thresholds from K - 2 interval widths, centred on the criterion; R twin: .assemble_thresholds()
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
    // Stan for-loops only count up; walk k = lo-1 down to 1 via ascending j
    for (j in 1:(lo - 1)) {
      int k = lo - j;
      thresholds[k] = thresholds[k + 1] - gaps[k];
    }
  }
  return thresholds;
}

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

// log_ratio thresholds (Paulewicz & Blaut, 2022); R twin: .log_ratio_widths()
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

// thresh_type follows .threshold_types in R/helpers-thresholds.R; R twin: .make_thresholds()
vector make_thresholds(real criterion, real spacing, array[] real deltas,
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
