// Shared confidence-rating SDT helper functions
// Used by rating SDT, DPSDT, and meta-d' custom families.

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
  int mid = (K_full - 1) / 2 + 1;
  vector[n_thresh] thresholds;
  for (k in 1:n_thresh) {
    thresholds[k] = criterion + (k - mid) * exp(spacing);
  }
  return thresholds;
}

vector sdt_thresholds_log_distance_rating(real criterion,
                                          array[] real deltas,
                                          int K_full) {
  int n_thresh = K_full - 1;
  int mid = (K_full - 1) / 2 + 1;
  vector[n_thresh] thresholds;
  thresholds[mid] = criterion;

  if (mid < n_thresh) {
    for (k in (mid + 1):n_thresh) {
      int delta_idx = k <= mid ? k : k - 1;
      thresholds[k] = thresholds[k - 1] + exp(deltas[delta_idx]);
    }
  }

  if (mid > 1) {
    for (k in (mid - 1):1) {
      thresholds[k] = thresholds[k + 1] - exp(deltas[k]);
    }
  }

  return thresholds;
}

vector sdt_thresholds_log_ratio_rating(real criterion,
                                       array[] real deltas,
                                       int K_full) {
  int n_thresh = K_full - 1;
  int mid = (K_full - 1) / 2 + 1;
  vector[n_thresh] thresholds;
  thresholds[mid] = criterion;

  if (mid < n_thresh) {
    int spread_idx = mid;
    real spread_above = exp(deltas[spread_idx]);
    thresholds[mid + 1] = criterion + spread_above;

    if (mid > 1) {
      real spread_below = exp(deltas[mid - 1]) * spread_above;
      thresholds[mid - 1] = criterion - spread_below;

      if (mid - 2 >= 1) {
        for (k in (mid - 2):1) {
          thresholds[k] = thresholds[k + 1] -
            exp(deltas[k]) * spread_below;
        }
      }
    }

    if (mid + 2 <= n_thresh) {
      for (k in (mid + 2):n_thresh) {
        int delta_idx = k - 1;
        thresholds[k] = thresholds[k - 1] +
          exp(deltas[delta_idx]) * spread_above;
      }
    }
  }

  return thresholds;
}

vector sdt_thresholds_softmax_rating(real criterion, real spacing,
                                     array[] real deltas,
                                     int K_full) {
  int n_thresh = K_full - 1;
  int mid = (K_full - 1) / 2 + 1;
  int n_intervals = K_full - 2;
  vector[n_thresh] thresholds;

  thresholds[mid] = criterion;

  if (n_intervals == 1) {
    real interval = exp(spacing);
    thresholds[mid - 1] = criterion - interval;
    thresholds[mid + 1] = criterion + interval;
    return thresholds;
  }

  vector[n_intervals] logits;
  vector[n_intervals] intervals;
  if (n_intervals > 1) {
    for (j in 1:(n_intervals - 1)) {
      logits[j] = deltas[j];
    }
  }
  logits[n_intervals] = 0;
  intervals = softmax(logits) * n_intervals * exp(spacing);

  if (mid < n_thresh) {
    thresholds[mid + 1] = criterion + intervals[mid];
    if (mid + 2 <= n_thresh) {
      for (k in (mid + 2):n_thresh) {
        thresholds[k] = thresholds[k - 1] + intervals[k - 1];
      }
    }
  }

  if (mid > 1) {
    thresholds[mid - 1] = criterion - intervals[mid - 1];
    if (mid - 2 >= 1) {
      for (k in (mid - 2):1) {
        thresholds[k] = thresholds[k + 1] - intervals[k];
      }
    }
  }

  return thresholds;
}

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

real sdt_rating_logprob_category(int category, int stimulus,
                                 vector thresholds,
                                 real dprime, real criterion,
                                 real sdratio,
                                 int dist_type) {
  int K_full = num_elements(thresholds) + 1;
  real shift = dprime / 2.0 * (2 * stimulus - 1);
  real scale = stimulus == 1 ? exp(sdratio) : 1.0;

  if (category == 1) {
    return sdt_log_cumprob((thresholds[1] - shift) / scale, dist_type);
  }
  if (category == K_full) {
    return sdt_log_one_minus_cumprob(
      (thresholds[K_full - 1] - shift) / scale, dist_type
    );
  }

  return log_diff_exp(
    sdt_log_cumprob((thresholds[category] - shift) / scale, dist_type),
    sdt_log_cumprob((thresholds[category - 1] - shift) / scale, dist_type)
  );
}

real sdt_rating_row_lpmf(int y, int category, int stimulus,
                         int n_ratings, int dist_type, int thresh_type,
                         real dprime, real criterion, real spacing,
                         real sdratio,
                         array[] real deltas) {
  if (y == 0) return 0;

  vector[n_ratings - 1] thresholds = sdt_make_thresholds_rating(
    criterion, spacing, deltas, n_ratings, thresh_type
  );
  return y * sdt_rating_logprob_category(
    category, stimulus, thresholds, dprime, criterion, sdratio, dist_type
  );
}

