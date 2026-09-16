  // EZ-Diffusion Model likelihood for aggregated data, symmetric start point.
  //
  // The response counts follow a binomial; the two RT summaries follow the
  // sampling distribution of the mean and variance of `trials` decision times,
  // matched to the exact first four cumulants (issue #407). Decision times are
  // right-skewed, so the independent normal and scaled chi-square terms this
  // replaced understated the sampling variance of `vrt` several-fold and
  // ignored its correlation with `mrt`.
  //
  // With a symmetric start point the decision time is independent of which
  // boundary is hit, so all `trials` responses inform one set of cumulants.
  //
  // mu is a dummy dpar required by brms and is not used.
  real ezdm_3par_lpdf(real mrt, real mu, real drift, real bound, real ndt,
                      real s, real vrt, int hits, int trials) {
    vector[4] cumulants = ezdm_cumulants(0.5 * bound, bound,
                                         square(drift) / pow(s, 4), s);
    vector[4] rt = ezdm_rt_terms(cumulants, trials);
    real pC = ezdm_pc(0.5 * bound, bound, drift / square(s));

    return binomial_lpmf(hits | trials, pC)
           + gamma_lpdf(vrt | rt[1], rt[2])
           + normal_lpdf(mrt | ndt + cumulants[1] + rt[3] * (vrt - cumulants[2]),
                         rt[4]);
  }
