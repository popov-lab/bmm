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
    real s_sq = square(s);
    real k = drift / s_sq;
    real w = square(k);
    real x0 = square(bound) * w;
    int series = x0 < 0.49;
    real h1;
    real h2;
    real h3;
    real h4;
    if (series) {
      real scale = 2 * square(bound) / s_sq;
      real scale_sq = square(scale);
      h1 = scale * ezdm_log_sinhc_d1(x0);
      h2 = scale_sq * ezdm_log_sinhc_d2(x0);
      h3 = scale_sq * scale * ezdm_log_sinhc_d3(x0);
      h4 = square(scale_sq) * ezdm_log_sinhc_d4(x0);
    } else {
      real t0 = sqrt(x0);
      real e0 = exp(-2 * t0);
      real p0 = ezdm_coth_term(t0, e0);
      real q0 = ezdm_csch_term(t0, e0);
      h1 = p0;
      h2 = -p0 - q0;
      h3 = ezdm_cgf_d3(p0, q0);
      h4 = ezdm_cgf_d4(p0, q0, t0);
    }

    return binomial_lpmf(hits | trials, ezdm_pc(0.5 * bound, bound, k))
           + ezdm_boundary_lpdf(mrt | vrt, trials, ndt, 0.5 * bound, w, s_sq,
                                series, h1, h2, h3, h4);
  }
