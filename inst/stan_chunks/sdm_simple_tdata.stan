  // precompute chebyshev points
  matrix[200,200] COSN;
  for (m in 1:200) {
    for (i in 1:m) {
      COSN[i,m] = cos((2*i-1)*pi()/(2*m))-1;
    }
  }
  // fail fast if the precomputed run metadata does not describe the data Stan
  // received (it is computed in R and can drift out of sync with the data)
  {
    int sdm_run_total = 0;
    for (g in 1:G_sdm_runs) {
      if (sdm_run_start[g] < 1 || sdm_run_start[g] > N) {
        reject("bmm error: sdm_run_start[", g, "] = ", sdm_run_start[g],
               " is outside the data range [1, ", N, "]. The SDM run metadata ",
               "does not match the model data. Please report this at ",
               "https://github.com/popov-lab/bmm/issues");
      }
      if (g > 1 && sdm_run_start[g] <= sdm_run_start[g - 1]) {
        reject("bmm error: sdm_run_start must be strictly increasing. The SDM ",
               "run metadata does not match the model data. Please report ",
               "this at https://github.com/popov-lab/bmm/issues");
      }
      sdm_run_total += sdm_run_count[g];
    }
    if (sdm_run_total != N) {
      reject("bmm error: sum(sdm_run_count) = ", sdm_run_total, " but N = ", N,
             ". The SDM run metadata does not match the model data. Please ",
             "report this at https://github.com/popov-lab/bmm/issues");
    }
  }
