  // fail fast if the cell tables computed in R do not match the data Stan received
  {
    int sdm_cell_total = 0;
    for (g in 1:G_sdm_cells) {
      if (sdm_cell_first[g] < 1 || sdm_cell_first[g] > N
          || (g > 1 && sdm_cell_first[g] <= sdm_cell_first[g - 1])
          || sdm_cell_set[g] < 0 || sdm_cell_set[g] > S_sdm_sets) {
        reject("bmm error: cell ", g, " of the SDM option tables does not match the model ",
               "data. Please report this at https://github.com/popov-lab/bmm/issues");
      }
      sdm_cell_total += sdm_cell_count[g];
    }
    if (sdm_cell_total != N) {
      reject("bmm error: sum(sdm_cell_count) = ", sdm_cell_total, " but N = ", N,
             ". Please report this at https://github.com/popov-lab/bmm/issues");
    }
    for (s in 1:S_sdm_sets) {
      if (sdm_set_size[s] < 1 || sdm_set_size[s] > M_sdm_sets) {
        reject("bmm error: option set ", s, " has ", sdm_set_size[s], " options. Please ",
               "report this at https://github.com/popov-lab/bmm/issues");
      }
      for (j in 1:sdm_set_size[s]) {
        if (sdm_set_off[s, j] <= -sdm_K %/% 2 || sdm_set_off[s, j] > sdm_K %/% 2) {
          reject("bmm error: option set ", s, " holds an offset outside the wheel. Please ",
                 "report this at https://github.com/popov-lab/bmm/issues");
        }
      }
    }
    for (n in 1:N) {
      if (Y[n] <= -sdm_K %/% 2 || Y[n] > sdm_K %/% 2) {
        reject("bmm error: response ", n, " is ", Y[n], " steps from the target, outside ",
               "the wheel of ", sdm_K, " options. Please report this at ",
               "https://github.com/popov-lab/bmm/issues");
      }
    }
  }
