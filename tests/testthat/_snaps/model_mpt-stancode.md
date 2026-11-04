# the Stan code of the two-high-threshold model is unchanged

    Code
      golden_stancode(bmf(D ~ 1 + (1 | id), g ~ 1), mpt_2htm_data(n_id = 2), mpt(
        mpt_2htm_trees(), tree_id = "item_type"))
    Output
      functions {
        /* multinomial-logit log-PMF
         * Args:
         *   y: array of integer response values
         *   mu: vector of category logit probabilities
         * Returns:
         *   a scalar to be added to the log posterior
         */
         real multinomial_logit2_lpmf(array[] int y, vector mu) {
           return multinomial_lpmf(y | softmax(mu));
         }
      }
      data {
        int<lower=1> N;  // total number of observations
        int<lower=2> ncat;  // number of categories
        array[N, ncat] int Y;  // response array
        array[N] int trials;  // number of trials
        int<lower=1> K_D;  // number of population-level effects
        matrix[N, K_D] X_D;  // population-level design matrix
        int<lower=1> K_g;  // number of population-level effects
        matrix[N, K_g] X_g;  // population-level design matrix
        // covariates for non-linear functions
        array[N] int C_old_1;
        array[N] int C_old_2;
        // covariates for non-linear functions
        array[N] int C_new_1;
        array[N] int C_new_2;
        // data for group-level effects of ID 1
        int<lower=1> N_1;  // number of grouping levels
        int<lower=1> M_1;  // number of coefficients per level
        array[N] int<lower=1> J_1;  // grouping indicator per observation
        // group-level predictor values
        vector[N] Z_1_D_1;
        int prior_only;  // should the likelihood be ignored?
      }
      transformed data {
      }
      parameters {
        vector[K_D] b_D;  // regression coefficients
        vector[K_g] b_g;  // regression coefficients
        vector<lower=0>[M_1] sd_1;  // group-level standard deviations
        array[M_1] vector[N_1] z_1;  // standardized group-level effects
      }
      transformed parameters {
        vector[N_1] r_1_D_1;  // actual group-level effects
        // prior contributions to the log posterior
        real lprior = 0;
        r_1_D_1 = (sd_1[1] * (z_1[1]));
        lprior += logistic_lpdf(b_D[1] | 0, 1);
        lprior += logistic_lpdf(b_g[1] | 0, 1);
        lprior += exponential_lpdf(sd_1 | 1);
      }
      model {
        // likelihood including constants
        if (!prior_only) {
          // initialize linear predictor term
          vector[N] nlp_D = rep_vector(0.0, N);
          // initialize linear predictor term
          vector[N] nlp_g = rep_vector(0.0, N);
          // initialize non-linear predictor term
          vector[N] nlp_old;
          // initialize non-linear predictor term
          vector[N] nlp_new;
          // initialize non-linear predictor term
          vector[N] muold;
          // initialize non-linear predictor term
          vector[N] munew;
          // linear predictor matrix
          array[N] vector[ncat] mu;
          nlp_D += X_D * b_D;
          nlp_g += X_g * b_g;
          for (n in 1:N) {
            // add more terms to the linear predictor
            nlp_D[n] += r_1_D_1[J_1[n]] * Z_1_D_1[n];
          }
          for (n in 1:N) {
            // compute non-linear predictor values
            nlp_old[n] = (C_old_1[n] * (inv_logit(nlp_D[n]) + (1 - inv_logit(nlp_D[n])) * inv_logit(nlp_g[n])) + C_old_2[n] * ((1 - inv_logit(nlp_D[n])) * inv_logit(nlp_g[n])));
          }
          for (n in 1:N) {
            // compute non-linear predictor values
            nlp_new[n] = (C_new_1[n] * ((1 - inv_logit(nlp_D[n])) * (1 - inv_logit(nlp_g[n]))) + C_new_2[n] * (inv_logit(nlp_D[n]) + (1 - inv_logit(nlp_D[n])) * (1 - inv_logit(nlp_g[n]))));
          }
          for (n in 1:N) {
            // compute non-linear predictor values
            muold[n] = (log(nlp_old[n]));
          }
          for (n in 1:N) {
            // compute non-linear predictor values
            munew[n] = (log(nlp_new[n]));
          }
          for (n in 1:N) {
            mu[n] = transpose([muold[n], munew[n]]);
          }
          for (n in 1:N) {
            target += multinomial_logit2_lpmf(Y[n] | mu[n]);
          }
        }
        // priors including constants
        target += lprior;
        target += std_normal_lpdf(z_1[1]);
      }
      generated quantities {
      }

# the Stan code of the pair-clustering model is unchanged

    Code
      golden_stancode(bmf(cp ~ 1 + (1 | id), rp ~ 1), data.frame(id = factor(1:2), C = 20L,
      E = 12L, U = 8L), mpt(pcm_tree))
    Output
      functions {
        /* multinomial-logit log-PMF
         * Args:
         *   y: array of integer response values
         *   mu: vector of category logit probabilities
         * Returns:
         *   a scalar to be added to the log posterior
         */
         real multinomial_logit2_lpmf(array[] int y, vector mu) {
           return multinomial_lpmf(y | softmax(mu));
         }
      }
      data {
        int<lower=1> N;  // total number of observations
        int<lower=2> ncat;  // number of categories
        array[N, ncat] int Y;  // response array
        array[N] int trials;  // number of trials
        int<lower=1> K_cp;  // number of population-level effects
        matrix[N, K_cp] X_cp;  // population-level design matrix
        int<lower=1> K_rp;  // number of population-level effects
        matrix[N, K_rp] X_rp;  // population-level design matrix
        // data for group-level effects of ID 1
        int<lower=1> N_1;  // number of grouping levels
        int<lower=1> M_1;  // number of coefficients per level
        array[N] int<lower=1> J_1;  // grouping indicator per observation
        // group-level predictor values
        vector[N] Z_1_cp_1;
        int prior_only;  // should the likelihood be ignored?
      }
      transformed data {
      }
      parameters {
        vector[K_cp] b_cp;  // regression coefficients
        vector[K_rp] b_rp;  // regression coefficients
        vector<lower=0>[M_1] sd_1;  // group-level standard deviations
        array[M_1] vector[N_1] z_1;  // standardized group-level effects
      }
      transformed parameters {
        vector[N_1] r_1_cp_1;  // actual group-level effects
        // prior contributions to the log posterior
        real lprior = 0;
        r_1_cp_1 = (sd_1[1] * (z_1[1]));
        lprior += logistic_lpdf(b_cp[1] | 0, 1);
        lprior += logistic_lpdf(b_rp[1] | 0, 1);
        lprior += exponential_lpdf(sd_1 | 1);
      }
      model {
        // likelihood including constants
        if (!prior_only) {
          // initialize linear predictor term
          vector[N] nlp_cp = rep_vector(0.0, N);
          // initialize linear predictor term
          vector[N] nlp_rp = rep_vector(0.0, N);
          // initialize non-linear predictor term
          vector[N] nlp_C;
          // initialize non-linear predictor term
          vector[N] nlp_E;
          // initialize non-linear predictor term
          vector[N] nlp_U;
          // initialize non-linear predictor term
          vector[N] muC;
          // initialize non-linear predictor term
          vector[N] muE;
          // initialize non-linear predictor term
          vector[N] muU;
          // linear predictor matrix
          array[N] vector[ncat] mu;
          nlp_cp += X_cp * b_cp;
          nlp_rp += X_rp * b_rp;
          for (n in 1:N) {
            // add more terms to the linear predictor
            nlp_cp[n] += r_1_cp_1[J_1[n]] * Z_1_cp_1[n];
          }
          for (n in 1:N) {
            // compute non-linear predictor values
            nlp_C[n] = ((inv_logit(nlp_cp[n]) + (1 - inv_logit(nlp_cp[n])) * inv_logit(nlp_rp[n]) * inv_logit(nlp_rp[n])));
          }
          for (n in 1:N) {
            // compute non-linear predictor values
            nlp_E[n] = ((2 * (1 - inv_logit(nlp_cp[n])) * inv_logit(nlp_rp[n]) * (1 - inv_logit(nlp_rp[n]))));
          }
          for (n in 1:N) {
            // compute non-linear predictor values
            nlp_U[n] = (((1 - inv_logit(nlp_cp[n])) * (1 - inv_logit(nlp_rp[n])) * (1 - inv_logit(nlp_rp[n]))));
          }
          for (n in 1:N) {
            // compute non-linear predictor values
            muC[n] = (log(nlp_C[n]));
          }
          for (n in 1:N) {
            // compute non-linear predictor values
            muE[n] = (log(nlp_E[n]));
          }
          for (n in 1:N) {
            // compute non-linear predictor values
            muU[n] = (log(nlp_U[n]));
          }
          for (n in 1:N) {
            mu[n] = transpose([muC[n], muE[n], muU[n]]);
          }
          for (n in 1:N) {
            target += multinomial_logit2_lpmf(Y[n] | mu[n]);
          }
        }
        // priors including constants
        target += lprior;
        target += std_normal_lpdf(z_1[1]);
      }
      generated quantities {
      }

# the Stan code of the Oberauer and Lewandowsky (2019) model is unchanged

    Code
      golden_stancode(bmf(Pi ~ 0 + cond + (0 + cond || ID), Pb ~ 0 + cond + (0 +
      cond || ID), Pd ~ 1 + (1 | ID)), dat, mpt(trees, tree_id = "tree"))
    Output
      functions {
        /* multinomial-logit log-PMF
         * Args:
         *   y: array of integer response values
         *   mu: vector of category logit probabilities
         * Returns:
         *   a scalar to be added to the log posterior
         */
         real multinomial_logit2_lpmf(array[] int y, vector mu) {
           return multinomial_lpmf(y | softmax(mu));
         }
      }
      data {
        int<lower=1> N;  // total number of observations
        int<lower=2> ncat;  // number of categories
        array[N, ncat] int Y;  // response array
        array[N] int trials;  // number of trials
        int<lower=1> K_Pi;  // number of population-level effects
        matrix[N, K_Pi] X_Pi;  // population-level design matrix
        int<lower=1> K_Pb;  // number of population-level effects
        matrix[N, K_Pb] X_Pb;  // population-level design matrix
        int<lower=1> K_Pd;  // number of population-level effects
        matrix[N, K_Pd] X_Pd;  // population-level design matrix
        // covariates for non-linear functions
        array[N] int C_corr_1;
        array[N] int C_corr_2;
        // covariates for non-linear functions
        array[N] int C_other_1;
        array[N] int C_other_2;
        // covariates for non-linear functions
        array[N] int C_dist_1;
        array[N] int C_dist_2;
        // covariates for non-linear functions
        array[N] int C_npl_1;
        array[N] int C_npl_2;
        // covariates for non-linear functions
        array[N] int C_mudist_1;
        // data for group-level effects of ID 1
        int<lower=1> N_1;  // number of grouping levels
        int<lower=1> M_1;  // number of coefficients per level
        array[N] int<lower=1> J_1;  // grouping indicator per observation
        // group-level predictor values
        vector[N] Z_1_Pi_1;
        vector[N] Z_1_Pi_2;
        vector[N] Z_1_Pi_3;
        // data for group-level effects of ID 2
        int<lower=1> N_2;  // number of grouping levels
        int<lower=1> M_2;  // number of coefficients per level
        array[N] int<lower=1> J_2;  // grouping indicator per observation
        // group-level predictor values
        vector[N] Z_2_Pb_1;
        vector[N] Z_2_Pb_2;
        vector[N] Z_2_Pb_3;
        // data for group-level effects of ID 3
        int<lower=1> N_3;  // number of grouping levels
        int<lower=1> M_3;  // number of coefficients per level
        array[N] int<lower=1> J_3;  // grouping indicator per observation
        // group-level predictor values
        vector[N] Z_3_Pd_1;
        int prior_only;  // should the likelihood be ignored?
      }
      transformed data {
      }
      parameters {
        vector[K_Pi] b_Pi;  // regression coefficients
        vector[K_Pb] b_Pb;  // regression coefficients
        vector[K_Pd] b_Pd;  // regression coefficients
        vector<lower=0>[M_1] sd_1;  // group-level standard deviations
        array[M_1] vector[N_1] z_1;  // standardized group-level effects
        vector<lower=0>[M_2] sd_2;  // group-level standard deviations
        array[M_2] vector[N_2] z_2;  // standardized group-level effects
        vector<lower=0>[M_3] sd_3;  // group-level standard deviations
        array[M_3] vector[N_3] z_3;  // standardized group-level effects
      }
      transformed parameters {
        vector[N_1] r_1_Pi_1;  // actual group-level effects
        vector[N_1] r_1_Pi_2;  // actual group-level effects
        vector[N_1] r_1_Pi_3;  // actual group-level effects
        vector[N_2] r_2_Pb_1;  // actual group-level effects
        vector[N_2] r_2_Pb_2;  // actual group-level effects
        vector[N_2] r_2_Pb_3;  // actual group-level effects
        vector[N_3] r_3_Pd_1;  // actual group-level effects
        // prior contributions to the log posterior
        real lprior = 0;
        r_1_Pi_1 = (sd_1[1] * (z_1[1]));
        r_1_Pi_2 = (sd_1[2] * (z_1[2]));
        r_1_Pi_3 = (sd_1[3] * (z_1[3]));
        r_2_Pb_1 = (sd_2[1] * (z_2[1]));
        r_2_Pb_2 = (sd_2[2] * (z_2[2]));
        r_2_Pb_3 = (sd_2[3] * (z_2[3]));
        r_3_Pd_1 = (sd_3[1] * (z_3[1]));
        lprior += logistic_lpdf(b_Pi | 0, 1);
        lprior += logistic_lpdf(b_Pb | 0, 1);
        lprior += logistic_lpdf(b_Pd[1] | 0, 1);
        lprior += exponential_lpdf(sd_1 | 1);
        lprior += exponential_lpdf(sd_2 | 1);
        lprior += exponential_lpdf(sd_3 | 1);
      }
      model {
        // likelihood including constants
        if (!prior_only) {
          // initialize linear predictor term
          vector[N] nlp_Pi = rep_vector(0.0, N);
          // initialize linear predictor term
          vector[N] nlp_Pb = rep_vector(0.0, N);
          // initialize linear predictor term
          vector[N] nlp_Pd = rep_vector(0.0, N);
          // initialize non-linear predictor term
          vector[N] nlp_corr;
          // initialize non-linear predictor term
          vector[N] nlp_other;
          // initialize non-linear predictor term
          vector[N] nlp_dist;
          // initialize non-linear predictor term
          vector[N] nlp_npl;
          // initialize non-linear predictor term
          vector[N] mucorr;
          // initialize non-linear predictor term
          vector[N] muother;
          // initialize non-linear predictor term
          vector[N] mudist;
          // initialize non-linear predictor term
          vector[N] munpl;
          // linear predictor matrix
          array[N] vector[ncat] mu;
          nlp_Pi += X_Pi * b_Pi;
          nlp_Pb += X_Pb * b_Pb;
          nlp_Pd += X_Pd * b_Pd;
          for (n in 1:N) {
            // add more terms to the linear predictor
            nlp_Pi[n] += r_1_Pi_1[J_1[n]] * Z_1_Pi_1[n] + r_1_Pi_2[J_1[n]] * Z_1_Pi_2[n] + r_1_Pi_3[J_1[n]] * Z_1_Pi_3[n];
          }
          for (n in 1:N) {
            // add more terms to the linear predictor
            nlp_Pb[n] += r_2_Pb_1[J_2[n]] * Z_2_Pb_1[n] + r_2_Pb_2[J_2[n]] * Z_2_Pb_2[n] + r_2_Pb_3[J_2[n]] * Z_2_Pb_3[n];
          }
          for (n in 1:N) {
            // add more terms to the linear predictor
            nlp_Pd[n] += r_3_Pd_1[J_3[n]] * Z_3_Pd_1[n];
          }
          for (n in 1:N) {
            // compute non-linear predictor values
            nlp_corr[n] = (C_corr_1[n] * (inv_logit(nlp_Pi[n]) * inv_logit(nlp_Pb[n]) + inv_logit(nlp_Pi[n]) * (1 - inv_logit(nlp_Pb[n])) * (1 - inv_logit(nlp_Pd[n])) * 0.2 + (1 - inv_logit(nlp_Pi[n])) * 0.0666666666666667) + C_corr_2[n] * (inv_logit(nlp_Pi[n]) * inv_logit(nlp_Pb[n]) + inv_logit(nlp_Pi[n]) * (1 - inv_logit(nlp_Pb[n])) * 0.2 + (1 - inv_logit(nlp_Pi[n])) * 0.0666666666666667));
          }
          for (n in 1:N) {
            // compute non-linear predictor values
            nlp_other[n] = (C_other_1[n] * (inv_logit(nlp_Pi[n]) * (1 - inv_logit(nlp_Pb[n])) * (1 - inv_logit(nlp_Pd[n])) * 0.8 + (1 - inv_logit(nlp_Pi[n])) * 0.266666666666667) + C_other_2[n] * (inv_logit(nlp_Pi[n]) * (1 - inv_logit(nlp_Pb[n])) * 0.8 + (1 - inv_logit(nlp_Pi[n])) * 0.266666666666667));
          }
          for (n in 1:N) {
            // compute non-linear predictor values
            nlp_dist[n] = (C_dist_1[n] * (inv_logit(nlp_Pi[n]) * (1 - inv_logit(nlp_Pb[n])) * inv_logit(nlp_Pd[n]) + (1 - inv_logit(nlp_Pi[n])) * 0.333333333333333) + C_dist_2[n] * (1));
          }
          for (n in 1:N) {
            // compute non-linear predictor values
            nlp_npl[n] = (C_npl_1[n] * ((1 - inv_logit(nlp_Pi[n])) * 0.333333333333333) + C_npl_2[n] * ((1 - inv_logit(nlp_Pi[n])) * 0.666666666666667));
          }
          for (n in 1:N) {
            // compute non-linear predictor values
            mucorr[n] = (log(nlp_corr[n]));
          }
          for (n in 1:N) {
            // compute non-linear predictor values
            muother[n] = (log(nlp_other[n]));
          }
          for (n in 1:N) {
            // compute non-linear predictor values
            mudist[n] = (C_mudist_1[n] * log(nlp_dist[n]) + (1 - C_mudist_1[n]) * ( - 100));
          }
          for (n in 1:N) {
            // compute non-linear predictor values
            munpl[n] = (log(nlp_npl[n]));
          }
          for (n in 1:N) {
            mu[n] = transpose([mucorr[n], muother[n], mudist[n], munpl[n]]);
          }
          for (n in 1:N) {
            target += multinomial_logit2_lpmf(Y[n] | mu[n]);
          }
        }
        // priors including constants
        target += lprior;
        target += std_normal_lpdf(z_1[1]);
        target += std_normal_lpdf(z_1[2]);
        target += std_normal_lpdf(z_1[3]);
        target += std_normal_lpdf(z_2[1]);
        target += std_normal_lpdf(z_2[2]);
        target += std_normal_lpdf(z_2[3]);
        target += std_normal_lpdf(z_3[1]);
      }
      generated quantities {
      }

