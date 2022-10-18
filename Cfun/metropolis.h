#ifndef METROPOLIS_H
#define METROPOLIS_H

void metropolis(struct str_state* last,   // IN last known values of generated parameters
                double* hyperparameter,   // IN hyperparameters
                double* Y,                // IN [*N]
                double* X,                // IN [N*(#regr)]    regressors
                int* spec,            // [5]                class-specific parameters
                // order:   [0] tau_num, 
                //          [1] c_ord,
                //          [2] InvSigma, 
                //          [3] InvQ, 
                //          [4] naY
                int* dims,            // [33]:              what parameters are saved
                // order:   [0-1]   beta_num_fix, beta_num,
                //          [2-4]   tau_num, sd_num, var_num,
                //          [5-6]   beta_poi_fix, beta_poi,
                //          [7-8]   beta_bin_fix, beta_bin,
                //          [9-10]  beta_ord_fix, beta_ord,
                //          [11-13] c_ord, a_ord, pi_ord,
                //          [14-15] beta_cat_fix, beta_cat,
                //          [16-20] InvSigma, Sigma, sdSigma, corSigma, detInvSigma,
                //          [21-23] InvQ, Q, detInvQ,
                //          [24]    b,
                //          [25-26] w, ng,
                //          [27-29] loglik, pUig, U,
                //          [30,31] Gplus, e0,
                //          [32]    naY
                double* predictor_num,
                double* predictor_poi,
                double* predictor_bin,
                double* predictor_ord,
                double* predictor_cat,
                int* pred_skip,         // IN [1] Which part of predictor is to be skipped
                int* N,                 // IN [1] number of observations
                int* n,                 // IN [1] number of subjects
                int* g,                 // IN [1] *G will be passed for FIXED, g for group-specific data
                int* Kord,              // IN [1] total number of categories of Ord outcome -1 (useless)
                int* Kcat,              // IN [1] total number of categories of Cat outcome -1 (useless)
                int* is_cat,            // IN [1] Is it categorical variable?
                int* kspec_bi_cat,      // IN [1] TRUE = each cat var has different set of bi for each of K-1 levels
                //        FALSE = all levels of cat var have one set of bi - compares last category with others
                int* nfix,              // IN [1] number of FIXED  regressors or each response
                int* cumnfix,           // IN [1] where to start in FormulaF
                int* totnfix,           // IN [1] total dimension of fixed parameter (useless)
                int* FormulaF,          // IN [nfix]  numbers of columns of X that should be used for FIXED  effects of modelled responses
                int* n_i,               // IN [n] number of observations dedicated to j-th subject
                int* i_j,               // IN [n * max_n_j] indeces dedicated to j-th subject 
                int* nUg,               // IN [G] number of subjects in classes 1, ..., G
                int* listUi,            // IN [n*G] row by row indices of subjects in group g
                double* theta,          // IN+OUT [dims[.]] theta from which we start the random walk
                //double* loglik,         // IN+OUT [1] the value of maximized log-likelihood
                double* proposal_chol,  // IN [dims[.] * (dims[.]+1) / 2]
                // Functions used
                void (*fun)(struct str_state*, double*, double*, double*, int*, int*, 
                      double*, double*, double*, double*, double*, int*,
                      int*, int*, int*, int*, int*, int*, 
                      int*, int*, int*, int*, int*, int*, int*, int*,
                      double*, double*),
                      // Tuning parameters
                      struct str_tuning* tuning, // tuning parameters
                      double* const_proposal   // IN [1] pointer to one of tuning parameters
                  
);

#endif