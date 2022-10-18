#ifndef GIBBS_INVQ_H
#define GIBBS_INVQ_H

void gibbs_InvQ(struct str_state* last,  // OUT+IN last known values of generated parameters
                struct str_param* param, // IN hyperparameters
                int* spec,            // [5]                class-specific parameters
                // order:   [0] tau_num, 
                //          [1] c_ord,
                //          [2] InvSigma, 
                //          [3] InvQ, 
                //          [4] naY
                int* whatsave,        // [33]:              what parameters are saved
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
                int* dimswithG,       // [33]:              what parameters are saved
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
                int* G,                 // IN [1] number of classes
                int* totnran,           // IN [1] total number of random regressors
                int* m,                 // IN [1] iteration number
                double* Q,              // [(B+M)*totnran*(totnran+1)/2(* G)] Prior scale matrix of InvSigma
                double* detInvQ         // [(B+M)(* G)] determinant of InvQ
); 

#endif