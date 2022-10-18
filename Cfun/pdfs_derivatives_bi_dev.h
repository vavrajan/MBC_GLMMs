#ifndef PDFS_DERIVATIVES_BI_DEV_H
#define PDFS_DERIVATIVES_BI_DEV_H

void logpYi_dev(double* Y,                // IN [totnY*N]
                int* isYna,               // [N*sum(nY)]  0 = Y value present, 1 = Y value is NA
                double* X,                // IN [N*(#regr)]    regressors
                double* pbeta_num_fix,     // pointer to beta_num_fix
                double* pbeta_num,         // pointer to beta_num
                double* ptau_num,          // pointer to tau_num
                double* pbeta_poi_fix,     // pointer to beta_poi_fix
                double* pbeta_poi,         // pointer to beta_poi
                double* pbeta_bin_fix,     // pointer to beta_bin_fix
                double* pbeta_bin,         // pointer to beta_bin
                double* pbeta_ord_fix,     // pointer to beta_ord_fix
                double* pbeta_ord,         // pointer to beta_ord
                double* pc_ord,            // pointer to c_ord
                double* pbeta_cat_fix,     // pointer to beta_cat_fix
                double* pbeta_cat,         // pointer to beta_cat
                double* pInvSigma,         // pointer to InvSigma
                int* N,                   // IN [1] number of observations
                int* n,                   // IN [1] number of subjects
                int* nY,                  // IN [4]
                int* i,                   // IN [1] number of subject 
                int* Kord,                // IN [nY[2]] total number of categories of Ord outcomes -1
                int* Kcat,                // IN [nY[3]] total number of categories of Cat outcomes -1
                int* nfix,                 // IN [totnY] number of FIXED  regressors or each response
                int* cumnfix,              // IN [totnY] where to start in FormulaF
                int* FormulaF,             // IN [totnfix]  numbers of columns of X that should be used for FIXED  effects of modelled responses
                int* ngrp,                 // IN [totnY] number of GROUP-SPECIFIC  regressors or each response
                int* cumngrp,              // IN [totnY] where to start in FormulaG
                int* FormulaG,             // IN [totnfix]  numbers of columns of X that should be used for GROUP-SPECIFIC  effects of modelled responses
                int* nran,                 // IN [totnY] number of RANDOM  regressors or each response
                int* cumnran,              // IN [totnY] where to start in FormulaR
                int* totnran,              // IN [1] total dimension of b_i
                int* FormulaR,             // IN [totnran]  numbers of columns of X that should be used for RANDOM  effects of modelled responses
                int* noff,                 // IN [totnY] number of OFFSET  regressors or each response
                int* cumnoff,              // IN [totnY] where to start in FormulaO
                int* FormulaO,             // IN [totnran]  numbers of columns of X that should be used for OFFSET  effects of modelled responses
                int* ij,                  // IN [1] number of first ij observations of j-th subject to be used
                int* i_j,                 // IN [n * max_n_j] indeces dedicated to j-th subject 
                int* kspec_bi_cat,        // IN [1] TRUE = each cat var has different set of bi for each of K-1 levels
                //        FALSE = all levels of cat var have one set of bi - compares last category with others
                int* print,
                double* theta,            // IN [dims[.]] the parameter in which the function is evaluated
                double* value             // OUT [1]
);

void logpbi_dev(double* Y,              // IN [totnY*N]
                int* isYna,             // [N*sum(nY)]  0 = Y value present, 1 = Y value is NA
                double* X,              // IN [N*(#regr)]    regressors
                double* pbeta_num_fix,     // pointer to beta_num_fix
                double* pbeta_num,         // pointer to beta_num
                double* ptau_num,          // pointer to tau_num
                double* pbeta_poi_fix,     // pointer to beta_poi_fix
                double* pbeta_poi,         // pointer to beta_poi
                double* pbeta_bin_fix,     // pointer to beta_bin_fix
                double* pbeta_bin,         // pointer to beta_bin
                double* pbeta_ord_fix,     // pointer to beta_ord_fix
                double* pbeta_ord,         // pointer to beta_ord
                double* pc_ord,            // pointer to c_ord
                double* pbeta_cat_fix,     // pointer to beta_cat_fix
                double* pbeta_cat,         // pointer to beta_cat
                double* pInvSigma,         // pointer to InvSigma
                int* N,                 // IN [1] number of observations
                int* n,                 // IN [1] number of subjects
                int* nY,                // IN [4]
                int* i,                 // IN [1] number of subject 
                int* Kord,              // IN [nY[2]] total number of categories of Ord outcomes -1
                int* Kcat,              // IN [nY[3]] total number of categories of Cat outcomes -1
                int* nfix,                 // IN [totnY] number of FIXED  regressors or each response
                int* cumnfix,              // IN [totnY] where to start in FormulaF
                int* FormulaF,             // IN [totnfix]  numbers of columns of X that should be used for FIXED  effects of modelled responses
                int* ngrp,                 // IN [totnY] number of GROUP-SPECIFIC  regressors or each response
                int* cumngrp,              // IN [totnY] where to start in FormulaG
                int* FormulaG,             // IN [totnfix]  numbers of columns of X that should be used for GROUP-SPECIFIC  effects of modelled responses
                int* nran,                 // IN [totnY] number of RANDOM  regressors or each response
                int* cumnran,              // IN [totnY] where to start in FormulaR
                int* totnran,              // IN [1] total dimension of b_i
                int* FormulaR,             // IN [totnran]  numbers of columns of X that should be used for RANDOM  effects of modelled responses
                int* noff,                 // IN [totnY] number of OFFSET  regressors or each response
                int* cumnoff,              // IN [totnY] where to start in FormulaO
                int* FormulaO,             // IN [totnran]  numbers of columns of X that should be used for OFFSET  effects of modelled responses
                int* ij,                // IN [1] number of first ij observations of j-th subject to be used
                int* i_j,               // IN [n * max_n_j] indeces dedicated to j-th subject 
                int* kspec_bi_cat,      // IN [1] TRUE = each cat var has different set of bi for each of K-1 levels
                //        FALSE = all levels of cat var have one set of bi - compares last category with others
                int* print,
                double* theta,          // IN [dims[.]] the parameter in which the function is evaluated
                double* value           // OUT [1]
);

void d_d2_logpbi_dev(double* Y,             // IN [totnY*N]
                     int* isYna,            // [N*sum(nY)]  0 = Y value present, 1 = Y value is NA
                     double* X,             // IN [N*(#regr)]    regressors
                     double* pbeta_num_fix,     // pointer to beta_num_fix
                     double* pbeta_num,         // pointer to beta_num
                     double* ptau_num,          // pointer to tau_num
                     double* pbeta_poi_fix,     // pointer to beta_poi_fix
                     double* pbeta_poi,         // pointer to beta_poi
                     double* pbeta_bin_fix,     // pointer to beta_bin_fix
                     double* pbeta_bin,         // pointer to beta_bin
                     double* pbeta_ord_fix,     // pointer to beta_ord_fix
                     double* pbeta_ord,         // pointer to beta_ord
                     double* pc_ord,            // pointer to c_ord
                     double* pbeta_cat_fix,     // pointer to beta_cat_fix
                     double* pbeta_cat,         // pointer to beta_cat
                     double* pInvSigma,         // pointer to InvSigma
                     int* N,                // IN [1] number of observations
                     int* n,                // IN [1] number of subjects
                     int* nY,               // IN [4]
                     int* i,                // IN [1] number of subject 
                     int* Kord,             // IN [nY[2]] total number of categories of Ord outcomes -1
                     int* Kcat,             // IN [nY[3]] total number of categories of Cat outcomes -1
                     int* nfix,                 // IN [totnY] number of FIXED  regressors or each response
                     int* cumnfix,              // IN [totnY] where to start in FormulaF
                     int* FormulaF,             // IN [totnfix]  numbers of columns of X that should be used for FIXED  effects of modelled responses
                     int* ngrp,                 // IN [totnY] number of GROUP-SPECIFIC  regressors or each response
                     int* cumngrp,              // IN [totnY] where to start in FormulaG
                     int* FormulaG,             // IN [totnfix]  numbers of columns of X that should be used for GROUP-SPECIFIC  effects of modelled responses
                     int* nran,                 // IN [totnY] number of RANDOM  regressors or each response
                     int* cumnran,              // IN [totnY] where to start in FormulaR
                     int* totnran,              // IN [1] total dimension of b_i
                     int* FormulaR,             // IN [totnran]  numbers of columns of X that should be used for RANDOM  effects of modelled responses
                     int* noff,                 // IN [totnY] number of OFFSET  regressors or each response
                     int* cumnoff,              // IN [totnY] where to start in FormulaO
                     int* FormulaO,             // IN [totnran]  numbers of columns of X that should be used for OFFSET  effects of modelled responses
                     int* ij,               // IN [1] number of first ij observations of j-th subject to be used
                     int* i_j,              // IN [n * max_n_j] indeces dedicated to j-th subject 
                     int* kspec_bi_cat,     // IN [1] TRUE = each cat var has different set of bi for each of K-1 levels
                     //        FALSE = all levels of cat var have one set of bi - compares last category with others
                     double* theta,         // IN [dims[.]] the parameter in which the function is evaluated
                     double* grad,          // OUT [totnran]
                     double* hess           // OUT [totnran*(totnran+1)/2]
);
#endif