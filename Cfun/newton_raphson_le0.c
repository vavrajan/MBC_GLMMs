/*
 * Newton-Raphson method for maximization of function dependant on vector paramter theta=b_i
 * using first and second order derivatives
 * Only for random effects that are a bit specific
 */
#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "structures.h"
#include "pdfs_derivatives_le0.h"

void newton_raphson_le0(struct str_state* last,    // OUT+IN last known values of generated parameters
                       struct str_param* param,   // IN hyperparameters
                       int* G,                    // IN [1] number of all clusters
                       int* Gplus,                // IN [1] number of non-empty clusters sum(nUg>0)
                       int* n,                   // IN [1] total number of subjects
                       int* nUg,                  // IN [G] number of subjects currently within g-th cluster 
                       // Proposal distribution parameters
                       double* proposal_prec,      // IN+OUT[1] precision of proposal distribution for log(e0)
                       // Tuning parameters
                       struct str_tuning* tuning,
                       int* iter,                   // OUT [1] number of iterations performed before convergence
                       int* converged,              // OUT [1] indicator of convergence
                       double* max_value            // OUT [1] the value of maximized log-likelihood
){
  
  // Auxiliary parameters
  int j;
  double loglik;
  
  double le0 = log(*((*last).e0));
  double grad = 0.0;
  double shift = 0.0;
  double norm_shift = 0.0;
  
  logple0(last, param, 
          G, Gplus, n, nUg,
          &le0, max_value);
  
  //printf("\nnewton_raphson_le0: max_value = %f", *max_value);
  //fflush(stdout);
  *converged = 0;
  
  for(j = 0; ((j < *((*tuning).maxiter)) & (!(*converged))); j++){
    //printf("newton_raphson_le0: j = %d \n", j);
    //fflush(stdout);
    
    // First compute the first and second derivative
    d_d2_logple0(last, param, 
                 G, Gplus, n, nUg,
                 &le0, &grad, proposal_prec);
    //printf("\nnewton_raphson_le0: first = %f, second = %f", grad, *proposal_prec);
    //fflush(stdout);
    
    // Compute the shift
    shift = grad/(*proposal_prec);
    
    // Did we converge?
    norm_shift = ((shift > 0) ? shift : -shift);
    *converged = ( norm_shift < *((*tuning).tolerance));
    //printf("\nnewton_raphson_e0: normshift = %f", norm_shift);
    //fflush(stdout);
    
    // Update new value of theta by adding shift to the current one
    le0 += shift;
    //printf("\nnewton_raphson_e0: le0 = %f", le0);
    //fflush(stdout);
    
    // How about new loglik value
    logple0(last, param, 
            G, Gplus, n, nUg,
            &le0, &loglik);
    //printf("\nnewton_raphson_e0: prev_loglik = %f, loglik = %f, dif = %f", *max_value, loglik, loglik-*max_value);
    //fflush(stdout);
    *max_value = loglik;
    
  }
  
  // Cycle has ended by either convergence at j-th step or at maxiter;
  *iter = j;
  // Update proposal_prec to be at the last value of le0
  d_d2_logple0(last, param, 
               G, Gplus, n, nUg,
               &le0, &grad, proposal_prec);
  
}