/*
 * Probability density functions evaluated purely in the last state --> for calculating Puig 
 */


#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "structures.h"

#include "my_math.h"
#include "matrmult.h"
#include "cholesky.h"


void logpYi_g(struct str_state* last,  // IN last known values of generated parameters
              struct str_param* param, // IN hyperparameters
              double* Y,                // IN [totnY*N]
              int* isYna,               // IN [*N * totnY]
              int* isYna_inv,           // IN [*N * totnY]
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
              int* g,                 // IN [1] group=cluster number
              int* N,                 // IN [1] number of observations
              int* n,                 // IN [1] number of subjects
              int* nY,                // IN [4]
              int* i,                 // IN [1] number of subject 
              int* Kord,              // IN [nY[2]] total number of categories of Ord outcomes -1
              int* Kcat,              // IN [nY[3]] total number of categories of Cat outcomes -1
              int* ngrp,              // IN [1] number of GROUP-SPECIFIC  regressors or each response
              int* cumngrp,           // IN [1] where to start in FormulaG
              int* FormulaG,          // IN [totngrp]  numbers of columns of X that should be used for GROUP-SPECIFIC  effects of modelled responses
              int* n_i,               // IN [n] number of observations dedicated to i-th subject
              int* i_j,               // IN [n * max_n_j] indeces dedicated to i-th subject 
              int* print,             // IN [1] T/F Should we print sequentially what is going on?
              double* value           // OUT [1]
){
  // Auxiliary variables
  int j, k, l, y, t, t2, auxk;
  int coord, Ycoord;
  double eta, auxY;
  int shift_nY;
  double* ptau;
  double* pbeta_num; 
  double* pbeta_poi; 
  double* pbeta_bin;
  double* pbeta_ord;
  double* pbeta_cat;
  double* pc;
  double etak, etak_1, qk;
  double pk, pk_1;
  int add_to_c;
  int cumKcat;
  int maxKcat = 0;
  for(y = 0; y < nY[4]; y++){
    if(maxKcat < Kcat[y]){
      maxKcat = Kcat[y];
    }
  }
  
  double* etas;
  etas = (double*)malloc(maxKcat * sizeof(double));
  //double etas[maxKcat];
  double log_sum_exp;
  
  
  *value = 0.0;
  
  
  if(spec[0]){
    // tau_num is g-specific
    ptau = (*last).tau_num + *g*dims[2];
  }else{
    ptau = (*last).tau_num;
  }
  
  // group-specific beta coeficients
  pbeta_num = (*last).beta_num + *g*dims[1];
  pbeta_poi = (*last).beta_poi + *g*dims[6];
  pbeta_bin = (*last).beta_bin + *g*dims[8];
  pbeta_ord = (*last).beta_ord + *g*dims[10];
  pbeta_cat = (*last).beta_cat + *g*dims[15];
  
  if(spec[1]){
    // a_ord and c_ord are g-specific
    pc = (*last).c_ord + *g*dims[11];
  }else{
    pc = (*last).c_ord;
  }
  
  for(j = 0; j < n_i[*i]; j++){
    coord = i_j[*i + *n * j]; // index of l-th observation of j-th individual within g-th cluster
    // Numeric variables
    shift_nY = 0;   // index within outcomes Y
    t = 0;          // index within beta_num
    for(y = 0; y < nY[0]; y++){
      Ycoord = coord + *N * shift_nY;
      // Update predictor
      eta  = predictor_num[4*(coord + *N * y)];     // fixed part of predictor, which is stable
      eta += predictor_num[2 + 4*(coord + *N * y)]; // random part
      eta += predictor_num[3 + 4*(coord + *N * y)]; // offset part
      // group-specific part
      for(l = 0; l < ngrp[shift_nY]; l++){
        // group-specific part of predictor, that needs to be updated wrt beta_num of different g
        eta += pbeta_num[t] * X[coord + *N * FormulaG[t]]; 
        t++;
      }
      
      if(isYna[Ycoord] & spec[4]){
        auxY = (*last).naY[isYna_inv[Ycoord] + *g * dims[32]] - eta;
      }else{
        auxY = Y[Ycoord] - eta;
      }
      
      *value -= 0.5 * ptau[y] * auxY * auxY; // subtracting the square
      *value += 0.5 * log(ptau[y]); // adding log of precision
      
      if(*print){
        printf("\nlogpYi_g: Nums[%d], i = %d, j = %d, value = %f", y, *i, j, *value);
        fflush(stdout);
      }
    
      shift_nY++;
    }
    
    t2 = t; // where numeric ended in FormulaG
    t = 0;          // index within beta_poi
    for(y = 0; y < nY[1]; y++){
      Ycoord = coord + *N * shift_nY;
      // Update predictor
      eta  = predictor_poi[4*(coord + *N * y)];     // fixed part of predictor, which is stable
      eta += predictor_poi[2 + 4*(coord + *N * y)]; // random part
      eta += predictor_poi[3 + 4*(coord + *N * y)]; // offset part
      // group-specific part
      for(l = 0; l < ngrp[shift_nY]; l++){
        // group-specific part of predictor, that needs to be updated wrt beta_poi of different g
        eta += pbeta_poi[t] * X[coord + *N * FormulaG[t2]]; 
        t++; t2++;
      }
      
      if(isYna[Ycoord] & spec[4]){
        auxY = (*last).naY[isYna_inv[Ycoord] + *g * dims[32]];
      }else{
        auxY = Y[Ycoord];
      }
      
      *value += auxY*eta - exp(eta);
      
      if(*print){
        printf("\nlogpYi_g: Pois[%d], i = %d, j = %d, value = %f", y, *i, j, *value);
        fflush(stdout);
      }
      
      shift_nY++;
    }
    
    
    // Binary variables
    t = 0;
    for(y = 0; y < nY[2]; y++){
      Ycoord = coord + *N * shift_nY;
      // Update predictor
      eta  = predictor_bin[4*(coord + *N * y)];     // fixed part of predictor, which is stable
      eta += predictor_bin[2 + 4*(coord + *N * y)]; // random part
      eta += predictor_bin[3 + 4*(coord + *N * y)]; // offset part
      // group-specific part
      for(l = 0; l < ngrp[shift_nY]; l++){
        // group-specific part of predictor, that needs to be updated wrt beta_bin of different g
        eta += pbeta_bin[t] * X[coord + *N * FormulaG[t2]]; 
        t++; t2++;
      }
      
      if(isYna[Ycoord] & spec[4]){
        auxY = (*last).naY[isYna_inv[Ycoord] + *g * dims[32]];
      }else{
        auxY = Y[Ycoord];
      }
      
      *value += auxY * eta - log(1+exp(eta));
      
      if(*print){
        printf("\nlogpYi_g: Bins[%d], i = %d, j = %d, value = %f", y, *i, j, *value);
        fflush(stdout);
      }
      
      shift_nY++;
    }
    
    // Ordinal variables
    add_to_c = 0;
    t = 0;
    for(y = 0; y < nY[3]; y++){
      Ycoord = coord + *N * shift_nY;
      // Update predictor
      eta  = predictor_ord[4*(coord + *N * y)];     // fixed part of predictor, which is stable
      eta += predictor_ord[2 + 4*(coord + *N * y)]; // random part
      eta += predictor_ord[3 + 4*(coord + *N * y)]; // offset part
      // group-specific part
      for(l = 0; l < ngrp[shift_nY]; l++){
        // group-specific part of predictor, that needs to be updated wrt beta_ord of different g
        eta += pbeta_ord[t] * X[coord + *N * FormulaG[t2]]; 
        t++; t2++;
      }
      
      if(isYna[Ycoord] & spec[4]){
        k = (*last).naY[isYna_inv[Ycoord] + *g * dims[32]];
      }else{
        k = Y[Ycoord]; // double to int, requires values in {0, ..., Kord[y]}
      }

      // Subtract the intercept for the category within
      if(k == 0){
        pk_1 = 1.0;
      }else{
        etak_1 = eta - pc[add_to_c + k-1];
        pk_1 = logit_inv(etak_1);
      }
      if(k == Kord[y]){
        pk = 0.0;
      }else{
        etak = eta - pc[add_to_c + k];
        pk = logit_inv(etak);
      }
      qk = pk_1 - pk;
      *value += log(qk);
      
      if(*print){
        printf("\nlogpYi_g: Ords[%d], i = %d, j = %d, value = %f", y, *i, j, *value);
        fflush(stdout);
      }
      
      shift_nY++;
      add_to_c += Kord[y];
    }
    
    // Categorical variables
    cumKcat = 0;
    t = 0;
    for(y = 0; y < nY[4]; y++){
      Ycoord = coord + *N * shift_nY;
      // Update predictor
      log_sum_exp = 1.0;
      for(k = 0; k < Kcat[y]; k++){
        etas[k]  = predictor_cat[4*(coord + *N * cumKcat)];     // fixed part of predictor, which is stable
        etas[k] += predictor_cat[2 + 4*(coord + *N * cumKcat)]; // random part
        etas[k] += predictor_cat[3 + 4*(coord + *N * cumKcat)]; // offset part
        // group-specific part
        for(l = 0; l < ngrp[shift_nY]; l++){
          etas[k] += pbeta_cat[t] * X[coord + *N * FormulaG[l + cumngrp[shift_nY]]];
          t++;
        }
        log_sum_exp += exp(etas[k]);
        cumKcat++; // updates by 1 to another set of predictors
      }
      
      if(isYna[Ycoord] & spec[4]){
        auxk = (*last).naY[isYna_inv[Ycoord] + *g * dims[32]];
      }else{
        auxk = Y[Ycoord]; // double to int, requires values in {0, ..., Kcat[y]}
      }
      
      for(k = 0; k < Kcat[y]; k++){
        if(auxk == k){
          *value += etas[k]; // if one of first K-1 categories, add corresponding predictor
        }
      }
      // in ALL cases, even Y==K-1 (last category), subtract log_sum_exp
      *value -= log(log_sum_exp);
      
      if(*print){
        printf("\nlogpYi_g: Cats[%d], i = %d, j = %d, value = %f", y, *i, j, *value);
        fflush(stdout);
      }
      
      shift_nY++;
    } // end of categorical variables
    
  } // end of for j in 1:n_i
  
  free(etas);
   
}