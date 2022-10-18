/*
 * Mathematical functions needed
 */

#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "structures.h"

double logit_inv(double x){
  double expx;
  
  if(x < 0){
    expx = exp(x);
    return(expx/(1.0+expx));
  }else{
    expx = exp(-x);
    return(1.0/(1.0+expx));
  }
  
  
}