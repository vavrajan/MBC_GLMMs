# MBC_GLMMs
Model-based clustering applied to several jointly modelled longitudinal outcomes of a mixed type by generalized linear mixed-effects models. 

### Clusterwise multivariate regression of mixed-type panel data 
#### J. Vávra, A. Komárek, B. Grün, and G. Malsiner-Walli 
#### Submitted, available as preprint
#### doi: https://doi.org/10.21203/rs.3.rs-1882841/v1

Implementation of MCMC sampler as a combination of C routines called from R.

Subdirectories:
* Cfun: implemented ".c" functions, ".h" headers, ".o" and ".dll" (for Windows)
* Rfun: implemented ".R" functions and scripts (including tutorials)
* Figures: saved figures
* RData: saved ".RData" files

Run "tutorial_0....R" in "Rfun" directory to know how to use the method.
But first you may need to compile the C function using "compile_....R" functions.


