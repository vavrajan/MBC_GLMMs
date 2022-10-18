###-----------------------------------------------------------------------------
###     Compiling C functions needed for function_Metropolis_within_Gibbs_MBC_NumPoiBinOrdCat.R
###-----------------------------------------------------------------------------

### Purpose of this script is to compile used C functions
### into ".dll" file that is loaded within R functions


# ROOT <- "C:/Users/Vavra/Desktop/SKOLA/PhD studium/Disertacka/MBC_GLMMs/"
ROOT <- "path/to/your/directory/MBC_GLMMs/"
### set the paths to subdirectories
CROOT <- paste0(ROOT, "Cfun/")  # C functions
RROOT <- paste0(ROOT, "Rfun/")  # R functions
RDATA <- paste0(ROOT, "RData/") # for saved RData files
FIG <- paste0(ROOT, "Figures/") # for saved figures
TAB <- paste0(ROOT, "Tables/")  # for saved tables

setwd(CROOT)

### Compiling C functions
### libraries "withr" and "callr" needed
library("withr")
library("callr")

dyn.unload("Metropolis_within_Gibbs_MBC_NumPoiBinOrdCat.dll")

out <- rcmd(cmd = "SHLIB", cmdargs = c("Metropolis_within_Gibbs_MBC_NumPoiBinOrdCat.c",
                                       "calculate_predictor.c",
                                       "cholesky.c",
                                       "matrmult.c",
                                       "my_math.c",
                                       "myrdirichlet.c",
                                       "myrwishart.c",
                                       "gibbs_naY.c",
                                       "gibbs_beta_num.c",
                                       "gibbs_tau_num.c",
                                       "gibbs_InvQ.c",
                                       "gibbs_InvSigma.c",
                                       "gibbs_pUig.c",
                                       "metrgibbs_beta_poi.c",
                                       "metrgibbs_beta_bin.c",
                                       "metrgibbs_beta_ord.c",
                                       "metrgibbs_a_ord.c",
                                       "metrgibbs_beta_cat.c",
                                       "metrgibbs_bi.c",
                                       "metrgibbs_e0.c",
                                       "metropolis.c",
                                       "metropolis_le0.c",
                                       "newton_raphson.c",
                                       "newton_raphson_le0.c",
                                       "pdfs_derivatives_beta_poi.c",
                                       "pdfs_derivatives_beta_bin.c",
                                       "pdfs_derivatives_beta_ord.c",
                                       "pdfs_derivatives_api_ord.c",
                                       "pdfs_derivatives_beta_cat.c",
                                       "pdfs_derivatives_bi.c",
                                       "pdfs_derivatives_le0.c",
                                       "pdfs_last.c",
                                       "structures.h"))
cat(out$stderr)
cat(out$stdout)
out

dyn.load("Metropolis_within_Gibbs_MBC_NumPoiBinOrdCat.dll")
dyn.unload("Metropolis_within_Gibbs_MBC_NumPoiBinOrdCat.dll")
