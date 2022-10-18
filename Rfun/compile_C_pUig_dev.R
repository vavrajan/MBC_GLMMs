###-----------------------------------------------------------------------------
###     Compiling C functions needed for pUig_dev.R
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

dyn.unload("pUig_dev.dll")

out <- rcmd(cmd = "SHLIB", cmdargs = c("pUig_dev.c",
                                       "cholesky.c",
                                       "matrmult.c",
                                       "my_math.c",
                                       "newton_raphson_bi_dev.c",
                                       "pdfs_derivatives_bi_dev.c"))
cat(out$stderr)
cat(out$stdout)
out

dyn.load("pUig_dev.dll")
dyn.unload("pUig_dev.dll")
