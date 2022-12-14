################################################
###         Tutorial use of MBC_GLMMs        ###
###------------------------------------------###

### Purpose of this script is to demonstrate the use of implemented functions
### for estimating model for longitudinal data of 
##    numeric, 
##    Poisson count, 
##    binary, 
##    ordinal,
##    nominal nature
### for details of methodology, first read 
### J. Vávra, A. Komárek, B. Grün, and G. Malsiner-Walli. Clusterwise multivariate
### regression of mixed-type panel data. Submitted, available as preprint.
### doi: https://doi.org/10.21203/rs.3.rs-1882841/v1.

### First set up the working directory
### ROOT should contain the path to directory where (ending with "/")
### the following directories have been saved to
### - ./Rfun - contains all R functions
### - ./Cfun - contains all C functions

# ROOT <- "C:/Users/Vavra/Desktop/SKOLA/PhD studium/Disertacka/MBC_GLMMs/"
ROOT <- "path/to/your/directory/MBC_GLMMs/"
### set the paths to subdirectories
CROOT <- paste0(ROOT, "Cfun/")  # C functions
RROOT <- paste0(ROOT, "Rfun/")  # R functions
RDATA <- paste0(ROOT, "RData/") # for saved RData files
FIG <- paste0(ROOT, "Figures/") # for saved figures
TAB <- paste0(ROOT, "Tables/")  # for saved tables

setwd(ROOT)

###----------------------------------------------------
### Loading needed libraries, R and C functions
###----------------------------------------------------

### needed R libraries
library("coda")
library("mvtnorm")
library("HDInterval")
library("colorspace")

### loading R functions
source(paste0(RROOT, "GenerateData.R"))
source(paste0(RROOT, "function_Metropolis_within_Gibbs_MBC_NumPoiBinOrdCat.R"))
source(paste0(RROOT, "FromCtoList.R"))
source(paste0(RROOT, "FromListtoC.R"))
source(paste0(RROOT, "FromCtoMatrix.R"))
source(paste0(RROOT, "FromMatrixtoC.R"))
source(paste0(RROOT, "FromListtoMatrix.R"))

### loading C functions (in Windows ".dll")
# if needed, you can recompile it using "compile_Cfun.R"
dyn.load(paste0(CROOT, "Metropolis_within_Gibbs_MBC_NumPoiBinOrdCat.dll"))

### loading C functions (".so" version)
#dyn.load(paste0(CROOT, "Metropolis_within_Gibbs_MBC_NumPoiBinOrdCat.so"))



###-------------------------------------------------------
### Generating sample longitudinal data
###-------------------------------------------------------
# One outcome per each type
# will be generated by function GenerateData() from "GenerateData.R"
# which was also used to generate data in our simulation study (see the paper for details)
# joined through random intercepts for each outcome distributed from joint normal distribution
# given by matrix Sigma

## Parameters for simulation of the dataset
n <- 1000       # number of sample units
n_i <- 4        # number of observations per unit
G <- 2          # number of clusters (GenerateSimData is adjusted only for cases G = 2 or 3)
timepar <- "mix" # different parametrization types for different types of outcome
# other options:  "no" - no effect of time
#                 "parallel" - groups will be described by parallel lines
#                 "cross" - crossing lines
#                 "spline" - different spline parametrizations
#                 "mix" - different combinations of the previous ones for different outcomes
#                       - Num - spline, Poi - cross, Bin - parallel, Ord - spline, Cat - cross
# see the body of the function on how the betas are chosen...
catpar <- "no"  # effect of categorical covariate
# other options:  "fixed" - the same effect (0.7) for all groups
#                 "group" - differing among groups (-0.7, 0, 0.7)
catsubjfix = T  # should the values of categorical covariate be fixed for the same unit or alternate among observations?
xin1 = 1        # parameter xi, 1/xi describing the width of observational window
Kord = 5        # number of categories for Ord outcome
Kcat = 4        # number of categories for Cat outcome (careful, only default Kcat=4 implemented! It would require to set more betas...)
kspec_bi_cat = F# should the random effects be specific to each level of categorical outcome?
# If so, the the dimension of Sigma should be larger

## Covariance matrix(ces) Sigma
corr <- matrix(c(1, -0.5, -0.5, -0.4, -0.4,
                 -0.5, 1, 0.3, 0.4, 0.5,
                 -0.5, 0.3, 1, 0.2, 0.3,
                 -0.4, 0.4, 0.2, 1, 0.1,
                 -0.4, 0.5, 0.3, 0.1, 1), byrow = T, ncol = 5)
sdSigmas <- c(0.5, 0.5, 0.5, 0.5, 0.5)
Sigma <- diag(sdSigmas) %*% corr %*% diag(sdSigmas)

# If group-specific Sigma is desired, then create a list of matrices
# e.g.:
# 
# Sigma <- sdSigmas <- list()
# for(g in 1:G){
#   sdSigmas[[g]] <- rep(0.8 + (g-2)*0.4, 5)
#   Sigma[[g]] <- diag(sdSigmas[[g]]) %*% corr %*% diag(sdSigmas[[g]])

set.seed(31415)
ldata <- GenerateData(n = n, n_i = n_i,
                     G = G, timepar = timepar,
                     catpar = catpar, catsubjfix = catsubjfix,
                     xin1 = xin1, Sigma = Sigma,
                     Kord = Kord, Kcat = Kcat, kspec_bi_cat = kspec_bi_cat)
data <- ldata$data
head(data)

sd_num <- ldata$sd_num
beta_num <- ldata$beta_num
beta_poi <- ldata$beta_poi
beta_bin <- ldata$beta_bin
beta_ord <- ldata$beta_ord
c_ord <- ldata$c_ord
beta_cat <- ldata$beta_cat

# Outcomes:
Nums <- "ynum"
Pois <- "ypoi"
Bins <- "ybin"
Ords <- "yord"
Cats <- "ycat"
Ys <- c(Nums, Pois, Bins, Ords, Cats)
# pred_type - predictor for particular type
# b_type - random intercept for particular type

## Histogram of predictors
{
  layout(matrix(c(1,1,1,1,2,2,2,2,3,3,3,3,
                  4,4,4,5,5,5,6,6,6,7,7,7,
                  8,8,8,8,9,9,9,9,10,10,10,10), 3, 12, byrow = TRUE))
  par(mar = c(4,4,2,1))
  hist(data$pred_num, xlab = "predictor", main = "predictor num")
  hist(data$pred_poi, xlab = "predictor", main = "predictor poi")
  hist(data$pred_bin, xlab = "predictor", main = "predictor bin")
  hist(data$pred_ord - data$c_ord1, xlab = "predictor", main = "predictor ord - c1")
  hist(data$pred_ord - data$c_ord2, xlab = "predictor", main = "predictor ord - c2")
  hist(data$pred_ord - data$c_ord3, xlab = "predictor", main = "predictor ord - c3")
  hist(data$pred_ord - data$c_ord4, xlab = "predictor", main = "predictor ord - c4")
  hist(data$pred_cat1, xlab = "predictor", main = "predictor cat1")
  hist(data$pred_cat2, xlab = "predictor", main = "predictor cat2")
  hist(data$pred_cat3, xlab = "predictor", main = "predictor cat3")
}
## Outcomes vs group
{
  layout(matrix(c(1,1,1,2,2,2,
                  3,3,4,4,5,5), 2, 6, byrow = TRUE))
  par(mar = c(4,4,2,1))
  plot(ynum ~ factor(g), data = data, main = "ynum")
  plot(log(ypoi+1) ~ factor(g), data = data, main = "log(ypoi+1)")
  plot(factor(ybin) ~ factor(g), data = data, main = "ybin")
  plot(factor(yord) ~ factor(g), data = data, main = "yord")
  plot(factor(ycat) ~ factor(g), data = data, main = "ycat")
}
## Outcome vs. covariate x
{
  layout(matrix(c(1,1,1,2,2,2,
                  3,3,4,4,5,5), 2, 6, byrow = TRUE))
  par(mar = c(4,4,2,1))
  plot(ynum ~ x, data = data, pch = 16, col = data$g)
  plot(ypoi ~ x, data = data, pch = 16, col = data$g)
  plot(factor(ybin) ~ x, data = data)
  plot(factor(yord) ~ x, data = data)
  plot(factor(ycat) ~ x, data = data)
}
## Outcome vs. covariate x grouped
{
  layout(matrix(c(1,1,1,2,2,2,
                  3,3,4,4,5,5), 2, 6, byrow = TRUE))
  par(mar = c(4,4,2,1))
  plot(ynum ~ x, data = data, pch = 16, col = data$g)
  plot(ypoi ~ x, data = data, pch = 16, col = data$g)
  PlotCatVsxGrouped(y = "ybin", x = "x", g = "g", dat = data,
                        G = G, xbreaks = seq(0,1,by=0.1))
  PlotCatVsxGrouped(y = "yord", x = "x", g = "g", dat = data,
                        G = G, xbreaks = seq(0,1,by=0.1))
  PlotCatVsxGrouped(y = "ycat", x = "x", g = "g", dat = data,
                        G = G, xbreaks = seq(0,1,by=0.1))
}
## predictor vs. covariate x
{
  layout(matrix(c(1,1,1,1,2,2,2,2,3,3,3,3,
                  4,4,4,5,5,5,6,6,6,7,7,7,
                  8,8,8,8,9,9,9,9,10,10,10,10), 3, 12, byrow = TRUE))
  par(mar = c(4,4,2,1))
  plot(pred_num ~ x, data = data, pch = 16, col = data$g)
  plot(pred_poi ~ x, data = data, pch = 16, col = data$g)
  plot(pred_bin ~ x, data = data, pch = 16, col = data$g)
  plot(pred_ord-c_ord1 ~ x, data = data, pch = 16, col = data$g)
  plot(pred_ord-c_ord2 ~ x, data = data, pch = 16, col = data$g)
  plot(pred_ord-c_ord3 ~ x, data = data, pch = 16, col = data$g)
  plot(pred_ord-c_ord4 ~ x, data = data, pch = 16, col = data$g)
  plot(pred_cat1 ~ x, data = data, pch = 16, col = data$g)
  plot(pred_cat2 ~ x, data = data, pch = 16, col = data$g)
  plot(pred_cat3 ~ x, data = data, pch = 16, col = data$g)
}
## Spaghetti plot - numeric+count outcome
{
  par(mfrow = c(1,2), mar = c(4,4,2,1))
  plot(ynum ~ x, data = data, pch = 16, col = data$g, type = "n")
  for(i in 1:n){
    datai <- data[data$i == i,]
    lines(ynum ~ x, data = datai, lty = 1, lwd = 1, col = datai$g[1])
  }
  
  plot(ypoi ~ x, data = data, pch = 16, col = data$g, type = "n")
  for(i in 1:n){
    datai <- data[data$i == i,]
    lines(ypoi ~ x, data = datai, lty = 1, lwd = 1, col = datai$g[1])
  }
}
## Spaghetti plot - predictors
{
  layout(matrix(c(1,1,1,1,2,2,2,2,3,3,3,3,
                  4,4,4,5,5,5,6,6,6,7,7,7,
                  8,8,8,8,9,9,9,9,10,10,10,10), 3, 12, byrow = TRUE))
  par(mar = c(4,4,2,1))
  formulas <- sapply(c("pred_num ~ x",
                       "pred_poi ~ x",
                       "pred_bin ~ x",
                       "pred_ord-c_ord1 ~ x",
                       "pred_ord-c_ord2 ~ x",
                       "pred_ord-c_ord3 ~ x",
                       "pred_ord-c_ord4 ~ x",
                       "pred_cat1 ~ x",
                       "pred_cat2 ~ x",
                       "pred_cat3 ~ x"), as.formula)
  for(j in 1:10){
    plot(formulas[[j]], data = data, pch = 16, col = data$g, type = "n")
    for(i in 1:n){
      datai <- data[data$i == i,]
      lines(formulas[[j]], data = datai, lty = 1, lwd = 1, col = datai$g[1])
    }
  }
}

# spline details
degree <- 2
knots <- c(0, 0.25, 0.5, 0.75, 1)
inner <- knots[-c(1, length(knots))]
bound <- knots[c(1, length(knots))]
xgrid <- seq(0,1,by=0.01)
splB <- bs(xgrid, knots = inner, Boundary.knots = bound,
           degree = degree, intercept = FALSE)
spldim <- dim(splB)[2]



###-----------------------------------------------------------------------------
###   Preparing the input for MCMC sampler
###-----------------------------------------------------------------------------

## Id variable, matrix of outcomes, matrix of covariates
Id <- "i"
Y <- data[,c(Id, Ys)]
X <- data[,c(Id, "x", paste0("bs",1:spldim), "f")]

## Create missingness (just 2 values per outcome)
set.seed(123456789)
naij <- sample(1:dim(Y)[1], 2*length(Ys))

natrue <- Y[naij[1:2], "ynum"]
Y[naij[1:2], "ynum"] <- NA

natrue <- c(natrue, Y[naij[3:4], "ypoi"])
Y[naij[3:4], "ypoi"] <- NA

natrue <- c(natrue, Y[naij[5:6], "ybin"])
Y[naij[5:6], "ybin"] <- NA

natrue <- c(natrue, Y[naij[7:8], "yord"])
Y[naij[7:8], "yord"] <- NA

natrue <- c(natrue, Y[naij[9:10], "ycat"])
Y[naij[9:10], "ycat"] <- NA

## How do you wish to estimate your model?
timemodel <- "mix" # different parametrization types for different types of outcome
# other options:  "no" - no effect of time
#                 "parallel" - groups will be described by parallel lines
#                 "cross" - crossing lines
#                 "spline" - different spline parametrizations
#                 "mix" - different combinations of the previous ones for different outcomes
#                       - Num - spline, Poi - cross, Bin - parallel, Ord - spline, Cat - cross
# see the body of the function on how the betas are chosen...
catmodel <- "no"  # effect of categorical covariate
# other options:  "fixed" - the same effect (0.7) for all groups
#                 "group" - differing among groups (-0.7, 0, 0.7)

## Formula will be a list of formulas for each outcome and contains
# $fixed  - formula for effects common to all clusters
# $group  - formula for cluster-specific effects
# $random - formula for Id-specific effects (units are given by the Id variable)
# $offset - a colname from X to be used as an offset, "" means no offset

# if some effects are (accidentally) included both in fixed and group --> group-specificity is assumed

Formula <- list()
for(y in Ys){
  Formula[[y]] <- list()
  
  ## beta parameters fixed for all groups
  # betas for time variable x
  fixedtime <- ifelse(timemodel == "parallel", "x", "")
  if(timemodel == "mix"){
    fixedtime <- ifelse(y == "ybin", "x", "")
  }
  # betas for categorical covariate f
  fixedcat <- ifelse(catmodel == "fixed", "f", "")
  # all together
  form <- paste0(fixedtime,
                 ifelse(fixedtime!="" & fixedcat != "", " + ", ""),
                 fixedcat)
  Formula[[y]]$fixed <- as.formula(paste0("~ ",
                                          ifelse(form == "", "1", form)))
  
  ## beta parameters different among the groups
  # betas for time variable x
  grouptime <- switch(timemodel,
                      no = "",
                      parallel = "",
                      cross = "x",
                      spline = paste0("bs", 1:spldim, collapse = " + "),
                      mix = switch(y,
                                   "ynum" = paste0("bs", 1:spldim, collapse = " + "),
                                   "ypoi" = "x",
                                   "ybin" = "",
                                   "yord" = paste0("bs", 1:spldim, collapse = " + "),
                                   "ycat" = "x"))
  # betas for categorical covariate f
  groupcat <- ifelse(catmodel == "group", "f", "")
  # all together
  form <- paste0(grouptime,
                 ifelse(grouptime!="" & groupcat != "", " + ", ""),
                 groupcat)
  Formula[[y]]$group <- as.formula(paste0("~ ",
                                          ifelse(form == "", "1", form)))
  
  # random effect (just random intercept here) 
  Formula[[y]]$random <- as.formula("~ 1")
  # offset setting (just column name given, "" stands for no offset)
  Formula[[y]]$offset <- ""
}
Formula

### Tuning parameters, divided into integers and doubles
tuning <- list()
tuning$integer <- list(
  freq_proposal_update = 10, # how often to update the proposal distributions
  times_proposal = 10, # how many times to repeat the proposal/acceptance step
  maxiter = 25, # maximal number of iterations of Newton-Raphson method
  maxnrep = 100, # maximal number of repetitions of N-R in case it fails (starting from a different point)
  kspec_bi_cat = kspec_bi_cat # should the random effects for Cats be specific to each category level?
)
tuning$double <- list(
  const_proposal_beta_poi_fix = 1, # proposal variance inflation factor
  const_proposal_beta_poi = 1, # proposal variance inflation factor
  const_proposal_beta_bin_fix = 1, # proposal variance inflation factor
  const_proposal_beta_bin = 1, # proposal variance inflation factor
  const_proposal_beta_ord_fix = 0.5, # proposal variance inflation factor
  const_proposal_beta_ord = 0.5, # proposal variance inflation factor
  const_proposal_beta_cat_fix = 1, # proposal variance inflation factor
  const_proposal_beta_cat = 1, # proposal variance inflation factor
  const_proposal_a_ord = 0.5, # proposal variance inflation factor
  const_proposal_b = 0.5, # proposal variance inflation factor
  const_proposal_e0 = 1, # proposal variance inflation factor
  tolerance = 1e-7 # tolerance for convergence of Newton-Raphson
)

### Group-specificity of model parameters in specific order
spec <- c(T, T,
          F, F, # common covariance matrix of random effects
          F) # missing outcome values
names(spec) <- c("tau_num", "c_ord", "InvSigma", "InvQ", "naY")

### What parameters should be saved from the sampled chain
### ... and which should be disregarded
whatsave <- c(T, T,        # beta_num_fix, beta_num,
              T, T, T,     # tau_num, sd_num, var_num
              T, T,        # beta_poi_fix, beta_poi
              T, T,        # beta_bin_fix, beta_bin
              T, T,        # beta_ord_fix, beta_ord
              T, F, T,     # c_ord, a_ord, pi_ord, 
              T, T,        # beta_cat_fix, beta_cat
              T,T,T,T,F,   # InvSigma, Sigma, sdSigma, corSigma, detInvSigma
              F, F, F,     # InvQ, Q, detInvQ
              F,           # b
              T, T,        # w, ng
              T, F, T,     # loglik, pUig, U
              T, T,        # Gplus, e0
              T            # naY
)
names(whatsave) <- c("beta_num_fix", "beta_num",
                     "tau_num", "sd_num", "var_num",
                     "beta_poi_fix", "beta_poi",
                     "beta_bin_fix", "beta_bin",
                     "beta_ord_fix", "beta_ord",
                     "c_ord", "a_ord", "pi_ord",
                     "beta_cat_fix", "beta_cat",
                     "InvSigma", "Sigma", "sdSigma", "corSigma", "detInvSigma",
                     "InvQ", "Q", "detInvQ",
                     "b",
                     "w", "ng",
                     "loglik", "pUig", "U",
                     "Gplus", "e0",
                     "naY")
## How the chains should be saved?
# "matrix" large matrix, row = 1 iteration, column = 1 univariate parameter
# "list" parameters will be saved in a structured list in which it is easy to orient
howsave = "matrix"

### Hyperparameters of prior distributions
totnran <- dim(Sigma)[1]
param <- list(nu_0 = totnran + 1, # for I-W prior of Sigma
              nu_1 = totnran + 1, # for I-W prior of Q
              gamma_a=1, gamma_b=10, # gamma prior for taus
              # standard deviations for beta parameters depending on type and group-specificity
              sd_beta_num_fix = 1, 
              sd_beta_num = 1,
              sd_beta_poi_fix = 1,
              sd_beta_poi = 1,
              sd_beta_bin_fix = 1,
              sd_beta_bin = 1,
              sd_beta_ord_fix = 1,
              sd_beta_ord = 1,
              api_prior = 1,
              sd_beta_cat_fix = 1,
              sd_beta_cat = 1,
              # standard deviations for random effects, used only for sampling the initial values
              init_sd_b = 0.01,
              # gamma prior hyperparameters for e_0
              ae = 1, be = 100,
              # scale matrix for I-W prior of Q
              InvV = diag(0.01, totnran))



## function for parallel computation of different chains
SampleChainsParallel <- function(chain = 1, 
                                 Y = Y, XX = X, G = G,
                                 spec = spec, whatsave = whatsave,
                                 Nums = Nums, Pois = Pois, Bins = Bins, Ords = Ords, Cats = Cats, 
                                 Id = Id, Formula = Formula, 
                                 inits = inits, param = param,
                                 M = M, B = B, howsave = howsave,
                                 tuning = tuning,
                                 ROOT = ROOT, RROOT = RROOT, CROOT = CROOT){
  
  source(paste0(RROOT, "function_Metropolis_within_Gibbs_MBC_NumPoiBinOrdCat.R"))
  source(paste0(RROOT, "FromCtoList.R"))
  source(paste0(RROOT, "FromListtoC.R"))
  source(paste0(RROOT, "FromCtoMatrix.R"))
  source(paste0(RROOT, "FromMatrixtoC.R"))
  
  library("mvtnorm")
  
  dyn.load(paste0(CROOT, "Metropolis_within_Gibbs_MBC_NumPoiBinOrdCat.dll"))
  #dyn.unload(paste0(CROOT, "Metropolis_within_Gibbs_MBC_NumPoiBinOrdCat.dll"))
  #dyn.load(paste0(CROOT, "Metropolis_within_Gibbs_MBC_NumPoiBinOrdCat.so"))
  
  print(paste0("Printing chain ", chain))
  
  RET <- list()
  
  if(missing(inits)){
    RET <- Metropolis_within_Gibbs_MBC_NumPoiBinOrdCat(Y = Y, X = XX, G = G,
                                                       spec = spec, whatsave = whatsave,
                                                       Nums = Nums, Pois = Pois, Bins = Bins, Ords = Ords, Cats = Cats, 
                                                       Id = Id, Formula = Formula, 
                                                       param = param, howsave = howsave,
                                                       tuning = tuning,
                                                       M = M, B = B, Nchains = 1)
  }else{
    chinits <- list()
    chinits[[1]] <- inits[[chain]]
    #chinits[[1]] <- inits[[ifelse(chain<4, chain, 4)]]
    
    RET <- Metropolis_within_Gibbs_MBC_NumPoiBinOrdCat(Y = Y, X = XX, G = G,
                                                       spec = spec, whatsave = whatsave,
                                                       Nums = Nums, Pois = Pois, Bins = Bins, Ords = Ords, Cats = Cats,
                                                       Id = Id, Formula = Formula, 
                                                       inits = chinits,
                                                       param = param, howsave = howsave,
                                                       tuning = tuning,
                                                       M = M, B = B, Nchains = 1)
  }
  
  return(RET)
}

library("parallel")

### First inicialization, just to see whether it works
B <- 0         # length of the burn-in period
M <- 500       # length of the final sample
# in the end (B+M) iterations is recorded, but only the last M of them will be saved
Gmax <- 10     # maximal number of mixture components to be considered
Nchains = 4    # number of chains to be sampled 

# ## To do it in serial - chain by chain - simply use the following command
# ## benefit: you can see the percent progress bar + potential errors or problems
mcmc <- Metropolis_within_Gibbs_MBC_NumPoiBinOrdCat(Y = Y, X = X, G = Gmax,
                                            spec = spec, whatsave = whatsave,
                                            Nums = Nums, Pois = Pois, Bins = Bins, Ords = Ords, Cats = Cats,
                                            Id = Id, Formula = Formula,
                                            param = param, howsave = howsave,
                                            tuning = tuning,
                                            M = M, B = B, Nchains = Nchains)


## In parallel
# Starting cluster
myCluster <- makeCluster(Nchains, nnodes = Nchains) 

# Sampling
sampletime <-
  system.time(
    mcmcs <- parLapply(cl = myCluster,
                       X = 1:Nchains,
                       fun = SampleChainsParallel,
                       Y = Y, XX = X, G = Gmax,
                       spec = spec, whatsave = whatsave,
                       Nums = Nums, Pois = Pois, Bins = Bins, Ords = Ords, Cats = Cats,
                       Id = Id, Formula = Formula,
                       # by not giving the inits, the procedure will find its own
                       #inits = inits, 
                       param = param,
                       M = M, B = B, howsave = howsave,
                       tuning = tuning,
                       ROOT = ROOT, RROOT = RROOT, CROOT = CROOT)
  )
# Ending cluster
stopCluster(myCluster)

# the results are in lists (mcmcs), to merge them apply the following commands
# to have the same structure as mcmc above does:
if(howsave == "matrix"){
  mcmc <- mcmcs[[1]]
  if(Nchains > 1){
    for(chain in 2:Nchains){
      mcmcs[[chain]]$all[,"chain"] <- rep(chain, dim(mcmcs[[chain]]$all)[1])
      mcmc$all <- rbind(mcmc$all, mcmcs[[chain]]$all)
      mcmc$inits[[chain]] <- mcmcs[[chain]]$inits[[1]]
      mcmc$last[[chain]] <- mcmcs[[chain]]$last[[1]]
    }
  }
  mcmc$Nchains <- Nchains
}

if(howsave == "list"){
  mcmc <- mcmcs[[1]]
  if(Nchains > 1){
    for(chain in 2:Nchains){
      mcmc[[chain]] <- mcmcs[[chain]][[1]]
      mcmc$inits[[chain]] <- mcmcs[[chain]]$inits[[1]]
      mcmc$last[[chain]] <- mcmcs[[chain]]$last[[1]]
    }
  }
  mcmc$Nchains <- Nchains
}

# Saving last state for other sampling (to be used as initial values some next time)
inits <- list()
for(ch in 1:Nchains){
  inits[[ch]] <- mcmc$last[[ch]]
  inits[[ch]]$U <- inits[[ch]]$U+1
}
save(inits, file = paste0(RDATA, "tutorial_inits_burnin.RData"))

# Saving sampled chains
save(mcmc, file = paste0(RDATA, "tutorial_mcmc_burnin.RData"))

# Saving the whole image
save.image(file = paste0(RDATA, "tutorial_image01_burnin.RData"))



### Continue in sampling
B = 0
M = 10000

## In parallel
# Starting cluster
myCluster <- makeCluster(Nchains, nnodes = Nchains) 

# Sampling
sampletime <-
  system.time(
    mcmcs <- parLapply(cl = myCluster,
                       X = 1:Nchains,
                       fun = SampleChainsParallel,
                       Y = Y, XX = X, G = Gmax,
                       spec = spec, whatsave = whatsave,
                       Nums = Nums, Pois = Pois, Bins = Bins, Ords = Ords, Cats = Cats,
                       Id = Id, Formula = Formula,
                       # Now we supply the initial values
                       inits = inits, 
                       param = param,
                       M = M, B = B, howsave = howsave,
                       tuning = tuning,
                       ROOT = ROOT, RROOT = RROOT, CROOT = CROOT)
  )
# Ending cluster
stopCluster(myCluster)

# the results are in lists (mcmcs), to merge them apply the following commands
# to have the same structure as mcmc above does:
if(howsave == "matrix"){
  mcmc <- mcmcs[[1]]
  if(Nchains > 1){
    for(chain in 2:Nchains){
      mcmcs[[chain]]$all[,"chain"] <- rep(chain, dim(mcmcs[[chain]]$all)[1])
      mcmc$all <- rbind(mcmc$all, mcmcs[[chain]]$all)
      mcmc$inits[[chain]] <- mcmcs[[chain]]$inits[[1]]
      mcmc$last[[chain]] <- mcmcs[[chain]]$last[[1]]
    }
  }
  mcmc$Nchains <- Nchains
}

if(howsave == "list"){
  mcmc <- mcmcs[[1]]
  if(Nchains > 1){
    for(chain in 2:Nchains){
      mcmc[[chain]] <- mcmcs[[chain]][[1]]
      mcmc$inits[[chain]] <- mcmcs[[chain]]$inits[[1]]
      mcmc$last[[chain]] <- mcmcs[[chain]]$last[[1]]
    }
  }
  mcmc$Nchains <- Nchains
}

# Saving last state for other sampling (to be used as initial values some next time)
inits <- list()
for(ch in 1:Nchains){
  inits[[ch]] <- mcmc$last[[ch]]
  inits[[ch]]$U <- inits[[ch]]$U+1
}
save(inits, file = paste0(RDATA, "tutorial_inits.RData"))




# Saving sampled chains
save(mcmc, file = paste0(RDATA, "tutorial_mcmc.RData"))

# Saving the whole image
save.image(file = paste0(RDATA, "tutorial_image01.RData"))




 

