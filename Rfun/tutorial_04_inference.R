###-----------------------------------------------------------------------------
###     Inference based on the post-processed MCMC
###-----------------------------------------------------------------------------

# ROOT <- "C:/Users/Vavra/Desktop/SKOLA/PhD studium/Disertacka/MBC_GLMMs/"
ROOT <- "path/to/your/directory/MBC_GLMMs/"
### set the paths to subdirectories
CROOT <- paste0(ROOT, "Cfun/")  # C functions
RROOT <- paste0(ROOT, "Rfun/")  # R functions
RDATA <- paste0(ROOT, "RData/") # for saved RData files
FIG <- paste0(ROOT, "Figures/") # for saved figures
TAB <- paste0(ROOT, "Tables/")  # for saved tables

setwd(ROOT)

### Libraries
library("coda")
library("HDInterval")

### Source implemented R functions
source(paste0(RROOT, "FromListtoMatrix.R"))

### Loading the post-processed image
trial <- "" # or "_burnin" for shorter 
load(file = paste0(RDATA, "tutorial_image03", trial, ".RData"))

### Look at mcmc content
# if howsave="matrix" then matrix $all contains all chains and all sampled parameters
# if howsave="list" then it is saved as a structured list
# In the following it would be better to work with $all matrix
# If you do not have the all matrix, do the following command
if(mcmc$howsave=="list"){
  # Then transfer to howsave="matrix" -- create mcmc$all matrix
  mcmc <- FromListtoMatrix_settings(mcmc)
}

dim(mcmc$all)
mcmc$all[1:5, 1:10]

# Other interesting content of mcmc
mcmc$Formula
mcmc$lfixnames
mcmc$lgrpnames
mcmc$lrannames
mcmc$G          # Gmax actually - the maximal number of considered mixture components
mcmc$spec       # group-specificity of model parameters
mcmc$whatsave   # what quantities were and which were not saved
mcmc$totnran    # dimension of the Sigma matrix (total number of random effects)
mcmc$Kord       # number of categories for ordinal outcomes
mcmc$Kcat       # number of categories for categorical outcomes
mcmc$settings   # detailed description of model parameters and their specific attributes
sum(mcmc$isYna) # indicators of unobserved outcome values

## Outputs of post-processing
mcmc$modeGplus  # posterior mode of Gplus = number of non-empty components
mcmc$Ms         # number of states (from total M) which have Gplus equal to the posterior mode
mcmc$inds[[1]]  # indexes of those states (for chain 1)
mcmc$Mrhos      # number of states (from total Ms) which do not yield permutation after k-means procedure
dim(mcmc$clusterspec[[1]]) # matrix of cluster-specific parameters only (chain 1)
# rows with NA values are the ones not yielding a permutation after k-means (not suitable for inference)
# Other group-related parameters post-processed
head(mcmc$wGplus[[1]])      # mixture weights
head(mcmc$ngGplus[[1]])     # counts of units in each cluster
mcmc$UGplus[[1]][1:5,1:15]  # cluster allocations
mcmc$pUigGplus              # not saved in mcmc$all --> not in here


### Coda compatibility
## Un-processed $all data matrix
## This includes even parameters for empty components
codamcmc <- list()
lsumcodamcmc <- list()
thin = 1 # thinning parameter
# There is slight possibility that it will produce some NA values
# We also omit n allocation indicators U
for(ch in 1:Nchains){
  codamcmc[[ch]] <- mcmc$all[mcmc$all$chain == ch,
                             setdiff(2:dim(mcmc$all)[2], grep("U", colnames(mcmc$all)))]
  codamcmc[[ch]] <- mcmc(codamcmc[[ch]], thin = thin, start = mcmc$B+1, end = mcmc$BM)
  lsumcodamcmc[[ch]] <- summary(codamcmc[[ch]])
}
codamcmc <- as.mcmc.list(codamcmc)
sumcodamcmc <- summary(codamcmc)
## !! This summary combines all 4 chains together!!
## !!!Moreover!!! cluster labels do not have unified meaning
## Permute them first!!!
## Otherwise cluster specific parameters (those with "(g)" in the name) are 
## wrongly estimated
sumcodamcmc
lsumcodamcmc[[1]]
lsumcodamcmc[[2]]
lsumcodamcmc[[3]]
lsumcodamcmc[[4]]

save(sumcodamcmc, file = paste0(RDATA, "sumcodamcmc_all.RData"))
save(lsumcodamcmc, file = paste0(RDATA, "lsumcodamcmc_all.RData"))

## Example: focus on precision tau for numeric outcome
sumcodamcmc$statistics[grep("tau_num", rownames(sumcodamcmc$statistics)),]
lsumcodamcmc[[1]]$statistics[grep("tau_num", rownames(lsumcodamcmc[[1]]$statistics)),]
lsumcodamcmc[[2]]$statistics[grep("tau_num", rownames(lsumcodamcmc[[2]]$statistics)),]
lsumcodamcmc[[3]]$statistics[grep("tau_num", rownames(lsumcodamcmc[[3]]$statistics)),]
lsumcodamcmc[[4]]$statistics[grep("tau_num", rownames(lsumcodamcmc[[4]]$statistics)),]
# See? 3.7 precision corresponds to labels 5, 4, 6, 1 for chains 1,2,3,4, respectively
# Which is all WRONGLY merged in the total sumcodamcmc summary
# The same for quantile summaries:
sumcodamcmc$quantiles[grep("tau_num", rownames(sumcodamcmc$quantiles)),]
lsumcodamcmc[[1]]$quantiles[grep("tau_num", rownames(lsumcodamcmc[[1]]$quantiles)),]
lsumcodamcmc[[2]]$quantiles[grep("tau_num", rownames(lsumcodamcmc[[2]]$quantiles)),]
lsumcodamcmc[[3]]$quantiles[grep("tau_num", rownames(lsumcodamcmc[[3]]$quantiles)),]
lsumcodamcmc[[4]]$quantiles[grep("tau_num", rownames(lsumcodamcmc[[4]]$quantiles)),]
# true values:
1/sd_num[[2]]^2

# Solution:
# a) use the processed data (below)
# b) permute the meaning in each chain so that it corresponds
# c) use only summaries from one chain

## But the overall summary still works for non-specific parameters
sumcodamcmc$quantiles[grep("sdSigma", rownames(sumcodamcmc$quantiles)),]
lsumcodamcmc[[1]]$quantiles[grep("sdSigma", rownames(lsumcodamcmc[[1]]$quantiles)),]
lsumcodamcmc[[2]]$quantiles[grep("sdSigma", rownames(lsumcodamcmc[[2]]$quantiles)),]
lsumcodamcmc[[3]]$quantiles[grep("sdSigma", rownames(lsumcodamcmc[[3]]$quantiles)),]
lsumcodamcmc[[4]]$quantiles[grep("sdSigma", rownames(lsumcodamcmc[[4]]$quantiles)),]
# true values:
sdSigmas

sumcodamcmc$quantiles[grep("corSigma", rownames(sumcodamcmc$quantiles)),]
lsumcodamcmc[[1]]$quantiles[grep("corSigma", rownames(lsumcodamcmc[[1]]$quantiles)),]
lsumcodamcmc[[2]]$quantiles[grep("corSigma", rownames(lsumcodamcmc[[2]]$quantiles)),]
lsumcodamcmc[[3]]$quantiles[grep("corSigma", rownames(lsumcodamcmc[[3]]$quantiles)),]
lsumcodamcmc[[4]]$quantiles[grep("corSigma", rownames(lsumcodamcmc[[4]]$quantiles)),]
# true values:
corr
corrmeds <- sumcodamcmc$quantiles[grep("corSigma", rownames(sumcodamcmc$quantiles)),"50%"]
corrmed <- matrix(0, nrow = mcmc$totnran, ncol = mcmc$totnran)
corrmed[upper.tri(corrmed)] <- corrmeds
corrmed <- t(corrmed)
corrmed[upper.tri(corrmed)] <- corrmeds
diag(corrmed) <- rep(1, mcmc$totnran)
corrmed

## Equal-tailed (ET) confidence intervals are the 2.5% and 97.5% columns from $quantiles
qs <- paste0(c("2.5", "50", "97.5"), "%")
posteriorcorr <- list()
for(q in qs){
  corrq <- sumcodamcmc$quantiles[grep("corSigma", rownames(sumcodamcmc$quantiles)), q]
  posteriorcorr[[q]] <- matrix(0, nrow = mcmc$totnran, ncol = mcmc$totnran)
  posteriorcorr[[q]][upper.tri(posteriorcorr[[q]])] <- corrq
  posteriorcorr[[q]] <- t(posteriorcorr[[q]])
  posteriorcorr[[q]][upper.tri(posteriorcorr[[q]])] <- corrq
  diag(posteriorcorr[[q]]) <- rep(1, mcmc$totnran)
}
posteriorcorr
# true values:
corr

## Highest-posterior density (HPD) intervals
# use function hdi from library("HDInterval")
HPDs <- hdi(codamcmc)
HPDs[,grep("tau_num", colnames(HPDs))]
HPDs[,grep("corSigma", colnames(HPDs))]


### Post-processed group-specific parameters
clusterspecall <- matrix(0, nrow = 0, ncol = dim(mcmc$clusterspec[[1]])[2]-1)
codaclusterspec <- list()
lsumcodaclusterspec <- list()
thin = 1 # thinning parameter
# There is slight possibility that it will produce some NA values
# We also omit n allocation indicators U
for(ch in 1:Nchains){
  cnames <- colnames(mcmc$clusterspec[[ch]])[-1]
  codaclusterspec[[ch]] <- matrix(unlist(mcmc$clusterspec[[ch]][,2:dim(mcmc$clusterspec[[ch]])[2]]),
                                  ncol = dim(mcmc$clusterspec[[ch]])[2]-1, byrow = F)
  colnames(codaclusterspec[[ch]]) <- cnames
  clusterspecall <- rbind(clusterspecall, codaclusterspec[[ch]])
  codaclusterspec[[ch]] <- mcmc(codaclusterspec[[ch]], 
                                thin = thin, start = 1, end = dim(mcmc$clusterspec[[ch]])[1])
  lsumcodaclusterspec[[ch]] <- summary(codaclusterspec[[ch]])
}
# The following might not work due to different lengths of the chains
codaclusterspec <- as.mcmc.list(codaclusterspec)
sumcodaclusterspec <- summary(codaclusterspec)
# Hence, merge all chains together
codaclusterspecall <- mcmc(clusterspecall, 
                           thin = thin, start = 1, end = dim(clusterspecall)[1])
sumcodaclusterspec <- summary(codaclusterspecall)

## Example: focus on precision tau for numeric outcome
# Posterior mean and MC error
sumcodaclusterspec$statistics[grep("tau_num", rownames(sumcodaclusterspec$statistics)),]
lsumcodaclusterspec[[1]]$statistics[grep("tau_num", rownames(lsumcodaclusterspec[[1]]$statistics)),]
lsumcodaclusterspec[[2]]$statistics[grep("tau_num", rownames(lsumcodaclusterspec[[2]]$statistics)),]
lsumcodaclusterspec[[3]]$statistics[grep("tau_num", rownames(lsumcodaclusterspec[[3]]$statistics)),]
lsumcodaclusterspec[[4]]$statistics[grep("tau_num", rownames(lsumcodaclusterspec[[4]]$statistics)),]

# Posterior quantiles
sumcodaclusterspec$quantiles[grep("tau_num", rownames(sumcodaclusterspec$quantiles)),]
lsumcodaclusterspec[[1]]$quantiles[grep("tau_num", rownames(lsumcodaclusterspec[[1]]$quantiles)),]
lsumcodaclusterspec[[2]]$quantiles[grep("tau_num", rownames(lsumcodaclusterspec[[2]]$quantiles)),]
lsumcodaclusterspec[[3]]$quantiles[grep("tau_num", rownames(lsumcodaclusterspec[[3]]$quantiles)),]
lsumcodaclusterspec[[4]]$quantiles[grep("tau_num", rownames(lsumcodaclusterspec[[4]]$quantiles)),]

# HPD intervals
hdi(clusterspecall)[,paste0("tau_num_ynum(",1:2,")")]
hdi(codaclusterspec[[1]])[,paste0("tau_num_ynum(",1:2,")")]
hdi(codaclusterspec[[2]])[,paste0("tau_num_ynum(",1:2,")")]
hdi(codaclusterspec[[3]])[,paste0("tau_num_ynum(",1:2,")")]
hdi(codaclusterspec[[4]])[,paste0("tau_num_ynum(",1:2,")")]


### In the same way we could do even other group-related parameters
summary(mcmc(mcmc$wGplus[[1]]))
summary(mcmc(mcmc$ngGplus[[1]]))
# the following makes (sort of) sense only with G=Gplus=2
sumU <- summary(mcmc(mcmc$UGplus[[1]]))
sumU$statistics[,"Mean"]
