###-----------------------------------------------------------------------------
###     Clustering the units based on sampled allocation indicators U
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
library("mvtnorm")
library("gaussquad")
library("polynom")
library("orthopolynom")
library("colorspace")

### Source implemented R functions
source(paste0(RROOT, "GenerateData.R"))
source(paste0(RROOT, "pUig_dev.R"))
source(paste0(RROOT, "FromCtoList.R"))
source(paste0(RROOT, "FromListtoC.R"))
source(paste0(RROOT, "FromCtoMatrix.R"))
source(paste0(RROOT, "FromMatrixtoC.R"))
source(paste0(RROOT, "FromListtoMatrix.R"))
source(paste0(RROOT, "plotting_functions.R"))

### Source implemented C functions
dyn.load(paste0(CROOT, "pUig_dev.dll"))
#dyn.unload(paste0(CROOT, "pUig_dev.dll"))


### Loading the post-processed image
trial <- "" # or "_burnin" for shorter 
load(file = paste0(RDATA, "tutorial_image01", trial, ".RData"))


### Generate new dataset according to the same model as in tutorial_01
### ... or use the same dataset to explore posterior distribution of classification probabilities
#ldatanew <- ldata
nnew <- 20 # number of newly observed units
n_i <- 4        # number of observations per unit
G <- 2          # number of clusters (GenerateSimData is adjusted only for cases G = 2 or 3)
timepar <- "mix" 
catpar <- "no"  
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

set.seed(27818)
ldatanew <- GenerateData(n = nnew, n_i = n_i,
                         G = G, timepar = timepar,
                         catpar = catpar, catsubjfix = catsubjfix,
                         xin1 = xin1, Sigma = Sigma,
                         Kord = Kord, Kcat = Kcat, kspec_bi_cat = kspec_bi_cat)
datanew <- ldatanew$data
head(datanew)
dim(datanew)

## Id variable, matrix of outcomes, matrix of covariates
Id <- "i"
Ynew <- datanew[,c(Id, Ys)]
Xnew <- datanew[,c(Id, "x", paste0("bs",1:spldim), "f")]

## Create missingness (just 1 value per outcome)
set.seed(123456789)
naijnew <- sample(1:dim(Ynew)[1], 1*length(Ys))

natruenew <- Ynew[naij[1], "ynum"]
Ynew[naijnew[1], "ynum"] <- NA

natruenew <- c(natruenew, Ynew[naij[2], "ypoi"])
Ynew[naijnew[2], "ypoi"] <- NA

natruenew <- c(natruenew, Ynew[naij[3], "ybin"])
Ynew[naijnew[3], "ybin"] <- NA

natruenew <- c(natruenew, Ynew[naij[4], "yord"])
Ynew[naijnew[4], "yord"] <- NA

natruenew <- c(natruenew, Ynew[naij[5], "ycat"])
Ynew[naijnew[5], "ycat"] <- NA

Ynew

### Serial code:

system.time(
  probs <- pUig_dev(mcmc, Ynew, Xnew, Id,
                    dynamic_prob = F, # should the probabilities be calculated dynamically
                    start = 1,
                    end = 10000,
                    thin = 1000,
                    ChainsToUse = 1:4,
                    tolerance = 1e-7, # tolerance when computing norm of the shift within Newton-Raphson
                    maxiter = 100, # maximum iterations allowed during Newton-Raphson update of proposal distribution
                    maxnrep = 20, # maximum number of repetitions when N-R fails
                    NGQ = 1 # number of points for Adaptive Gaussian Quadrature
                    # default value 1 corresponds to Laplacian approximation
                    # do not use values higher than 7, since there are NGQ^totnran summands
  )
)

### Different chains in parallel:
pUigChainsParallel <- function(chain = 1, 
                               mcmc = mcmc,
                               Ynew = Ynew, Xnew = Xnew, Id = Id,
                               dynamic_prob = F,
                               start = mcmc$B+1, end = mcmc$M, thin = 1,
                               tolerance = 1e-7,
                               maxiter = 100,
                               maxnrep = 20,
                               NGQ = 1,
                               RROOT = RROOT, CROOT = CROOT){
  
  ### Libraries
  library("mvtnorm")
  library("gaussquad")
  library("polynom")
  library("orthopolynom")
  
  ### Source implemented R functions
  source(paste0(RROOT, "pUig_dev.R"))
  source(paste0(RROOT, "FromCtoList.R"))
  source(paste0(RROOT, "FromListtoC.R"))
  source(paste0(RROOT, "FromCtoMatrix.R"))
  source(paste0(RROOT, "FromMatrixtoC.R"))
  source(paste0(RROOT, "FromListtoMatrix.R"))
  
  ### Source implemented C functions
  dyn.load(paste0(CROOT, "pUig_dev.dll"))
  #dyn.unload(paste0(CROOT, "pUig_dev.dll"))
  #dyn.load(paste0(CROOT, "pUig_devt.so"))
  
  print(paste0("Printing chain ", chain))

  RET <- pUig_dev(mcmc = mcmc,
                  Y = Ynew, X = Xnew, Id = Id,
                  dynamic_prob = dynamic_prob,
                  start = start,
                  end = end,
                  thin = thin,
                  ChainsToUse = chain,
                  tolerance = tolerance,
                  maxiter = maxiter,
                  maxnrep = maxnrep,
                  NGQ = NGQ)
  return(RET)
}

library("parallel")
## In parallel
# Starting cluster
myCluster <- makeCluster(mcmc$Nchains, nnodes = mcmc$Nchains) 

# Sampling
sampletime <-
  system.time(
    lprobs <- parLapply(cl = myCluster,
                       X = 1:mcmc$Nchains,
                       fun = pUigChainsParallel,
                       mcmc = mcmc,
                       Ynew = Ynew, Xnew = Xnew, Id = Id,
                       dynamic_prob = F,
                       start = 1, end = 10000, thin = 1,
                       tolerance = 1e-7,
                       maxiter = 100,
                       maxnrep = 20,
                       NGQ = 2,
                       RROOT = RROOT, CROOT = CROOT)
  )
# Ending cluster
stopCluster(myCluster)

# Stacking into one matrix
probs <- lprobs[[1]]
if(Nchains > 1){
  for(chain in 2:Nchains){
    lprobs[[chain]][,"chain"] <- rep(chain, dim(lprobs[[chain]])[1])
    probs <- rbind(probs, lprobs[[chain]])
  }
}

save(probs, file = paste0(RDATA, "tutorial_probs.RData"))



colnames(probs)

### Traceplots using implemented functions
pmcmc <- list()
pmcmc$all <- probs
pmcmc$all$logdev <- log(pmcmc$all$deviance)
pmcmc$howsave <- "matrix"
pmcmc$B = mcmc$B
pmcmc$M = mcmc$M
pmcmc$BM = mcmc$BM
pmcmc$Nchains = mcmc$Nchains
pmcmc$spec = mcmc$spec
pmcmc$settings = mcmc$settings
pmcmc$Nums = mcmc$Nums
pmcmc$Pois = mcmc$Pois
pmcmc$Bins = mcmc$Bins
pmcmc$Ords = mcmc$Ords
pmcmc$Cats = mcmc$Cats

figthin = max(round(M/1000), 1)
figthin = 1
chtplt = 1:mcmc$Nchains

## Deviance
cairo_pdf(paste0(FIG, "traceplots_deviance", trial, "_Gmax_", Gmax, "_nch_", Nchains, "_M_", M, "_B_", B, ".pdf"),
          width = 6, height = 5)
{
  par(mfrow = c(1,1), mar = c(4,4,1,1))
  plot.traceplots(mcmc = pmcmc, thin = figthin, B = 0, #M = 5,
                  what = "deviance", labcex = 1,
                  whatLAB = "Deviance", ChainsToPlot = chtplt)
}
dev.off()

## Log-deviance
cairo_pdf(paste0(FIG, "traceplots_logdev", trial, "_Gmax_", Gmax, "_nch_", Nchains, "_M_", M, "_B_", B, ".pdf"),
          width = 6, height = 5)
{
  par(mfrow = c(1,1), mar = c(4,4,1,1))
  plot.traceplots(mcmc = pmcmc, thin = figthin, B = 0, #M = 5,
                  what = "logdev", labcex = 1,
                  whatLAB = "log-deviance", ChainsToPlot = chtplt)
}
dev.off()

## Individual contirubitons to deviance
cairo_pdf(paste0(FIG, "traceplots_dev_i", trial, "_Gmax_", Gmax, "_nch_", Nchains, "_M_", M, "_B_", B, ".pdf"),
          width = 6, height = 5)
{
  par(mfrow = c(5,4), mar = c(4,4,1,1))
  for(i in 1:nnew){
    plot.traceplots(mcmc = pmcmc, thin = figthin, B = 0, #M = 5,
                    what = "dev_i", dimspec = i, labcex = 1, 
                    whatLAB = paste0("Dev[",i,"]"), ChainsToPlot = chtplt)
  }
}
dev.off()

### Traceplots
par(mfrow = c(1,2))
plot(log(probs$deviance), type = "l")
hist(log(probs$deviance))
par(mfrow = c(2,5))
for(i in 1:nnew){
  plot(probs[,paste0("dev_i[",i,"]")], type = "l",
       ylab = "Deviance", main = paste0("Contribution of i=",i))
}
par(mfrow = c(2,5))
for(i in 1:nnew){
  hist(probs[,paste0("dev_i[",i,"]")],
       xlab = "Deviance", main = paste0("Contribution of i=",i))
}

### Probabilities
ps <- probs[,grep("pUig", colnames(probs))]
psum <- summary(mcmc(ps))
probmean <- matrix(psum$statistics[,"Mean"], ncol = mcmc$G, nrow = nnew, byrow = T)
probmean
rowSums(probmean)
i <- 8
par(mfrow = c(2,5))
for(g in 1:mcmc$G){
  hist(probs[,paste0("pUig_int[",i,",",g,"]")], xlim = c(0,1))
}

threshold <- 0.5
clustering <- matrix(-1, nrow = as.numeric(nnew), ncol = mcmc$Nchains)
certainty <- matrix(-1, nrow = as.numeric(nnew), ncol = mcmc$Nchains)
clustering[,2] <- apply(probmean, 1, which.max)
certainty[,2] <- apply(probmean, 1, max)
# units with frequency ratio smaller than threshold remain unclassified (group 0)
clustering[certainty < threshold] <- 0
datanew$clustering2 <- clustering[,2]

datanew1 <- datanew[datanew$j==1,]
table(datanew1$clustering2, datanew1$g)
# row 0: unclassified 





#### How to parallelize computations from one long chain?
# function for parallelization of computation one chain into several balanced chunks
pUigStatesParallel <- function(run = 1, 
                               ncores = 8,
                               mcmc = mcmc,
                               Ynew = Ynew, Xnew = Xnew, Id = Id,
                               dynamic_prob = F,
                               start = mcmc$B+1, end = mcmc$M, thin = 1,
                               chain = 1,
                               tolerance = 1e-7,
                               maxiter = 25,
                               maxnrep = 10,
                               NGQ = 1,
                               RROOT = RROOT, CROOT = CROOT){
  
  ### Libraries
  library("mvtnorm")
  library("gaussquad")
  library("polynom")
  library("orthopolynom")
  
  ### Source implemented R functions
  source(paste0(RROOT, "pUig_dev.R"))
  source(paste0(RROOT, "FromCtoList.R"))
  source(paste0(RROOT, "FromListtoC.R"))
  source(paste0(RROOT, "FromCtoMatrix.R"))
  source(paste0(RROOT, "FromMatrixtoC.R"))
  source(paste0(RROOT, "FromListtoMatrix.R"))
  
  ### Source implemented C functions
  dyn.load(paste0(CROOT, "pUig_dev.dll"))
  #dyn.unload(paste0(CROOT, "pUig_dev.dll"))
  #dyn.load(paste0(CROOT, "pUig_devt.so"))
  
  print(paste0("Printing chain ", chain))
  
  states <- seq(start, end, by = thin)
  chunklength <- floor(length(states)/ncores)
  remainders <- length(states)-chunklength*ncores
  chunklengths <- rep(chunklength, ncores) + c(rep(1,remainders), rep(0, ncores-remainders))
  chunkcumsum <- c(0,cumsum(chunklengths))
  run_start <- states[chunkcumsum[run]+1]
  run_end <- states[chunkcumsum[run+1]]
  
  RET <- pUig_dev(mcmc = mcmc,
                  Y = Ynew, X = Xnew, Id = Id,
                  dynamic_prob = dynamic_prob,
                  start = run_start,
                  end = run_end,
                  thin = thin,
                  ChainsToUse = chain,
                  tolerance = tolerance,
                  maxiter = maxiter,
                  maxnrep = maxnrep,
                  NGQ = NGQ)
  return(RET)
}



library("parallel")
## In parallel
ncores <- 8
chain <- 1
# Starting cluster
myCluster <- makeCluster(ncores, nnodes = ncores) 

# Sampling
sampletime <-
  system.time(
    lprobs <- parLapply(cl = myCluster,
                        X = 1:ncores,
                        fun = pUigStatesParallel,
                        ncores = ncores,
                        mcmc = mcmc,
                        Ynew = Ynew, Xnew = Xnew, Id = Id,
                        dynamic_prob = F,
                        start = 1, end = 10000, thin = 1,
                        chain = chain,
                        tolerance = 1e-7,
                        maxiter = 25,
                        maxnrep = 10,
                        NGQ = 1,
                        RROOT = RROOT, CROOT = CROOT)
  )
# Ending cluster
stopCluster(myCluster)

# Stacking into one matrix
probs <- lprobs[[1]]
for(core in 2:ncores){
  lprobs[[core]][,"chain"] <- rep(chain, dim(lprobs[[core]])[1])
  probs <- rbind(probs, lprobs[[core]])
}
dim(probs)
probs[,1:10]
plot(log(probs$deviance), type = "l")




# #probsch <- probs
# par(mfrow = c(1,2))
# plot(log(probs$deviance), type = "l")
# plot(log(probsch$deviance[9001:10000]), type = "l")
