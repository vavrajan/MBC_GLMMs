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


### Source implemented R functions


### Loading the post-processed image
trial <- "" # or "_burnin" for shorter 
load(file = paste0(RDATA, "tutorial_image03", trial, ".RData"))


###-----------------------------------------------------------------------------
###     Not pre-processed data (only based on mcmc$all)
###-----------------------------------------------------------------------------
### Take the sampled U indicators
threshold <- 0.5 # if not exceeded, unit remains unclassified
n
n <- mcmc$settings["U", "dims"]

# posterior mode of Gplus=number of nonemtpy components
modeGplus <- as.numeric(tapply(mcmc$all$Gplus, mcmc$all$chain, 
                               function(x){names(table(x))[which.max(table(x))]}))
mcmc$modeGplus

clusters <- list()
for(ch in 1:Nchains){
  TAB <- table(unlist(mcmc$all[mcmc$all$chain == ch,grep("U\\[",colnames(mcmc$all))])+1)
  TAB <- TAB[order(TAB, decreasing = T)]
  clusters[[ch]] <- as.numeric(names(TAB)[1:modeGplus[ch]])
}
TAB # clusters are labelled from 0 to mcmc$G-1
clusters # the labels of non-empty clusters in each chain 

Us <- mcmc$all[,c("chain", paste0("U[",1:n,"]"))]
clusteringU <- matrix(-1, nrow = as.numeric(n), ncol = mcmc$Nchains)
certaintyU <- matrix(-1, nrow = as.numeric(n), ncol = mcmc$Nchains)
for(i in 1:n){
  TAB <- table(Us[,i+1]+1, Us$chain)
  # the cluster is the one with the highest frequency
  clusteringU[i,] <- as.numeric(rownames(TAB)[apply(TAB, 2, which.max)])
  # the highest frequency ratio
  certaintyU[i,] <- apply(TAB/mcmc$BM, 2, max)
}
# units with frequency ratio smaller than threshold remain unclassified (group 0)
clusteringU[certaintyU < threshold] <- 0

# proposed clustering by each chain separately
for(ch in 1:mcmc$Nchains){
  data[,paste0("clustering",ch)] <- clusteringU[data[,Id],ch]
  data[,paste0("certainty",ch)] <- certaintyU[data[,Id],ch]
}

data1 <- data[data$j==1,]
# Confusion matrices for each chain separately:
table(data1$clustering1, data1$g)
table(data1$clustering2, data1$g)
table(data1$clustering3, data1$g)
table(data1$clustering4, data1$g)
# We can also see the suitable permutation for labels 




###-----------------------------------------------------------------------------
###     Pre-processed data (based on mcmc$UGplus)
###-----------------------------------------------------------------------------
### Take the sampled and already permuted U indicators
dim(mcmc$UGplus[[1]])
mcmc$UGplus[[1]][1:5,1:10]

clusteringU <- matrix(-1, nrow = as.numeric(n), ncol = mcmc$Nchains)
certaintyU <- matrix(-1, nrow = as.numeric(n), ncol = mcmc$Nchains)
for(ch in 1:Nchains){
  for(i in 1:n){
    TAB <- table(mcmc$UGplus[[ch]][,i]+1)
    # the cluster is the one with the highest frequency
    clusteringU[i,ch] <- as.numeric(rownames(TAB)[which.max(TAB)])
    # the highest frequency ratio
    certaintyU[i,ch] <- max(TAB/mcmc$BM)
  }
}
# units with frequency ratio smaller than threshold remain unclassified (group 0)
clusteringU[certaintyU < threshold] <- 0

# And now continue the same way as before...
# proposed clustering by each chain separately
for(ch in 1:mcmc$Nchains){
  data[,paste0("clustering",ch)] <- clusteringU[data[,Id],ch]
  data[,paste0("certainty",ch)] <- certaintyU[data[,Id],ch]
}

data1 <- data[data$j==1,]
# Confusion matrices for each chain separately:
table(data1$clustering1, data1$g)
table(data1$clustering2, data1$g)
table(data1$clustering3, data1$g)
table(data1$clustering4, data1$g)

