###-----------------------------------------------------------------------------
###     Post-processing the sampled chains
###-----------------------------------------------------------------------------
# Label-switching check
# post-processing procedure by Fruhwirth, Malsiner-Walli

# ROOT <- "C:/Users/Vavra/Desktop/SKOLA/PhD studium/Disertacka/MBC_GLMMs/"
ROOT <- "path/to/your/directory/MBC_GLMMs/"
### set the paths to subdirectories
CROOT <- paste0(ROOT, "Cfun/")  # C functions
RROOT <- paste0(ROOT, "Rfun/")  # R functions
RDATA <- paste0(ROOT, "RData/") # for saved RData files
FIG <- paste0(ROOT, "Figures/") # for saved figures
TAB <- paste0(ROOT, "Tables/")  # for saved tables

setwd(ROOT)

### Source implemented R functions
source(paste0(RROOT, "plotting_functions.R"))
source(paste0(RROOT, "FromListtoMatrix.R"))
source(paste0(RROOT, "RestrictionToModeGplus.R"))
source(paste0(RROOT, "PermuteClusterSpec.R"))

# If you wish to inspect the burn-in period
load(file = paste0(RDATA, "tutorial_image01_burnin.RData"))
trial <- "_burnin"

# If you wish to inspect the final chains
load(file = paste0(RDATA, "tutorial_image01.RData"))
trial <- ""

### Gplus traceplot
#figthin = max(round(M/500), 1)
figthin = 1

par(mfrow = c(1,1))
plot.traceplots(mcmc = mcmc, thin = figthin, 
                what = "Gplus", labcex = 0.6, 
                whatLAB = "Gplus",
                ChainsToPlot = 1:mcmc$Nchains)

### First run the function RestrictionToModeGplus
# takes some time with M large....
mcmc <- RestrictionToModeGplus(mcmc)

# the most frequent Gplus value in each chain
mcmc$modeGplus
# number of states with modeGplus value
mcmc$Ms
# number of states (out of total Ms) not resulting in a permutation
mcmc$Mrhos
# (potentially reshuffled) cluster-specific parameters
# also reduced notation to 1:modeGplus
# different chains may have still permuted meaning compared to other chains
dim(mcmc$clusterspec[[1]])
head(mcmc$clusterspec[[1]])

# indexes of states with current Gplus=modeGplus
mcmc$inds[[1]]

# cluster related parameters permuted to fit clusterspec matrix
head(mcmc$wGplus[[1]])
head(mcmc$ngGplus[[1]])
head(mcmc$UGplus[[1]][,1:20])

### Additionally permute so that the meaning of the number of clusters
### is kept the same across all chains

## first we have to recognize which chain has which meaning, see the plots below

perm <- list()
#perm[[1]] <- c(2, 1)
#perm[[2]] <- c(1, 2)
#perm[[3]] <- c(2, 4, 1, 3)
perm[[1]] <- 1:mcmc$modeGplus[[ch]]
perm[[2]] <- 1:mcmc$modeGplus[[ch]]
perm[[3]] <- mcmc$modeGplus[[ch]]:1 # this chain has switched meaning of 1 and 2
perm[[4]] <- 1:mcmc$modeGplus[[ch]]

# function to permute meaning of cluster labels in mcmc according to perm
mcmc <- PermuteClusterSpec(mcmc, perm)



### You can check that the meanings of cluster labels are identical in all chains
### by plotting the traceplots or kernel density estimators

###---------
### FIGURES
###---------
figthin = max(round(M/500), 1)
#figthin = 1

## probabilites w
for(ch in 1:mcmc$Nchains){
  par(mfrow = c(mcmc$modeGplus[[ch]],4))
  for(g in 1:mcmc$modeGplus[[ch]]){
    plot.traceplots(mcmc = mcmc, thin = figthin, 
                    what = "w", dimspec = g, labcex = 0.6, 
                    whatLAB = paste0("w[",g,"]"),
                    ChainsToPlot = ch, restrictedtoGplus = T)
    plot.kerneldensity(mcmc = mcmc, thin = figthin, 
                       what = "w", dimspec = g, labcex = 0.6, 
                       whatLAB = paste0("w[",g,"]"),
                       ChainsToPlot = ch, restrictedtoGplus = T)
    plot.ECDF(mcmc = mcmc, thin = figthin, 
              what = "w", dimspec = g, labcex = 0.6, 
              whatLAB = paste0("w[",g,"]"),
              ChainsToPlot = ch, restrictedtoGplus = T)
    plot.ACF(mcmc = mcmc, thin = figthin, 
             what = "w", dimspec = g, labcex = 0.6, 
             whatLAB = paste0("w[",g,"]"), MaxLag = 30,
             ChainsToPlot = ch, restrictedtoGplus = T)
  }
}

totnran <- dim(mcmc$last[[1]]$InvSigma[[1]])[1]
totnran <- dim(mcmc$last[[1]]$InvSigma)[1]

## tau_num
for(ch in 1:mcmc$Nchains){
  for(y in Nums){
    for(g in 1:mcmc$modeGplus[[ch]]){
      plot.all(mcmc = mcmc, thin = figthin, gspec = g,
               what = "tau_num", yspec = y, MaxLag = 30,
               ChainsToPlot = ch, restrictedtoGplus = T)
    }
  }
}
# plot.all(mcmc = mcmc, thin = figthin, gspec = 2,
#          what = "tau_num", yspec = Nums[1], MaxLag = 30,
#          ChainsToPlot = 1:mcmc$Nchains, restrictedtoGplus = T)

## sd_num
for(ch in 1:mcmc$Nchains){
  for(y in Nums){
    for(g in 1:mcmc$modeGplus[[ch]]){
      plot.all(mcmc = mcmc, thin = figthin, gspec = g,
               what = "sd_num", yspec = y, MaxLag = 30,
               ChainsToPlot = ch, restrictedtoGplus = T)
    }
  }
}

## var_num
for(ch in 1:mcmc$Nchains){
  for(y in Nums){
    for(g in 1:mcmc$modeGplus[[ch]]){
      plot.all(mcmc = mcmc, thin = figthin, gspec = g,
               what = "var_num", yspec = y, MaxLag = 30,
               ChainsToPlot = ch, restrictedtoGplus = T)
    }
  }
}

## beta_num
nfr <- length(mcmc$lgrpnames$ynum)
for(ch in 1:mcmc$Nchains){
  for(g in 1:mcmc$modeGplus[[ch]]){
    # traceplots
    #pdf(paste0(WD, "Figures/traceplots_beta_num_Gplus_", mcmc$modeGplus[[ch]], "_", g, "_ch_", ch, ".pdf"), 
    #    width = 7, height = 6)
    par(mfrow = c(length(Nums),nfr))
    for (y in Nums){
      for(i in 1:nfr){
        plot.traceplots(mcmc = mcmc, thin = figthin, what = "beta_num", 
                        gspec = g, yspec = y, dimspec = i, labcex = 0.6,
                        ChainsToPlot = ch, restrictedtoGplus = T)
      }
    }
    #dev.off()
  }
}

## beta_num - classes comparison
for(ch in 1:Nchains){
  #pdf(paste0(WD, "Figures/classes_beta_num_Gplus_", mcmc$modeGplus[[ch]], "_ch_", ch, ".pdf"), 
  #    width = 7, height = 6)
  par(mfrow = c(length(Nums),nfr), mar = c(3.5,3.5,0.8,0.8))
  for(y in Nums){
    for(i in 1:nfr){
      plot.classes(mcmc = mcmc, thin = figthin, doECDF = F, ChainsToPlot = ch,
                   what = "beta_num", yspec = y, dimspec = i, parmfrow = F, 
                   restrictedtoGplus = T)
    }
  }
  #dev.off()
}

## beta_num
nfr <- length(mcmc$lgrpnames$ypoi)
for(ch in 1:mcmc$Nchains){
  for(g in 1:mcmc$modeGplus[[ch]]){
    # traceplots
    #pdf(paste0(WD, "Figures/traceplots_beta_poi_Gplus_", mcmc$modeGplus[[ch]], "_", g, "_ch_", ch, ".pdf"), 
    #    width = 7, height = 6)
    par(mfrow = c(length(Pois),nfr))
    for (y in Pois){
      for(i in 1:nfr){
        plot.traceplots(mcmc = mcmc, thin = figthin, what = "beta_poi", 
                        gspec = g, yspec = y, dimspec = i, labcex = 0.6,
                        ChainsToPlot = ch, restrictedtoGplus = T)
      }
    }
    #dev.off()
  }
}

## beta_num - classes comparison
for(ch in 1:Nchains){
  #pdf(paste0(WD, "Figures/classes_beta_poi_Gplus_", mcmc$modeGplus[[ch]], "_ch_", ch, ".pdf"), 
  #    width = 7, height = 6)
  par(mfrow = c(length(Pois),nfr), mar = c(3.5,3.5,0.8,0.8))
  for(y in Pois){
    for(i in 1:nfr){
      plot.classes(mcmc = mcmc, thin = figthin, doECDF = F, ChainsToPlot = ch,
                   what = "beta_poi", yspec = y, dimspec = i, parmfrow = F, 
                   restrictedtoGplus = T)
    }
  }
  #dev.off()
}

## beta_bin
nfr <- length(mcmc$lgrpnames$ybin)
for(ch in 1:mcmc$Nchains){
  for(g in 1:mcmc$modeGplus[[ch]]){
    # traceplots
    #pdf(paste0(WD, "Figures/traceplots_beta_bin_Gplus_", mcmc$modeGplus[[ch]], "_", g, "_ch_", ch, ".pdf"), 
    #    width = 7, height = 6)
    par(mfrow = c(length(Bins),nfr))
    for (y in Bins){
      for(i in 1:nfr){
        plot.traceplots(mcmc = mcmc, thin = figthin, what = "beta_bin", 
                        gspec = g, yspec = y, dimspec = i, labcex = 0.6,
                        ChainsToPlot = ch, restrictedtoGplus = T)
      }
    }
    #dev.off()
  }
}


## beta_bin - classes comparison
for(ch in 1:Nchains){
  #pdf(paste0(WD, "Figures/classes_beta_bin_Gplus_", mcmc$modeGplus[[ch]], "_ch_", ch, ".pdf"), 
  #    width = 7, height = 6)
  par(mfrow = c(length(Bins),nfr), mar = c(3.5,3.5,0.8,0.8))
  for(y in Bins){
    for(i in 1:nfr){
      plot.classes(mcmc = mcmc, thin = figthin, doECDF = F, ChainsToPlot = ch,
                   what = "beta_bin", yspec = y, dimspec = i, parmfrow = F, 
                   restrictedtoGplus = T)
    }
  }
  #dev.off()
}

## beta_ord
nfr <- length(mcmc$lgrpnames$yord)
for(ch in 1:mcmc$Nchains){
  for(g in 1:mcmc$modeGplus[[ch]]){
    # traceplots
    #pdf(paste0(WD, "Figures/traceplots_beta_ord_Gplus_", mcmc$modeGplus[[ch]], "_", g, "_ch_", ch, ".pdf"), 
    #    width = 7, height = 6)
    par(mfrow = c(length(Ords),nfr))
    for (y in Ords){
      for(i in 1:nfr){
        plot.traceplots(mcmc = mcmc, thin = figthin, what = "beta_ord", 
                        gspec = g, yspec = y, dimspec = i, labcex = 0.6,
                        ChainsToPlot = ch, restrictedtoGplus = T)
      }
    }
    #dev.off()
  }
}


## beta_ord - classes comparison
for(ch in 1:Nchains){
  #pdf(paste0(WD, "Figures/classes_beta_ord_Gplus_", mcmc$modeGplus[[ch]], "_ch_", ch, ".pdf"), 
  #    width = 7, height = 6)
  par(mfrow = c(length(Ords),nfr), mar = c(3.5,3.5,0.8,0.8))
  for(y in Ords){
    for(i in 1:nfr){
      plot.classes(mcmc = mcmc, thin = figthin, doECDF = F, ChainsToPlot = ch,
                   what = "beta_ord", yspec = y, dimspec = i, parmfrow = F, 
                   restrictedtoGplus = T)
    }
  }
  #dev.off()
}

## c_ord
names(Kord) <- Ords
par(mfrow = c(length(Ords),max(Kord)-1))
for (y in Ords){
  for(i in 1:(Kord[y]-1)){
    plot.classes(mcmc = mcmc, thin = figthin, doECDF = F, ChainsToPlot = ch,
                 what = "c_ord", yspec = y, dimspec = i, parmfrow = F, 
                 restrictedtoGplus = T)
  }
}


## beta_cat
names(Kcat) <- Cats
nfr <- length(mcmc$lgrpnames$ycat)
for(ch in 1:Nchains){
  for(g in 1:mcmc$modeGplus[[ch]]){
    # traceplots
    #pdf(paste0(WD, "Figures/traceplots_beta_cat_Gplus_", mcmc$modeGplus[[ch]], "_", g, "_ch_", ch, ".pdf"), 
    #    width = 7, height = 6)
    for(y in Cats){
      par(mfrow = c(Kcat[y],nfr))
      for(k in 1:Kcat[y]){
        for(i in 1:nfr){
          plot.traceplots(mcmc = mcmc, thin = figthin, gspec = g, 
                          what = "beta_cat", yspec = y, dimspec = c(i,k), labcex = 0.6,
                          ChainsToPlot = ch, restrictedtoGplus = T)
        }
      }
    }
    #dev.off()
  }
}

## beta_cat - classes comparison
for(ch in 1:Nchains){
  #pdf(paste0(WD, "Figures/classes_beta_cat_Gplus_", mcmc$modeGplus[[ch]], "_ch_", ch, ".pdf"), 
  #    width = 7, height = 6)
  for(y in Cats){
    par(mfrow = c(Kcat[y],nfr), mar = c(3.5,3.5,0.8,0.8))
    for(k in 1:Kcat[y]){
      for(i in 1:nfr){
        plot.classes(mcmc = mcmc, thin = figthin, doECDF = F, ChainsToPlot = ch,
                     what = "beta_cat", yspec = y, dimspec = c(i,k), parmfrow = F, 
                     restrictedtoGplus = T)
      }
    }
  }
  #dev.off()
}

# Saving the image
save.image(file = paste0(RDATA, "tutorial_image03", trial, ".RData"))

# Saving the mcmc list
save(mcmc, file = paste0(RDATA, "tutorial_mcmc03", trial, ".RData"))

