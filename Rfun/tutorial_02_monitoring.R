###-----------------------------------------------------------------------------
###     Monitoring of the sampled chains
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

### Source plotting functions
source(paste0(RROOT, "plotting_functions.R"))


# If you wish to inspect the burn-in period
load(file = paste0(RDATA, "tutorial_image01_burnin.RData"))
trial <- "_burnin"

# If you wish to inspect the final chains
load(file = paste0(RDATA, "tutorial_image01.RData"))
trial <- ""

### Thinning for the following plots
figthin = max(round(M/1000), 1)
figthin = 1
# chtplt = chains to plot
chtplt <- 1:Nchains


### Traceplots, ..., classes comparison
# Skip the cairo_pdf commands if you do not wish to save it
# but to see the plots inside RStudio

## Gplus
cairo_pdf(paste0(FIG, "traceplots_Gplus", trial, "_Gmax_", Gmax, "_nch_", Nchains, "_M_", M, "_B_", B, ".pdf"),
          width = 6, height = 5)
{
  par(mfrow = c(1,1), mar = c(4,4,1,1))
  plot.traceplots(mcmc = mcmc, thin = figthin, B = 0, #M = 5,
                  what = "Gplus", labcex = 1, 
                  whatLAB = "Gplus", ChainsToPlot = chtplt)
}
dev.off()

## e0
cairo_pdf(paste0(FIG, "traceplots_e0", trial, "_Gmax_", Gmax, "_nch_", Nchains, "_M_", M, "_B_", B, ".pdf"),
          width = 6, height = 5)
{
  par(mfrow = c(1,1), mar = c(4,4,1,1))
  plot.traceplots(mcmc = mcmc, thin = figthin, #B = 9, M = 10, 
                  what = "e0", labcex = 1, 
                  whatLAB = "e0", ChainsToPlot = chtplt)
}
dev.off()

## ng
cairo_pdf(paste0(FIG, "traceplots_ng", trial, "_Gmax_", Gmax, "_nch_", Nchains, "_M_", M, "_B_", B, ".pdf"),
          width = 9, height = 6)
{
  par(mfrow = c(2,5), mar = c(4,4,1,1))
  for(g in 1:Gmax){
    plot.traceplots(mcmc = mcmc, thin = figthin, dimspec = g, #B = 0, M = 10,
                    what = "ng", labcex = 0.6, 
                    whatLAB = paste0("nUg[",g,"]"), ChainsToPlot = chtplt)
  }
}
dev.off()

## w
cairo_pdf(paste0(FIG, "traceplots_w", trial, "_Gmax_", Gmax, "_nch_", Nchains, "_M_", M, "_B_", B, ".pdf"),
          width = 9, height = 6)
{
  par(mfrow = c(2,5), mar = c(4,4,1,1))
  for(g in 1:Gmax){
    plot.traceplots(mcmc = mcmc, thin = figthin, #B = 0, M = 100,
                    what = "w", dimspec = g, labcex = 0.6, 
                    whatLAB = paste0("w[",g,"]"), ChainsToPlot = chtplt)
  }
}
dev.off()

## sd Sigma
cairo_pdf(paste0(FIG, "traceplots_sdSigma", trial, "_Gmax_", Gmax, "_nch_", Nchains, "_M_", M, "_B_", B, ".pdf"),
          width = 8, height = 6)
{
  par(mfrow = c(2,3))
  for(j in 1:totnran){
    plot.traceplots(mcmc = mcmc, thin = figthin,
                    what = "sdSigma", dimspec = j, 
                    labcex = 0.6, ChainsToPlot = chtplt)
    abline(h = sdSigmas[j], col = "black", lty = 2)
    
  }
}
dev.off()

## corSigma
cairo_pdf(paste0(FIG, "traceplots_corSigma", trial, "_Gmax_", Gmax, "_nch_", Nchains, "_M_", M, "_B_", B, ".pdf"),
          width = 9, height = 8)
{
  par(mfrow = c(totnran-1, totnran-1), mar = c(4,4,1,1))
  for(i in 1:(totnran-1)){
    for(j in 2:(totnran)){
      if(i < j){
        plot.traceplots(mcmc = mcmc, thin = figthin,
                        what = "corSigma", dimspec = c(i,j),
                        labcex = 0.6, ChainsToPlot = chtplt)
        abline(h = corr[i,j], col = "black", lty = 2)
      }else{
        plot(x=c(0,1), y=c(0,1), type = "n", xlab = "", ylab = "",
             xaxt = "n", yaxt = "n", bty = "n")
      }
    }
  }
}
dev.off()

## Sigma
cairo_pdf(paste0(FIG, "traceplots_Sigma", trial, "_Gmax_", Gmax, "_nch_", Nchains, "_M_", M, "_B_", B, ".pdf"),
          width = 9, height = 8)
{
  par(mfrow = c(totnran,totnran), mar = c(4,4,1,1))
  for (i in 1:totnran){
    for(j in (1:totnran)){
      if(i <= j){
        plot.traceplots(mcmc = mcmc, thin = figthin,
                        what = "Sigma", dimspec = c(i,j),
                        labcex = ifelse(totnran ==3, 0.6, 0.5), ChainsToPlot = chtplt)
        abline(h = Sigma[i,j], col = "black", lty = 2)
      }else{
        pomdata <- FindAllData(mcmc=mcmc, thin = figthin, B = 0,
                               what = "Sigma", dimspec = c(j,i),
                               ChainsToPlot = 1:mcmc$Nchains)
        plot(x=c(0,1), y=c(0,1), type = "n", xlab = "", ylab = "",
             xaxt = "n", yaxt = "n", bty = "n")
        text(0.5,0.5, labels = format(mean(pomdata$y), digits = 2, nsmall = 2))
      }
    }
  }
}
dev.off()

## InvSigma
InvSigma <- solve(Sigma)
cairo_pdf(paste0(FIG, "traceplots_InvSigma", trial, "_Gmax_", Gmax, "_nch_", Nchains, "_M_", M, "_B_", B, ".pdf"),
          width = 9, height = 8)
{
  par(mfrow = c(totnran,totnran), mar = c(4,4,1,1))
  for (i in 1:totnran){
    for(j in (1:totnran)){
      if(i <= j){
        plot.traceplots(mcmc = mcmc, thin = figthin,
                        what = "InvSigma", dimspec = c(i,j),
                        labcex = ifelse(totnran ==3, 0.6, 0.5), ChainsToPlot = chtplt)
        abline(h = InvSigma[i,j], col = "black", lty = 2)
      }else{
        pomdata <- FindAllData(mcmc=mcmc, thin = figthin, B = 0,
                               what = "InvSigma", dimspec = c(j,i),
                               ChainsToPlot = 1:mcmc$Nchains)
        plot(x=c(0,1), y=c(0,1), type = "n", xlab = "", ylab = "",
             xaxt = "n", yaxt = "n", bty = "n")
        text(0.5,0.5, labels = format(mean(pomdata$y), digits = 2, nsmall = 2))
      }
    }
  }
}
dev.off()

## beta_num_fix
# in case some beta coefficients are common to all clusters
nfr <- length(mcmc$lfixnames[[Nums[1]]])
if(nfr > 0){
  cairo_pdf(paste0(FIG, "traceplots_beta_num_fix", trial, "_Gmax_", Gmax, "_nch_", Nchains, "_M_", M, "_B_", B, ".pdf"),
            width = 6, height = 6)
  {
    par(mfrow = c(length(Nums),nfr), mar = c(4,4,1,1))
    
    for(y in Nums){
      for(i in 1:nfr){
        plot.traceplots(mcmc = mcmc, thin = figthin, #B = 375, M = 10,
                        what = "beta_num_fix",
                        yspec = y, dimspec = i, labcex = 0.6,
                        ChainsToPlot = 1:mcmc$Nchains)
        # true values (needs to be updated if formula changed)
        abline(h = beta_num[[1]]["cat"], col = "black", lty = 2)
        
      }
    }
  }
  dev.off()
}


## beta_num
for(ch in 1:Nchains){
  # only the nonempty clusters (decided upon the last state only)
  TAB <- table(mcmc$last[[ch]]$U+1)
  nonempty <- as.numeric(names(TAB)[which(TAB>0)])
  
  cairo_pdf(paste0(FIG, "traceplots_beta_num", trial, "_ch_", ch, "_Gmax_", Gmax, "_nch_", Nchains, "_M_", M, "_B_", B, ".pdf"),
            width = 6, height = 9)
  {
    nfr <- length(mcmc$lgrpnames[[Nums[1]]])
    par(mfrow = c(length(nonempty),nfr), mar = c(4,4,1,1))
    
    for(g in nonempty){
      for(y in Nums){
        for(i in 1:nfr){
          plot.traceplots(mcmc = mcmc, thin = figthin, #B = 375, M = 10,
                          what = "beta_num",
                          gspec = g, yspec = y, dimspec = i, labcex = 0.6,
                          ChainsToPlot = ch)
          # true values (needs to be updated if formula changed)
          for(gg in 1:G){
            abline(h = beta_num[[gg]][i], col = "black", lty = 2)
          }
        }
      }
    }
  }
  dev.off()
}

# comparison between different groups by kernel density estimates
for(ch in 1:Nchains){
  TAB <- table(mcmc$last[[ch]]$U+1)
  nonempty <- as.numeric(names(TAB)[which(TAB>0)])
  
  cairo_pdf(paste0(FIG, "densities_beta_num", trial, "_ch_", ch, "_Gmax_", Gmax, "_nch_", Nchains, "_M_", M, "_B_", B, ".pdf"),
            width = 7, height = 6)
  {
    nfr <- length(mcmc$lgrpnames[[Nums[1]]])
    par(mfrow = c(length(Nums),nfr), mar = c(4,4,1,1))
    
    for(y in Nums){
      for(i in 1:nfr){
        plot.classes(mcmc = mcmc, thin = figthin, B = 0,
                     doECDF = F, ChainsToPlot = ch,
                     whichg = nonempty,
                     what = "beta_num", yspec = y, dimspec = i, parmfrow = F)
        for(gg in 1:G){
          abline(v = beta_num[[gg]][i], col = "black", lty = 2)
        }
      }
    }
  }
  dev.off()
}

## tau_num
for(ch in 1:Nchains){
  # again only parameters corresponding to non-empty clusters
  TAB <- table(mcmc$last[[ch]]$U+1)
  nonempty <- as.numeric(names(TAB)[which(TAB>0)])
  
  cairo_pdf(paste0(FIG, "traceplots_tau_sd_num", trial, "_ch_", ch, "_Gmax_", Gmax, "_nch_", Nchains, "_M_", M, "_B_", B, ".pdf"),
            width = 6, height = 9)
  {
    par(mfrow = c(length(nonempty),2), mar = c(4,4,1,1))
    
    for(g in nonempty){
      for(y in Nums){
        # precision (tau) on the left
        plot.traceplots(mcmc = mcmc, thin = figthin, #B = 375, M = 10,
                        what = "tau_num",
                        gspec = g, yspec = y, labcex = 0.6,
                        ChainsToPlot = ch)
        # true values
        for(gg in 1:G){
          abline(h = 1/(sd_num[[G]][gg])^2, col = "black", lty = 2)
        }
        
        # standard deviation on the right
        plot.traceplots(mcmc = mcmc, thin = figthin, #B = 375, M = 10,
                        what = "sd_num",
                        gspec = g, yspec = y, labcex = 0.6,
                        ChainsToPlot = ch)
        # true values
        for(gg in 1:G){
          abline(h = sd_num[[G]][gg], col = "black", lty = 2)
        }
      }
    }
  }
  dev.off()
}

for(ch in 1:Nchains){
  TAB <- table(mcmc$last[[ch]]$U+1)
  nonempty <- as.numeric(names(TAB)[which(TAB>0)])
  
  cairo_pdf(paste0(FIG, "densities_tau_sd_num", trial, "_ch_", ch, "_Gmax_", Gmax, "_nch_", Nchains, "_M_", M, "_B_", B, ".pdf"),
            width = 7, height = 6)
  {
    par(mfrow = c(length(Nums),2), mar = c(4,4,1,1))
    
    for(y in Nums){
      plot.classes(mcmc = mcmc, thin = figthin, B = 0,
                   doECDF = F, ChainsToPlot = ch,
                   whichg = nonempty,
                   what = "tau_num", yspec = y, parmfrow = F)
      for(gg in 1:G){
        abline(v = 1/(sd_num[[G]][gg])^2, col = "black", lty = 2)
      }
      
      plot.classes(mcmc = mcmc, thin = figthin, B = 0,
                   doECDF = F, ChainsToPlot = ch,
                   whichg = nonempty,
                   what = "sd_num", yspec = y, parmfrow = F)
      for(gg in 1:G){
        abline(v = sd_num[[G]][gg], col = "black", lty = 2)
      }
    }
  }
  dev.off()
}

## beta_poi_fix
nfr <- length(mcmc$lfixnames[[Pois[1]]])
if(nfr > 0){
  cairo_pdf(paste0(FIG, "traceplots_beta_poi_fix", trial, "_Gmax_", Gmax, "_nch_", Nchains, "_M_", M, "_B_", B, ".pdf"),
            width = 6, height = 6)
  {
    par(mfrow = c(length(Pois),nfr), mar = c(4,4,1,1))
    
    for(y in Pois){
      for(i in 1:nfr){
        plot.traceplots(mcmc = mcmc, thin = figthin, #B = 375, M = 10,
                        what = "beta_poi_fix",
                        yspec = y, dimspec = i, labcex = 0.6,
                        ChainsToPlot = 1:mcmc$Nchains)
        
        # true values (needs to be updated if formula changed)
        abline(h = beta_poi[[1]]["cat"], col = "black", lty = 2)
        
      }
    }
  }
  dev.off()
}

## beta_poi
for(ch in 1:Nchains){
  TAB <- table(mcmc$last[[ch]]$U+1)
  nonempty <- as.numeric(names(TAB)[which(TAB>0)])
  
  cairo_pdf(paste0(FIG, "traceplots_beta_poi", trial, "_ch_", ch, "_Gmax_", Gmax, "_nch_", Nchains, "_M_", M, "_B_", B, ".pdf"),
            width = 6, height = 9)
  {
    nfr <- length(mcmc$lgrpnames[[Pois[1]]])
    par(mfrow = c(length(nonempty),nfr), mar = c(4,4,1,1))
    
    for(g in nonempty){
      for(y in Pois){
        for(i in 1:nfr){
          plot.traceplots(mcmc = mcmc, thin = figthin, #B = 375, M = 10,
                          what = "beta_poi",
                          gspec = g, yspec = y, dimspec = i, labcex = 0.6,
                          ChainsToPlot = ch)
          for(gg in 1:G){
            abline(h = beta_poi[[gg]][i], col = "black", lty = 2)
          }
        }
      }
    }
  }
  dev.off()
}

for(ch in 1:Nchains){
  TAB <- table(mcmc$last[[ch]]$U+1)
  nonempty <- as.numeric(names(TAB)[which(TAB>0)])
  
  cairo_pdf(paste0(FIG, "densities_beta_poi", trial, "_ch_", ch, "_Gmax_", Gmax, "_nch_", Nchains, "_M_", M, "_B_", B, ".pdf"),
            width = 7, height = 6)
  {
    nfr <- length(mcmc$lgrpnames[[Pois[1]]])
    par(mfrow = c(length(Bins),nfr), mar = c(4,4,1,1))
    
    for(y in Pois){
      for(i in 1:nfr){
        plot.classes(mcmc = mcmc, thin = figthin, B = 0,
                     doECDF = F, ChainsToPlot = ch,
                     whichg = nonempty,
                     what = "beta_poi", yspec = y, dimspec = i, parmfrow = F)
        for(gg in 1:G){
          abline(v = beta_poi[[gg]][i], col = "black", lty = 2)
        }
      }
    }
  }
  dev.off()
}

## beta_bin_fix
nfr <- length(mcmc$lfixnames[[Bins[1]]])
if(nfr > 0){
  cairo_pdf(paste0(FIG, "traceplots_beta_bin_fix", trial, "_Gmax_", Gmax, "_nch_", Nchains, "_M_", M, "_B_", B, ".pdf"),
            width = 6, height = 6)
  {
    par(mfrow = c(length(Bins),nfr), mar = c(4,4,1,1))
    
    for(y in Bins){
      for(i in 1:nfr){
        plot.traceplots(mcmc = mcmc, thin = figthin, #B = 375, M = 10,
                        what = "beta_bin_fix",
                        yspec = y, dimspec = i, labcex = 0.6,
                        ChainsToPlot = 1:mcmc$Nchains)
        
        # true values (needs to be updated if formula changed)
        # here the x covariate was fixed
        abline(h = beta_bin[[1]]["time1"], col = "black", lty = 2)
        
      }
    }
  }
  dev.off()
}

## beta_bin
for(ch in 1:Nchains){
  TAB <- table(mcmc$last[[ch]]$U+1)
  nonempty <- as.numeric(names(TAB)[which(TAB>0)])
  
  cairo_pdf(paste0(FIG, "traceplots_beta_bin", trial, "_ch_", ch, "_Gmax_", Gmax, "_nch_", Nchains, "_M_", M, "_B_", B, ".pdf"),
            width = 6, height = 9)
  {
    nfr <- length(mcmc$lgrpnames[[Bins[1]]])
    par(mfrow = c(length(nonempty),nfr), mar = c(4,4,1,1))
    
    for(g in nonempty){
      for(y in Bins){
        for(i in 1:nfr){
          plot.traceplots(mcmc = mcmc, thin = figthin, #B = 375, M = 10,
                          what = "beta_bin",
                          gspec = g, yspec = y, dimspec = i, labcex = 0.6,
                          ChainsToPlot = ch)
          for(gg in 1:G){
            abline(h = beta_bin[[gg]][i], col = "black", lty = 2)
          }
        }
      }
    }
  }
  dev.off()
}

for(ch in 1:Nchains){
  TAB <- table(mcmc$last[[ch]]$U+1)
  nonempty <- as.numeric(names(TAB)[which(TAB>0)])
  
  cairo_pdf(paste0(FIG, "densities_beta_bin", trial, "_ch_", ch, "_Gmax_", Gmax, "_nch_", Nchains, "_M_", M, "_B_", B, ".pdf"),
            width = 7, height = 6)
  {
    nfr <- length(mcmc$lgrpnames[[Bins[1]]])
    par(mfrow = c(length(Bins),nfr), mar = c(4,4,1,1))
    
    for(y in Bins){
      for(i in 1:nfr){
        plot.classes(mcmc = mcmc, thin = figthin, B = 0,
                     doECDF = F, ChainsToPlot = ch,
                     whichg = nonempty,
                     what = "beta_bin", yspec = y, dimspec = i, parmfrow = F)
        for(gg in 1:G){
          abline(v = beta_bin[[gg]][i], col = "black", lty = 2)
        }
      }
    }
  }
  dev.off()
}

## beta_ord_fix
nfr <- length(mcmc$lfixnames[[Ords[1]]])
if(nfr > 0){
  cairo_pdf(paste0(FIG, "traceplots_beta_ord_fix", trial, "_Gmax_", Gmax, "_nch_", Nchains, "_M_", M, "_B_", B, ".pdf"),
            width = 6, height = 6)
  {
    par(mfrow = c(length(Ords),nfr), mar = c(4,4,1,1))
    
    for(y in Ords){
      for(i in 1:nfr){
        plot.traceplots(mcmc = mcmc, thin = figthin, #B = 375, M = 10,
                        what = "beta_ord_fix",
                        yspec = y, dimspec = i, labcex = 0.6,
                        ChainsToPlot = 1:mcmc$Nchains)
        
        # true values (needs to be updated if formula changed)
        abline(h = beta_ord[[1]]["cat"], col = "black", lty = 2)
        
      }
    }
  }
  dev.off()
}

## beta_ord
for(ch in 1:Nchains){
  TAB <- table(mcmc$last[[ch]]$U+1)
  nonempty <- as.numeric(names(TAB)[which(TAB>0)])
  
  nfr <- length(mcmc$lgrpnames[[Ords[1]]])
  if(nfr > 0){
    cairo_pdf(paste0(FIG, "traceplots_beta_ord", trial, "_ch_", ch, "_Gmax_", Gmax, "_nch_", Nchains, "_M_", M, "_B_", B, ".pdf"),
              width = 6, height = 9)
    {
      par(mfrow = c(length(nonempty),nfr), mar = c(4,4,1,1))
      
      for(g in nonempty){
        for(y in Ords){
          for(i in 1:nfr){
            plot.traceplots(mcmc = mcmc, thin = figthin, #B = 375, M = 10,
                            what = "beta_ord",
                            gspec = g, yspec = y, dimspec = i, labcex = 0.6,
                            ChainsToPlot = ch)
            # true values (careful, here the intercept is stand alone as "c_ord" parameter)
            for(gg in 1:G){
              abline(h = beta_ord[[gg]][i+1], col = "black", lty = 2)
            }
          }
        }
      }
    }
    dev.off()
  }
}

for(ch in 1:Nchains){
  TAB <- table(mcmc$last[[ch]]$U+1)
  nonempty <- as.numeric(names(TAB)[which(TAB>0)])
  
  nfr <- length(mcmc$lgrpnames[[Ords[1]]])
  if(nfr > 0){
    cairo_pdf(paste0(FIG, "densities_beta_ord", trial, "_ch_", ch, "_Gmax_", Gmax, "_nch_", Nchains, "_M_", M, "_B_", B, ".pdf"),
              width = 7, height = 6)
    {
      par(mfrow = c(length(Ords),nfr), mar = c(4,4,1,1))
      
      for(y in Ords){
        for(i in 1:nfr){
          plot.classes(mcmc = mcmc, thin = figthin, B = 0,
                       doECDF = F, ChainsToPlot = ch,
                       whichg = nonempty,
                       what = "beta_ord", yspec = y, dimspec = i, parmfrow = F)
          # true values (careful, here the intercept is stand alone as "c_ord" parameter)
          for(gg in 1:G){
            abline(v = beta_ord[[gg]][i+1], col = "black", lty = 2)
          }
        }
      }
    }
    dev.off()
  }
  
}

## c_ord
Kord <- numeric(length(Ords))
names(Kord) <- Ords
for(oo in Ords){
  Kord[oo] <- nlevels(as.factor(Y[,oo]))-1
}
for(ch in 1:Nchains){
  TAB <- table(mcmc$last[[ch]]$U+1)
  nonempty <- as.numeric(names(TAB)[which(TAB>0)])
  
  cairo_pdf(paste0(FIG, "traceplots_c_ord", trial, "_ch_", ch, "_Gmax_", Gmax, "_nch_", Nchains, "_M_", M, "_B_", B, ".pdf"),
            width = 9, height = 9)
  {
    par(mfrow = c(length(nonempty),max(Kord)), mar = c(4,4,1,1))
    
    for(g in nonempty){
      for(y in Ords){
        for(i in 1:Kord[y]){
          plot.traceplots(mcmc = mcmc, thin = figthin, #B = 375, M = 10,
                          what = "c_ord",
                          gspec = g, yspec = y, dimspec = i, labcex = 0.6,
                          ChainsToPlot = ch)
          # true values
          for(gg in 1:G){
            abline(h = c_ord[[gg]][i+1], col = "black", lty = 2)
          }
        }
      }
    }
  }
  dev.off()
}

for(ch in 1:Nchains){
  TAB <- table(mcmc$last[[ch]]$U+1)
  nonempty <- as.numeric(names(TAB)[which(TAB>0)])
  
  cairo_pdf(paste0(FIG, "densities_c_ord", trial, "_ch_", ch, "_Gmax_", Gmax, "_nch_", Nchains, "_M_", M, "_B_", B, ".pdf"),
            width = 9, height = 6)
  {
    par(mfrow = c(length(Ords),max(Kord)), mar = c(4,4,1,1))
    
    for(y in Ords){
      for(i in 1:Kord[y]){
        plot.classes(mcmc = mcmc, thin = figthin, B = 0,
                     doECDF = F, ChainsToPlot = ch,
                     whichg = nonempty,
                     what = "c_ord", yspec = y, dimspec = i, parmfrow = F)
        # true values
        for(gg in 1:G){
          abline(v = c_ord[[gg]][i+1], col = "black", lty = 2)
        }
      }
    }
  }
  dev.off()
}

## pi_ord
for(ch in 1:Nchains){
  TAB <- table(mcmc$last[[ch]]$U+1)
  nonempty <- as.numeric(names(TAB)[which(TAB>0)])
  
  cairo_pdf(paste0(FIG, "traceplots_pi_ord", trial, "_ch_", ch, "_Gmax_", Gmax, "_nch_", Nchains, "_M_", M, "_B_", B, ".pdf"),
            width = 9, height = 9)
  {
    par(mfrow = c(length(nonempty),max(Kord+1)), mar = c(4,4,1,1))
    
    for(g in nonempty){
      for(y in Ords){
        for(i in 1:(Kord[y]+1)){
          plot.traceplots(mcmc = mcmc, thin = figthin, #B = 375, M = 10,
                          what = "pi_ord",
                          gspec = g, yspec = y, dimspec = i, labcex = 0.6,
                          ChainsToPlot = ch)
        }
      }
    }
  }
  dev.off()
}

for(ch in 1:Nchains){
  TAB <- table(mcmc$last[[ch]]$U+1)
  nonempty <- as.numeric(names(TAB)[which(TAB>0)])
  
  cairo_pdf(paste0(FIG, "densities_pi_ord", trial, "_ch_", ch, "_Gmax_", Gmax, "_nch_", Nchains, "_M_", M, "_B_", B, ".pdf"),
            width = 9, height = 6)
  {
    par(mfrow = c(length(Ords),max(Kord+1)), mar = c(4,4,1,1))
    
    for(y in Ords){
      for(i in 1:(Kord[y]+1)){
        plot.classes(mcmc = mcmc, thin = figthin, B = 0,
                     doECDF = F, ChainsToPlot = ch,
                     whichg = nonempty,
                     what = "pi_ord", yspec = y, dimspec = i, parmfrow = F)
      }
    }
  }
  dev.off()
}

## beta_cat_fix
Kcat <- numeric(length(Cats))
names(Kcat) <- Cats
for(cc in Cats){
  Kcat[cc] <- nlevels(as.factor(Y[,cc]))-1
}
nfr <- length(mcmc$lfixnames[[Cats[1]]])
if(nfr > 0){
  cairo_pdf(paste0(FIG, "traceplots_beta_cat_fix", trial, "_Gmax_", Gmax, "_nch_", Nchains, "_M_", M, "_B_", B, ".pdf"),
            width = 6, height = 6)
  {
    par(mfrow = c(nfr, max(Kcat)), mar = c(4,4,1,1))
    
    for(y in Cats){
      for(i in 1:nfr){
        for(k in 1:Kcat[y]){
          plot.traceplots(mcmc = mcmc, thin = figthin, #B = 375, M = 10,
                          what = "beta_cat_fix",
                          yspec = y, dimspec = c(i,k), labcex = 0.6,
                          ChainsToPlot = 1:mcmc$Nchains)
          
          # true values (needs to be updated if formula changed)
          abline(h = beta_cat[[1]][k,"cat"], col = "black", lty = 2)
        }
      }
    }
  }
  dev.off()
}

## beta_cat
nfr <- length(mcmc$lgrpnames[[Cats[1]]])
for(ch in 1:Nchains){
  TAB <- table(mcmc$last[[ch]]$U+1)
  nonempty <- as.numeric(names(TAB)[which(TAB>0)])
  
  cairo_pdf(paste0(FIG, "traceplots_beta_cat", trial, "_ch_", ch, "_Gmax_", Gmax, "_nch_", Nchains, "_M_", M, "_B_", B, ".pdf"),
            width = 9, height = 9)
  {
    nfr <- length(mcmc$lgrpnames[[Cats[1]]])
    par(mfrow = c(length(nonempty),nfr*max(Kcat)), mar = c(4,4,1,1))
    
    for(g in nonempty){
      for(y in Cats){
        for(i in 1:nfr){
          for(k in 1:(Kcat[y])){
            plot.traceplots(mcmc = mcmc, thin = figthin, #B = 375, M = 10,
                            what = "beta_cat",
                            gspec = g, yspec = y, dimspec = c(i,k), labcex = 0.6,
                            ChainsToPlot = ch)
            for(gg in 1:G){
              abline(h = beta_cat[[gg]][k,i], col = "black", lty = 2)
            }
          }
        }
      }
    }
  }
  dev.off()
}

for(ch in 1:Nchains){
  TAB <- table(mcmc$last[[ch]]$U+1)
  nonempty <- as.numeric(names(TAB)[which(TAB>0)])
  
  cairo_pdf(paste0(FIG, "densities_beta_cat", trial, "_ch_", ch, "_Gmax_", Gmax, "_nch_", Nchains, "_M_", M, "_B_", B, ".pdf"),
            width = 9, height = 9)
  {
    nfr <- length(mcmc$lgrpnames[[Cats[1]]])
    par(mfrow = c(max(Kcat), nfr), mar = c(4,4,1,1))
    
    for(y in Cats){
      for(k in 1:(Kcat[y])){
        for(i in 1:nfr){
          plot.classes(mcmc = mcmc, thin = figthin, B = 0,
                       doECDF = F, ChainsToPlot = ch,
                       whichg = nonempty,
                       what = "beta_cat", yspec = y, dimspec = c(i,k), parmfrow = F)
          for(gg in 1:G){
            abline(v = beta_cat[[gg]][k,i], col = "black", lty = 2)
          }
        }
      }
    }
  }
  dev.off()
}

## loglik
cairo_pdf(paste0(FIG, "traceplots_loglik", trial, "_Gmax_", Gmax, "_nch_", Nchains, "_M_", M, "_B_", B, ".pdf"),
          width = 9, height = 9)
{
  par(mar = c(4,4,1,1))
  plot.all(mcmc = mcmc, thin = figthin, what = "loglik", MaxLag = 30,
           ChainsToPlot = chtplt)
}
dev.off()

## NA values
yindices <- rep(Y[,Id], sum(mcmc$nY))
ytypes <- rep(Ys, each = dim(Y)[1])
ynames <- paste0("naY_",
                 ytypes,
                 "[", rep(1:dim(Y)[1], sum(mcmc$nY)), "]")

nayindices <- yindices[mcmc$isYna]
naytypes <- ytypes[mcmc$isYna]
naynames <- ynames[mcmc$isYna]

# true values
naij
natrue

cairo_pdf(paste0(FIG, "traceplots_naY_", trial, "_Gmax_", Gmax, "_nch_", Nchains, "_M_", M, "_B_", B, ".pdf"),
          width = 9, height = 6)
{
  par(mfrow = c(2,5))
  for(j in 1:length(nayindices)){
    i = nayindices[j]
    if(is.element(naytypes[j], c(Nums, Pois))){
      plot.traceplots(mcmc = mcmc, thin = figthin,
                      what = "naY", whatLAB = naynames[j],
                      dimspec = j, labcex = 0.6,
                      ChainsToPlot = chtplt)
      # true missing value
      abline(h = natrue[j], col = "black", lty = 2)
      
    }else{
      par(mar = c(4,4,2.5,2.5))
      plot(factor(mcmc$all[,paste0("naY[",j,"]")]) ~ factor(mcmc$all$chain),
           main = paste0("true value = ", natrue[j]),
           xlab = "chain", ylab = naynames[j])
    }
  }
}
dev.off()

