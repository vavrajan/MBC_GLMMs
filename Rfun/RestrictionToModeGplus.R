### Post-processing 
### Finding the mode of Gplus
### Performing k-means clustering to find appropriate permutations
### Keeping only permuted k-specific parameter data

### mcmc contains traditional output of Metropolis_within_Gibbs_MBC_NumBinOrdCat
### mcmc HAS to contain sampled U's
### or at least counts of nUg for each sample


RestrictionToModeGplus <- function(mcmc, cols){
  
  ### 0) The type of saving results has to be "matrix" for the code below
  if(mcmc$howsave=="list"){
    # Then transfer to howsave="matrix" -- create mcmc$all matrix
    mcmc <- FromListtoMatrix_settings(mcmc)
  }
  
  ### 1) Finding the mode
  # Using mcmc$all$Gplus
  if(length(unique(mcmc$all$Gplus)) == 1){
    TAB <- list()
    for(ch in 1:mcmc$Nchains){
      TAB[[ch]] <- mcmc$BM
      names(TAB[[ch]]) <- as.character(mcmc$all$Gplus[1])
    }
  }else{
    TAB <- tapply(factor(mcmc$all$Gplus), factor(mcmc$all$chain), table)
  }
  maxs <- lapply(TAB, which.max)
  modeGplus <- lapply(maxs, function(x){as.numeric(names(x))})
  
  ### 2) Find the posterior draws with Gplus = modeGplus
  inds <- Ms <- list()
  for(ch in 1:mcmc$Nchains){
    submcmc <- mcmc$all[mcmc$all$chain==ch, c("chain", "Gplus")]
    inds[[ch]] <- which(submcmc$Gplus == modeGplus[[ch]])
    Ms[[ch]] <- length(inds[[ch]])
  }
  mcmc$modeGplus <- modeGplus
  mcmc$Ms <- Ms
  print("modeGplus")
  print(mcmc$modeGplus)
  print("Frequency of the modeGlus value")
  print(Ms)
  
  ### 3) Columns of mcmc$all that are cluster-specific parameters
  if(missing(cols)){
    cols <- colnames(mcmc$all)[grep("\\(\\d+\\)", colnames(mcmc$all))]
  }else{
    # use only selected cols for k-means clustering
  }
  print("columns used for k-means")
  print(cols)
  
  ### 4) Creating submatrices of sampled cluster-specific parameters
  clusterspec <- list()
  for(ch in 1:mcmc$Nchains){
    suball <- mcmc$all[mcmc$all$chain==ch, c("chain", "Gplus", paste0("ng[", 1:mcmc$G, "]"), cols)]
    suball <- suball[inds[[ch]], ]
    clusterspec[[ch]] <- matrix(0, nrow = 0, ncol = 1 + length(cols)/mcmc$G*modeGplus[[ch]])
    for(m in 1:Ms[[ch]]){
      whichG <- which(suball[m, paste0("ng[", 1:mcmc$G, "]")] > 0)
      print(whichG)
      expr <- paste(paste0("\\(", whichG, "\\)"), collapse = "|")
      subcols <- grep(expr, colnames(suball))
      clusterspec[[ch]] <- rbind(clusterspec[[ch]], c(inds[[ch]][m], suball[m,subcols]))
    }
    expr <- paste(paste0("\\(", 1:modeGplus[[ch]], "\\)"), collapse = "|")
    subcols <- colnames(suball)[grep(expr, colnames(suball))]
    cnames <- c("m", subcols)
    colnames(clusterspec[[ch]]) <- cnames
  }
  cat("Subset of states of cluster-specific parameters with Gplus=modeGplus 
created for each chain separately")
  
  
  ### 5) Prepare for k-means clustering
  forkmeans <- list()
  for(ch in 1:mcmc$Nchains){
    forkmeans[[ch]] <- matrix(0, nrow = 0, ncol = length(cols)/mcmc$G)
    for(g in 1:modeGplus[[ch]]){
      forkmeans[[ch]] <- rbind(forkmeans[[ch]], clusterspec[[ch]][,grep(paste0("\\(", g, "\\)"), 
                                                                        colnames(clusterspec[[ch]]))])
    }
    cnames <- colnames(clusterspec[[ch]])[grep("\\(1\\)", 
                                               colnames(clusterspec[[ch]]))]
    colnames(forkmeans[[ch]]) <- gsub("\\(1\\)", "", cnames)
  }
  print("Group-specific parameters stacked below each other, prepared for k-means")
  
  ### 6) Perform K-means clustering + checking the permutations
  Mrhos <- indmat <- rowuniques <- list()
  for(ch in 1:mcmc$Nchains){
    # standardization
    forkmeans[[ch]] <- matrix(as.numeric(forkmeans[[ch]]), ncol = dim(forkmeans[[ch]])[2])
    colnames(forkmeans[[ch]]) <- paste0("V", 1:dim(forkmeans[[ch]])[2])
    #class(forkmeans[[ch]])
    #dim(forkmeans[[ch]])
    #summary(forkmeans[[ch]])
    forkmeans[[ch]] <- scale(as.matrix(forkmeans[[ch]]))
    km <- kmeans(forkmeans[[ch]], centers = modeGplus[[ch]])
    # create index matrix
    indmat[[ch]] <- matrix(km$cluster, ncol = modeGplus[[ch]], byrow = F)
    head(indmat[[ch]])
    # number of not permutations = number of rows where is something repeated
    rowuniques[[ch]] <- apply(indmat[[ch]], 1, function(x){length(unique(x))})
    Mrhos[[ch]] <- sum(rowuniques[[ch]] < modeGplus[[ch]])
  }
  mcmc$Mrhos <- Mrhos
  print("K-means performed")
  print("Mrhos")
  print(Mrhos)
  
  ### 7) Update the clusterspec matrix according to suggested permutations
  for(ch in 1:mcmc$Nchains){
    head(clusterspec[[ch]])
    cnames <- colnames(clusterspec[[ch]])[grep("\\(1\\)", 
                                               colnames(clusterspec[[ch]]))]
    for(m in 1:Ms[[ch]]){
      if(rowuniques[[ch]][m] == modeGplus[[ch]]){
        # it is a permutation --> permute accordingly
        auxrow <- unlist(clusterspec[[ch]][m,])
        for(g in 1:modeGplus[[ch]]){
          gnames <- gsub("\\(1\\)", paste0("\\(",g,"\\)"), cnames)
          pnames <- gsub("\\(1\\)", paste0("\\(",indmat[[ch]][m,g],"\\)"), cnames)
          clusterspec[[ch]][m,pnames] <- auxrow[gnames]
        }
      }else{
        # not a permutation --> replace with NA
        clusterspec[[ch]][m,] <- NA
      }
    }
  }
  mcmc$clusterspec <- clusterspec
  print("Each row permuted according to the found permutation")
  
  ### 8) Switch w, U, pUig according to indmat
  wGplus <- ngGplus <- UGplus <- pUigGplus <- list()
  for(ch in 1:mcmc$Nchains){
    ## w
    if(mcmc$whatsave["w"]){
      suball <- mcmc$all[mcmc$all$chain==ch, c(paste0("w[", 1:mcmc$G, "]"),
                                               paste0("ng[",1:mcmc$G, "]"))]
      suball <- suball[inds[[ch]], ]
      wGplus[[ch]] <- matrix(0, nrow = 0, ncol = modeGplus[[ch]])
      ngGplus[[ch]] <- matrix(0, nrow = 0, ncol = modeGplus[[ch]])
      for(m in 1:Ms[[ch]]){
        whichG <- which(suball[m, paste0("ng[", 1:mcmc$G, "]")] > 0)
        wcols <- paste0("w[", whichG, "]")
        ngcols <- paste0("ng[", whichG, "]")
        wGplus[[ch]] <- rbind(wGplus[[ch]], suball[m,wcols])
        ngGplus[[ch]] <- rbind(ngGplus[[ch]], suball[m,ngcols])
      }
      colnames(wGplus[[ch]]) <- paste0("w[", 1:modeGplus[[ch]], "]")
      colnames(ngGplus[[ch]]) <- paste0("ng[", 1:modeGplus[[ch]], "]")
    
      for(m in 1:Ms[[ch]]){
        if(rowuniques[[ch]][m] == modeGplus[[ch]]){
          # it is a permutation --> permute accordingly
          auxrow <- unlist(wGplus[[ch]][m,])
          for(g in 1:modeGplus[[ch]]){
            gname <- paste0("w[",g,"]")
            pname <- paste0("w[",indmat[[ch]][m,g],"]")
            wGplus[[ch]][m,pname] <- auxrow[gname]
          }
          
          # it is a permutation --> permute accordingly
          auxrow <- unlist(ngGplus[[ch]][m,])
          for(g in 1:modeGplus[[ch]]){
            gname <- paste0("ng[",g,"]")
            pname <- paste0("ng[",indmat[[ch]][m,g],"]")
            ngGplus[[ch]][m,pname] <- auxrow[gname]
          }
        }else{
          # not a permutation --> replace with NA
          wGplus[[ch]][m,] <- NA
          ngGplus[[ch]][m,] <- NA
        }
      }
    }
    ## U
    if(mcmc$whatsave["U"]){
      suball <- mcmc$all[mcmc$all$chain==ch, c(paste0("U[", 1:mcmc$settings["U","dims"], "]"),
                                               paste0("ng[",1:mcmc$G, "]"))]
      suball <- suball[inds[[ch]], ]
      UGplus[[ch]] <- matrix(0, nrow = Ms[[ch]], ncol = mcmc$settings["U","dims"])
      for(m in 1:Ms[[ch]]){
        whichG <- which(suball[m, paste0("ng[", 1:mcmc$G, "]")] > 0)
        if(rowuniques[[ch]][m] == modeGplus[[ch]]){
          for(i in 1:mcmc$settings["U","dims"]){
            u <- 1+suball[m,paste0("U[",i,"]")] # +1 since U in 0,1,...,G-1
            whu <- which(whichG == u) # order within non-zero clusters
            UGplus[[ch]][m,i] <- indmat[[ch]][m,whu] # applying permutation
          }
        }else{
          UGplus[[ch]][m,] <- rep(NA, mcmc$settings["U","dims"])
        }
      }
    }
    
    ## pUig
    if(mcmc$whatsave["pUig"]){
      pUigGplus[[ch]] <- matrix(0, nrow = Ms[[ch]], ncol = 0)
      for(i in 1:mcmc$settings["U","dims"]){
        suball <- mcmc$all[mcmc$all$chain==ch, c(paste0("pUig[",i,",", 1:mcmc$G, "]"),
                                                 paste0("ng[",1:mcmc$G, "]"))]
        suball <- suball[inds[[ch]], ]
        pU <- matrix(0, nrow = 0, ncol = modeGplus[[ch]])
        for(m in 1:Ms[[ch]]){
          whichG <- which(suball[m, paste0("ng[", 1:mcmc$G, "]")] > 0)
          pUigcols <- paste0("pUig[",i,",", whichG, "]")
          pU <- rbind(pU, suball[m,pUigcols])
        }
        colnames(pU) <- paste0("pUig[",i,",", 1:modeGplus[[ch]], "]")
        for(m in 1:Ms[[ch]]){
          if(rowuniques[[ch]][m] == modeGplus[[ch]]){
            # it is a permutation --> permute accordingly
            auxrow <- unlist(pU[m,])
            for(g in 1:modeGplus[[ch]]){
              gname <- paste0("pUig[",i,",",g,"]")
              pname <- paste0("pUig[",i,",",indmat[[ch]][m,g],"]")
              pU[m,pname] <- auxrow[gname]
            }
          }else{
            # not a permutation --> replace with NA
            pU[m,] <- rep(NA, modeGplus[[ch]])
          }
        }
        pUigGplus[[ch]] <- cbind(pUigGplus[[ch]], pU)
      }
    }
  }
  cat("All cluster-related parameters are reorganized separately (omitting the empty clusters):
wGplus,
ngGplus,
UGplus,
pUigGplus")
  
  mcmc$wGplus <- wGplus
  mcmc$ngGplus <- ngGplus
  mcmc$UGplus <- UGplus
  mcmc$pUigGplus <- pUigGplus
  mcmc$inds <- inds
  
  return(mcmc)
}
