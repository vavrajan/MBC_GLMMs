PermuteClusterSpec <- function(mcmc, perm, restricted = T){
  # perm contains list of permutations - each for one chain
  # perm[[ch]][1] contains which old cluster should be the new first one (for chain ch)
  # ... and so on
  
  if(restricted){
    ## permute only the restricted parameters according to given list of permutations
    for(ch in 1:mcmc$Nchains){
      # permute w
      if(mcmc$whatsave["pUig"]){
        auxwGplus <- matrix(-1, nrow = dim(mcmc$wGplus[[ch]])[1], ncol = 0)
        for(g in 1:mcmc$modeGplus[[ch]]){
          auxwGplus <- cbind(auxwGplus, mcmc$wGplus[[ch]][,perm[[ch]][g]])
        }
        colnames(auxwGplus) <- colnames(mcmc$wGplus[[ch]])
        mcmc$wGplus[[ch]] <- auxwGplus
      }
      
      # permute ng
      if(mcmc$whatsave["ng"]){
        auxngGplus <- matrix(-1, nrow = dim(mcmc$ngGplus[[ch]])[1], ncol = 0)
        for(g in 1:mcmc$modeGplus[[ch]]){
          auxngGplus <- cbind(auxngGplus, mcmc$ngGplus[[ch]][,perm[[ch]][g]])
        }
        colnames(auxngGplus) <- colnames(mcmc$ngGplus[[ch]])
        mcmc$ngGplus[[ch]] <- auxngGplus
      }
      
      # permute U
      if(mcmc$whatsave["U"]){
        mcmc$UGplus[[ch]] <- apply(mcmc$UGplus[[ch]], c(1,2), function(g){which(perm[[ch]] == g)})
      }
      
      # permute pUig
      if(mcmc$whatsave["pUig"]){
        auxpU <- matrix(-1, nrow = dim(mcmc$pUigGplus[[ch]])[1], ncol = 0)
        for(i in 1:mcmc$settings["U","dims"]){
          auxpUig <- matrix(-1, nrow = dim(mcmc$pUigGplus[[ch]])[1], ncol = 0)
          auxpUi <- mcmc$pUigGplus[[ch]][,paste0("pUig[",i,",",1:mcmc$modeGplus[[ch]], "]")]
          for(g in 1:mcmc$modeGplus[[ch]]){
            auxpUig <- cbind(auxpUig, auxpUi[,perm[[ch]][g]])
          }
          colnames(auxpUig) <- colnames(auxpUi)
          auxpU <- cbind(auxpU, auxpUig)
        }
        mcmc$pUigGplus[[ch]] <- auxpU
      }
      
      # permute clusterspec
      aux <- mcmc$clusterspec[[ch]]
      cnames <- colnames(aux)[grep("\\(1\\)", colnames(aux))]
      for(g in 1:mcmc$modeGplus[[ch]]){
        pnames <- gsub("\\(1\\)", paste0("\\(",perm[[ch]][g],"\\)"), cnames)
        gnames <- gsub("\\(1\\)", paste0("\\(",g,"\\)"), cnames)
        aux[,gnames] <- mcmc$clusterspec[[ch]][,pnames]
      } 
      mcmc$clusterspec[[ch]] <- aux
      
    } # end of for ch in 1:Nchains
  }else{
    ## permute the raw sampled g-specific data
    
    if(mcmc$howsave == "matrix"){
      # data are save in matrix mcmc$all
      aux <- mcmc$all
      cnames <- colnames(aux)[grep("\\(1\\)", colnames(aux))]
      for(g in 1:mcmc$G){
        pnames <- gsub("\\(1\\)", paste0("\\(",perm[[ch]][g],"\\)"), cnames)
        gnames <- gsub("\\(1\\)", paste0("\\(",g,"\\)"), cnames)
        aux[,gnames] <- mcmc$all[,pnames]
      } 
      mcmc$all <- aux
    }
    
    if(mcmc$howsave == "list"){
      # data are saved in lists
      for(ch in 1:mcmc$Nchains){
        for(p in rownames(mcmc$settings)[as.logical(mcmc$settings$isspec)]){
          aux <- mcmc[[ch]][[p]]
          for(g in 1:mcmc$G){
            aux[[g]] <- mcmc[[ch]][[p]][[perm[[ch]][g]]]
          } 
          mcmc[[ch]][[p]] <- aux
        }
      }
    }
  }
    
  return(mcmc)
}
