FromListtoMatrix_settings <- function(mcmc){
  
  params <- unique(c(names(mcmc$whatsave$actual),
                     names(mcmc$whatsave$latent),
                     names(mcmc$whatsave$auxiliary),
                     names(mcmc$whatsave$sparse),
                     names(mcmc$calc)))
  settings <- mcmc$settings
  yspecd1 <- mcmc$yspecd1
  yspecd2 <- mcmc$yspecd2
  Nums <- mcmc$Nums
  Pois <- mcmc$Pois
  Bins <- mcmc$Bins
  Ords <- mcmc$Ords
  Cats <- mcmc$Cats
  cmcmc <- list()
  mcmc$all <- data.frame()
  
  # transfer from list to C
  for(ch in 1:mcmc$Nchains){
    for(p in params){
      if(settings[p, "save"]){
        if(settings[p,"ydepd1"] & settings[p,"ydepd2"]){
          cmcmc[[p]] <- FromListtoC_settings(chain = mcmc[[ch]],
                                             p = p,
                                             settings = settings,
                                             yspecd1 = yspecd1[[p]],
                                             yspecd2 = yspecd2[[p]],
                                             Nums = Nums, Pois = Pois, Bins = Bins, Ords = Ords, Cats = Cats)
        }
        if(settings[p,"ydepd1"] & !settings[p,"ydepd2"]){
          cmcmc[[p]] <- FromListtoC_settings(chain = mcmc[[ch]],
                                             p = p,
                                             settings = settings,
                                             yspecd1 = yspecd1[[p]],
                                             Nums = Nums, Pois = Pois, Bins = Bins, Ords = Ords, Cats = Cats)
        }
        if(!settings[p,"ydepd1"] & settings[p,"ydepd2"]){
          cmcmc[[p]] <- FromListtoC_settings(chain = mcmc[[ch]],
                                             p = p,
                                             settings = settings,
                                             yspecd2 = yspecd2[[p]],
                                             Nums = Nums, Pois = Pois, Bins = Bins, Ords = Ords, Cats = Cats)
        }
        if(!settings[p,"ydepd1"] & !settings[p,"ydepd2"]){
          cmcmc[[p]] <- FromListtoC_settings(chain = mcmc[[ch]],
                                             p = p,
                                             settings = settings,
                                             Nums = Nums, Pois = Pois, Bins = Bins, Ords = Ords, Cats = Cats)
        }
      }
    }
    
    AllData <- matrix(ch, nrow = B+M, ncol = 1)
    colnames(AllData) <- "chain"
    for(p in params){
      if(settings[p, "save"]){
        # transfer from C to matrix
        if(settings[p,"ydepd1"] & settings[p,"ydepd2"]){
          AllData <- cbind(AllData, FromCtoMatrix_settings(values = cmcmc[[p]],
                                                           p = p,
                                                           settings = settings,
                                                           yspecd1 = yspecd1[[p]],
                                                           yspecd2 = yspecd2[[p]],
                                                           Nums = Nums, Pois = Pois, Bins = Bins, Ords = Ords, Cats = Cats))
        }
        if(settings[p,"ydepd1"] & !settings[p,"ydepd2"]){
          AllData <- cbind(AllData, FromCtoMatrix_settings(values = cmcmc[[p]],
                                                           p = p,
                                                           settings = settings,
                                                           yspecd1 = yspecd1[[p]],
                                                           Nums = Nums, Pois = Pois, Bins = Bins, Ords = Ords, Cats = Cats))
        }
        if(!settings[p,"ydepd1"] & settings[p,"ydepd2"]){
          AllData <- cbind(AllData, FromCtoMatrix_settings(values = cmcmc[[p]],
                                                           p = p,
                                                           settings = settings,
                                                           yspecd2 = yspecd2[[p]],
                                                           Nums = Nums, Pois = Pois, Bins = Bins, Ords = Ords, Cats = Cats))
        }
        if(!settings[p,"ydepd1"] & !settings[p,"ydepd2"]){
          AllData <- cbind(AllData, FromCtoMatrix_settings(values = cmcmc[[p]],
                                                           p = p,
                                                           settings = settings,
                                                           Nums = Nums, Pois = Pois, Bins = Bins, Ords = Ords, Cats = Cats))
        }
      }
    }
    
    # chain ch completed
    mcmc$all <- rbind(mcmc$all, AllData)
  }
  
  return(mcmc)
}