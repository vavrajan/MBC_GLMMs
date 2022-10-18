### Function that takes subjects, their outcomes and regressors
# also takes previously generated parameters from posterior distribution
# and calculates probabilities of classification
# and most importantly even deviance of the model (based on supplied data)
# can be even applied to data not used for constructing MCMC

pUig_dev <- function(mcmc,      # MCMC generated parameters
                     #   output of Metropolis_within_Gibbs_MBC_NumBinOrdCat
                     Y,      # outcomes of subjects
                     X,      # regressors of subjects
                     Id,     # Identification numbers
                     # common to rows of Y and X
                     # when they should be taken together
                     dynamic_prob = F, # Should probabilities be calculated dynamically?
                     #dev = T,
                     start = mcmc$B+1,
                     end = mcmc$BM,
                     thin = 1,
                     ChainsToUse = 1:mcmc$Nchains,
                     tolerance = mcmc$tuning$double$tolerance, # tolerance when computing norm of the shift within Newton-Raphson
                     maxiter = mcmc$tuning$integer$maxiter, # maximum iterations allowed during Newton-Raphson update of proposal distribution
                     maxnrep = mcmc$tuning$integer$maxnrep,
                     NGQ = 1 # number of points for Adaptive Gaussian Quadrature
                     # default value 1 corresponds to Laplacian approximation
                     # do not use values higher than 7, since there are NGQ^totnran summands
){
  ### Every auxiliary values are stored within mcmc output
  # Such as whatsave, calc, spec, ...
  cat("\nR function pUig_dev succesfully called.\n")
  
  # Find roots and weights for Adaptive Gaussian Quadrature
  require("gaussquad")
  if(NGQ == 1){
    rules <- hermite.he.quadrature.rules(n = 2, normalized=FALSE)[[1]]
  }else{
    rules <- hermite.he.quadrature.rules(n = NGQ, normalized=FALSE)[[NGQ]]
  }
  
  # Iterations to be taken
  iterations <- seq(start, end, by = thin)
  niter <- length(iterations)
  
  # For probability calculation we need the following parameters:
  #   w, tau, beta, InvSigma
  # If those were not saved during the Metropolis_within_Gibbs_MBC_NumBinOrdCat function,
  # then we cannot use this function to calculate it.
  if(!(mcmc$whatsave["w"] & mcmc$whatsave["c_ord"] & mcmc$whatsave["tau_num"] 
       & mcmc$whatsave["beta_num_fix"] & mcmc$whatsave["beta_num"] 
       & mcmc$whatsave["beta_poi_fix"] & mcmc$whatsave["beta_poi"] 
       & mcmc$whatsave["beta_bin_fix"] & mcmc$whatsave["beta_bin"] 
       & mcmc$whatsave["beta_ord_fix"] & mcmc$whatsave["beta_ord"] 
       & mcmc$whatsave["beta_cat_fix"] & mcmc$whatsave["beta_cat"] 
       & mcmc$whatsave["InvSigma"])){
    stop("Some of the parameters w, c_ord, tau_num, betas, InvSigma was not saved.
Change the whatsave settings and use Metropolis_within_Gibbs_MBC_NumPoiBinOrdCat again.")
  }
  
  # Calculation will be done by C function, to speed up calculation
  # parameters need to be given to C function in a similar way to Metropolis_within_Gibbs_MBC_NumBinOrdCat
  
  n <- length(unique(Y[,Id]))
  N <- dim(Y)[1]
  UniqSubj <- unique(Y[,Id])
  n <- length(UniqSubj)
  NumUniqSubj <- c(1:n)
  names(NumUniqSubj) <- UniqSubj
  nsubj <- table(Y[,Id])
  max_n_i <- max(nsubj)
  
  Ys <- c(mcmc$Nums, mcmc$Pois, mcmc$Bins, mcmc$Ords, mcmc$Cats)
  nY <- sapply(list(mcmc$Nums, mcmc$Pois, mcmc$Bins, mcmc$Ords, mcmc$Cats), length)
  names(nY) <- c("Nums", "Pois", "Bins", "Ords", "Cats")
  
  
  if(mcmc$howsave == "matrix"){
    output <- data.frame()
  }
  
  if(mcmc$howsave == "list"){
    output <- list()
  }
  
  # Creating model matrix containing all needed columns for C
  finX <- data.frame(id = Y[,Id])
  Xcolnames <- numeric()
  fcols <- gcols <- rcols <- ocols <- list()
  for(y in Ys){
    fauxX <- model.matrix(mcmc$Formula[[y]]$fixed, X)
    gauxX <- model.matrix(mcmc$Formula[[y]]$group, X)
    rauxX <- model.matrix(mcmc$Formula[[y]]$random, X)
    if(is.element(y, mcmc$Ords)){
      fcols[[y]] <- colnames(fauxX)[-1]
      gcols[[y]] <- colnames(gauxX)[-1]
    }else{
      fcols[[y]] <- colnames(fauxX)
      gcols[[y]] <- colnames(gauxX)
    }
    fcols[[y]] <- setdiff(fcols[[y]], gcols[[y]])
    rcols[[y]] <- colnames(rauxX)
    
    notf <- setdiff(fcols[[y]], Xcolnames)
    Xcolnames <- c(Xcolnames, notf)
    addX <- data.frame(fauxX[,notf])
    colnames(addX) <- notf
    finX <- cbind(finX, addX)
    
    notg <- setdiff(gcols[[y]], Xcolnames)
    Xcolnames <- c(Xcolnames, notg)
    addX <- data.frame(gauxX[,notg])
    colnames(addX) <- notg
    finX <- cbind(finX, addX)
    
    notr <- setdiff(rcols[[y]], Xcolnames)
    Xcolnames <- c(Xcolnames, notr)
    addX <- data.frame(rauxX[,notr])
    colnames(addX) <- notr
    finX <- cbind(finX, addX)
    
    ocols[[y]] <- mcmc$Formula[[y]]$offset
    noff <- ifelse(mcmc$Formula[[y]]$offset == "", 0, 1)
    if(noff > 0){
      noto <- setdiff(ocols[[y]], Xcolnames)
      Xcolnames <- c(Xcolnames, noto)
      addX <- data.frame(X[,noto])
      colnames(addX) <- noto
      finX <- cbind(finX, addX)
    }
    
  }
  
  params <- c("beta_num_fix", "beta_num", "tau_num", 
              "beta_poi_fix", "beta_poi", 
              "beta_bin_fix", "beta_bin", 
              "beta_ord_fix", "beta_ord", "c_ord", 
              "beta_cat_fix", "beta_cat", 
              "InvSigma", 
              "w", "pUig_int", "deviance", "dev_i")
  ydepparams <- c(paste0("beta_", c("num", "poi", "bin", "ord", "cat"), "_fix"),
                  paste0("beta_", c("num", "poi", "bin", "ord", "cat")),
                  "c_ord")
  
  
  settings <- mcmc$settings
  rownames(settings)[rownames(settings)=="loglik"] <- "deviance"
  rownames(settings)[rownames(settings)=="pUig"] <- "pUig_int"
  settings["deviance", "save"] <- T
  settings["pUig_int", "save"] <- T
  settings["pUig_int",  c("d1", "dims", "dimswithG")] = c(ifelse(dynamic_prob,N,n), 
                                                          ifelse(dynamic_prob,N,n)*mcmc$G, 
                                                          ifelse(dynamic_prob,N,n)*mcmc$G)
  addset <- matrix(c(1,0,mcmc$G,0,0,0,0,0,0,0,0,niter,
              n, 0, 1, 0, 0, 1, 1, 
              n, n),
              nrow = 1)
  rownames(addset) <- "dev_i"
  colnames(addset) <- colnames(settings)
  settings <- rbind(settings, addset)
  
  settings <- settings[params,]
  settings$BM <- niter
  
  
  ### Preparation of parameters for Gibbs sampler in C
  ## Data
  # First column is going to be the id variable (0-th column)
  cId <- NumUniqSubj[as.character(Y[,Id])] - 1
  # -1 is there for C which works better with number beggining with 0
  # other columns  (beginning with 1st column)
  cY <- numeric()
  for(y in Ys){
    cY <- c(cY, as.numeric(as.character(Y[,y])))
    #Xcolnames <- unique(c(Xcolnames, Formula[[y]]$fr, Formula[[y]]$random))
  }
  cisYna <- is.na(cY)
  cY[cisYna] <- 0 # NA is not sent to C function
  
  cX <- numeric()
  # take only needed columns
  Xcolnums <- 1:length(Xcolnames)
  names(Xcolnums) <- Xcolnames
  for(x in Xcolnames){
    cX <- c(cX, finX[,x])
  } # no need for Id
  
  # Formula
  cFormulaF <- cFormulaG <- cFormulaR <- cFormulaO <- numeric()
  for(y in Ys){
    cFormulaF <- c(cFormulaF, Xcolnums[mcmc$lfixnames[[y]]])
    cFormulaG <- c(cFormulaG, Xcolnums[mcmc$lgrpnames[[y]]])
    cFormulaR <- c(cFormulaR, Xcolnums[mcmc$lrannames[[y]]])
    cFormulaO <- c(cFormulaO, Xcolnums[mcmc$loffnames[[y]]])
  }
  cFormulaF <- cFormulaF - 1 # index in C
  cFormulaG <- cFormulaG - 1 # index in C
  cFormulaR <- cFormulaR - 1 # index in C
  cFormulaO <- cFormulaO - 1 # index in C
  cnfix <- as.numeric(mcmc$nfix) # number of FIXED  regressors with variables y
  cngrp <- as.numeric(mcmc$ngrp) # number of GROUP-SPECIFIC  regressors with variables y
  cnran <- as.numeric(mcmc$nran) # number of RANDOM regressors with variables y
  cnoff <- as.numeric(mcmc$noff) # number of OFFSET regressors with variables y
  
  head(finX)
  
  for(ch in ChainsToUse){
    # arrays/fields where those variables are and will be stored
    cpUig_int <- double(niter*settings["pUig_int", "dimswithG"])
    cdeviance <- double(niter*settings["deviance", "dimswithG"])
    cdev_i <- double(niter*settings["dev_i", "dimswithG"])
    #if(dev){cdeviance <- double(niter*cdimswithK["deviance"])}else{cdeviance <- as.double(0)}
    
    if(mcmc$howsave == "list"){
      cbeta_num_fix <- FromListtoC_settings(mcmc[[ch]], iterations, "beta_num_fix", settings, yspecd1 = mcmc$yspecd1$beta_num_fix,
                                        Nums=mcmc$Nums, Pois=mcmc$Pois, Bins=mcmc$Bins, Ords=mcmc$Ords, Cats=mcmc$Cats)
      cbeta_num <- FromListtoC_settings(mcmc[[ch]], iterations, "beta_num", settings, yspecd1 = mcmc$yspecd1$beta_num,
                                        Nums=mcmc$Nums, Pois=mcmc$Pois, Bins=mcmc$Bins, Ords=mcmc$Ords, Cats=mcmc$Cats)
      ctau_num  <- FromListtoC_settings(mcmc[[ch]], iterations, "tau_num", settings,
                                        Nums=mcmc$Nums, Pois=mcmc$Pois, Bins=mcmc$Bins, Ords=mcmc$Ords, Cats=mcmc$Cats)
      cbeta_poi_fix <- FromListtoC_settings(mcmc[[ch]], iterations, "beta_poi_fix", settings, yspecd1 = mcmc$yspecd1$beta_poi_fix,
                                            Nums=mcmc$Nums, Pois=mcmc$Pois, Bins=mcmc$Bins, Ords=mcmc$Ords, Cats=mcmc$Cats)
      cbeta_poi <- FromListtoC_settings(mcmc[[ch]], iterations, "beta_poi", settings, yspecd1 = mcmc$yspecd1$beta_poi,
                                        Nums=mcmc$Nums, Pois=mcmc$Pois, Bins=mcmc$Bins, Ords=mcmc$Ords, Cats=mcmc$Cats)
      cbeta_bin_fix <- FromListtoC_settings(mcmc[[ch]], iterations, "beta_bin_fix", settings, yspecd1 = mcmc$yspecd1$beta_bin_fix,
                                        Nums=mcmc$Nums, Pois=mcmc$Pois, Bins=mcmc$Bins, Ords=mcmc$Ords, Cats=mcmc$Cats)
      cbeta_bin <- FromListtoC_settings(mcmc[[ch]], iterations, "beta_bin", settings, yspecd1 = mcmc$yspecd1$beta_bin,
                                        Nums=mcmc$Nums, Pois=mcmc$Pois, Bins=mcmc$Bins, Ords=mcmc$Ords, Cats=mcmc$Cats)
      cbeta_ord_fix <- FromListtoC_settings(mcmc[[ch]], iterations, "beta_ord_fix", settings, yspecd1 = mcmc$yspecd1$beta_ord_fix,
                                        Nums=mcmc$Nums, Pois=mcmc$Pois, Bins=mcmc$Bins, Ords=mcmc$Ords, Cats=mcmc$Cats)
      cbeta_ord <- FromListtoC_settings(mcmc[[ch]], iterations, "beta_ord", settings, yspecd1 = mcmc$yspecd1$beta_ord,
                                        Nums=mcmc$Nums, Pois=mcmc$Pois, Bins=mcmc$Bins, Ords=mcmc$Ords, Cats=mcmc$Cats)
      cc_ord    <- FromListtoC_settings(mcmc[[ch]], iterations, "c_ord", settings, yspecd1 = mcmc$yspecd1$c_ord,
                                        Nums=mcmc$Nums, Pois=mcmc$Pois, Bins=mcmc$Bins, Ords=mcmc$Ords, Cats=mcmc$Cats)
      cbeta_cat_fix <- FromListtoC_settings(mcmc[[ch]], iterations, "beta_cat_fix", settings, 
                                        yspecd1 = mcmc$yspecd1$beta_cat_fix, yspecd2 = mcmc$yspecd2$beta_cat_fix,
                                        Nums=mcmc$Nums, Pois=mcmc$Pois, Bins=mcmc$Bins, Ords=mcmc$Ords, Cats=mcmc$Cats)
      cbeta_cat <- FromListtoC_settings(mcmc[[ch]], iterations, "beta_cat", settings, 
                                        yspecd1 = mcmc$yspecd1$beta_cat, yspecd2 = mcmc$yspecd2$beta_cat,
                                        Nums=mcmc$Nums, Pois=mcmc$Pois, Bins=mcmc$Bins, Ords=mcmc$Ords, Cats=mcmc$Cats)
      cInvSigma <- FromListtoC_settings(mcmc[[ch]], iterations, "InvSigma", settings,
                                        Nums=mcmc$Nums, Pois=mcmc$Pois, Bins=mcmc$Bins, Ords=mcmc$Ords, Cats=mcmc$Cats)
      cw        <- FromListtoC_settings(mcmc[[ch]], iterations, "w", settings,
                                        Nums=mcmc$Nums, Pois=mcmc$Pois, Bins=mcmc$Bins, Ords=mcmc$Ords, Cats=mcmc$Cats)
    }
    if(mcmc$howsave == "matrix"){
      cbeta_num_fix <- FromMatrixtoC_settings(mcmc$all, chain = ch, iterations, "beta_num_fix", settings)
      cbeta_num <- FromMatrixtoC_settings(mcmc$all, chain = ch, iterations, "beta_num", settings)
      ctau_num  <- FromMatrixtoC_settings(mcmc$all, chain = ch, iterations, "tau_num", settings)
      cbeta_poi_fix <- FromMatrixtoC_settings(mcmc$all, chain = ch, iterations, "beta_poi_fix", settings)
      cbeta_poi <- FromMatrixtoC_settings(mcmc$all, chain = ch, iterations, "beta_poi", settings)
      cbeta_bin_fix <- FromMatrixtoC_settings(mcmc$all, chain = ch, iterations, "beta_bin_fix", settings)
      cbeta_bin <- FromMatrixtoC_settings(mcmc$all, chain = ch, iterations, "beta_bin", settings)
      cbeta_ord_fix <- FromMatrixtoC_settings(mcmc$all, chain = ch, iterations, "beta_ord_fix", settings)
      cbeta_ord <- FromMatrixtoC_settings(mcmc$all, chain = ch, iterations, "beta_ord", settings)
      cc_ord    <- FromMatrixtoC_settings(mcmc$all, chain = ch, iterations, "c_ord", settings)
      cbeta_cat_fix <- FromMatrixtoC_settings(mcmc$all, chain = ch, iterations, "beta_cat_fix", settings)
      cbeta_cat <- FromMatrixtoC_settings(mcmc$all, chain = ch, iterations, "beta_cat", settings)
      cInvSigma <- FromMatrixtoC_settings(mcmc$all, chain = ch, iterations, "InvSigma", settings)
      cw        <- FromMatrixtoC_settings(mcmc$all, chain = ch, iterations, "w", settings)
    }
    
    #dyn.load("Cfun/pUig_dev.dll")
    #dyn.unload("Cfun/pUig_dev.dll")
    
    print("About to start C computations.")
    # print(ch)
    # print(mcmc$G)
    # print(niter)
    # print(N)
    # print(n)
    # print(mcmc$nY)
    # print(cFormulaF)
    # print(cFormulaG)
    # print(cFormulaR)
    # print(cFormulaO)
    # print(cnfix)
    # print(cngrp)
    # print(cnran)
    # print(cnoff)
    # print(mcmc$totnran)
    # print(mcmc$Kord)
    # print(mcmc$Kcat)
    # print(settings$dims)
    # print(settings$dimswithG)
    # print(mcmc$tuning$integer$kspec_bi_cat)
    # print(cbeta_num_fix)
    # print(cbeta_num)
    # print(ctau_num)
    # print(cbeta_poi_fix)
    # print(cbeta_poi)
    # print(cbeta_bin_fix)
    # print(cbeta_bin)
    # print(cbeta_ord_fix)
    # print(cbeta_ord)
    # print(cc_ord)
    # print(cbeta_cat_fix)
    # print(cbeta_cat)
    # print(cInvSigma)
    # print(cw)
    # print(cpUig_int)
    # print(cdeviance)
    # print(cdev_i)
    # print(dynamic_prob)
    # print(NGQ)
    # print(rules$x)
    # print(rules$w)
    # print(tolerance)
    # print(maxiter)
    # print(maxnrep)
    
    #system.time(
    cmcmc <-
      .C("pUig_dev",
         Id        = as.integer(cId),
         Y         = as.double(cY),
         isYna     = as.integer(cisYna),
         X         = as.double(cX),
         spec      = as.integer(mcmc$spec),
         # parameters describing dimensions
         chain     = as.integer(ch), # number of the chain
         G         = as.integer(mcmc$G), # number of classes
         BM        = as.integer(niter), # total number of generated states
         N         = as.integer(N), # total number of observations
         n         = as.integer(n), # total number of subjects (different ids in the dataset)
         nY        = as.integer(mcmc$nY), # 3 numbers: counts of Nums, Ords and Bins variables
         FormulaF  = as.integer(cFormulaF), # numbers of columns of X that should be used for FIXED  effects of modelled responses
         FormulaG  = as.integer(cFormulaG), # numbers of columns of X that should be used for GROUP-SPECIFIC  effects of modelled responses
         FormulaR  = as.integer(cFormulaR), # numbers of columns of X that should be used for RANDOM effects of modelled responses
         FormulaO  = as.integer(cFormulaO), # numbers of columns of X that should be used for OFFSET effects of modelled responses
         nfix      = as.integer(cnfix), 
         ngrp      = as.integer(cngrp), 
         nran      = as.integer(cnran),
         noff      = as.integer(cnoff),
         totnran   = as.integer(mcmc$totnran),
         Kord      = as.integer(mcmc$Kord), # the counts of categories of ordinal variables (-1)
         Kcat      = as.integer(mcmc$Kcat), # the counts of categories of categorical variables (-1)
         dims      = as.integer(settings$dims), # the length of subarray that corresponds to one state (disected by various parameters)
         dimswithG = as.integer(settings$dimswithG), # the length of subarray that corresponds to one state
         kspec_bi_cat = as.integer(mcmc$tuning$integer$kspec_bi_cat),
         #predictor = double(N*sum(nY)),
         # the function should count totnran, totnfix and cumsum versions of nfix and nran
         # arrays to store generated states
         beta_num_fix= as.double(cbeta_num_fix),
         beta_num    = as.double(cbeta_num),
         tau_num     = as.double(ctau_num),
         beta_poi_fix= as.double(cbeta_poi_fix),
         beta_poi    = as.double(cbeta_poi),
         beta_bin_fix= as.double(cbeta_bin_fix),
         beta_bin    = as.double(cbeta_bin),
         beta_ord_fix= as.double(cbeta_ord_fix),
         beta_ord    = as.double(cbeta_ord),
         c_ord       = as.double(cc_ord),
         beta_cat_fix= as.double(cbeta_cat_fix),
         beta_cat    = as.double(cbeta_cat),
         InvSigma    = as.double(cInvSigma),
         w           = as.double(cw),
         pUig_int    = as.double(cpUig_int),
         deviance    = as.double(cdeviance),
         dev_i       = as.double(cdev_i),
         dynamic_prob = as.integer(dynamic_prob),
         NGQ = as.integer(NGQ),
         roots = as.double(rules$x),
         weights = as.double(rules$w),
         tolerance = as.double(tolerance),
         maxiter = as.integer(maxiter),
         maxnrep = as.integer(maxnrep)
      )

    if(mcmc$howsave == "list"){
      # results will be returned in structured list - the same way as original R function
      chain <- list()
      for(p in c("deviance", "dev_i", "pUig_int")){
        if(settings[p, "save"]){
          chain[[p]] <- FromCtoList_settings(cmcmc[[p]], p, settings,
                                             Nums=mcmc$Nums, Pois=mcmc$Pois, Bins=mcmc$Bins, Ords=mcmc$Ords, Cats=mcmc$Cats)
        }
      }
      chain[["m"]] <- iterations

      output[[ch]] <- chain
    }

    if(mcmc$howsave == "matrix"){
      # results will be returned in structured list - the same way as original R function
      AllData <- matrix(ch, nrow = niter, ncol = 1)
      AllData <- cbind(AllData, matrix(iterations, nrow = niter, ncol = 1))
      colnames(AllData) <- c("chain", "m")

      for(p in c("deviance", "dev_i", "pUig_int")){
        if(settings[p, "save"]){
          AllData <- cbind(AllData, FromCtoMatrix_settings(cmcmc[[p]], p, settings,
                                                           Nums=mcmc$Nums, Pois=mcmc$Pois, Bins=mcmc$Bins, Ords=mcmc$Ords, Cats=mcmc$Cats))
        }
      }
      

      output <- rbind(output, AllData)
    }
    
  }  # end of for ch in 1:Nchains
  
  return(output)
  
}
