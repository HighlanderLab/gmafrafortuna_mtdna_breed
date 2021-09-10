# Evaluation of the impact of accounting for mitochondrial variation on breeding values estimation for dairy cattle
# Gabriela Mafra Fortuna
# Highlander Lab
# The Roslin Institute 
# July 2020 - updated Aug 2021
# ------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------ Run Breeding Programme -----------------------------------------------
ids <- NULL

# Generate parameter file for BLUPF90
#preparePAR(paste0(program, model))

# If genomic model, generate initial snp file
if(program=="GEN"){
    file.remove("snp.dat")
    snpData(active)
}

for(i in (generation+1):((generation+1)+nBreeding)){
  cat("Breeding generation", i, "Replicate ", run, "...\n",
      "Program:", program, "Model:", model, "Selection:", selection, "Scenario:", traitScen, "\n")
  
  # define matings for each category:
  dams <- do.call(c, unlist(pop$eliteDams, recursive = FALSE))
  sires <- do.call(c, unlist(pop$eliteSires, recursive = FALSE))
  testing <- do.call(c, unlist(pop$waitingBulls, recursive = FALSE))
  multiplier <- do.call(c, unlist(pop$commercial, recursive = FALSE))
  
  # Reproduction happens every year
  offs0 <- randCross2(dams, sires, nCrosses = dams@nInd, nProgeny = concepRateC)
  offs1 <- randCross2(multiplier, testing, nCrosses = multiplier@nInd, nProgeny = concepRateC)
  offs2 <- randCross2(pop$heifers[[1]], testing, nCrosses = pop$heifers[[1]]@nInd, nProgeny = concepRateH)
  
  # replacement males came only from elite crossing. all other crosses generate only females
  offs1@sex[] <- "F" ; offs2@sex[] <- "F"
  
  # newborns are assign to first age category in the pop list
  pop$youngBulls <- c(pop$youngBulls, 
                      selectInd(offs0, sex="M", nYoungBulls, use="rand"))

  pop$heifers <- c(pop$heifers,
                   c(selectInd(c(offs1, offs2), sex = "F", nHeifers, use="rand"), # waiting bulls progeny to secure accuracy
                     selectInd(offs0, sex="F", nReplacement, use="rand"))) # elite crossing progeny
  
  cat("Waiting bull progeny =", mean(table(pop$heifers[[4]]@father)), "\n")
  
  # Genotype new pop
  if(program=="GEN"){ 
    snpData(c(pop$youngBulls[[3]], pop$heifers[[4]][pop$heifers[[4]]@mother %in% dams@id]))
    
    geno.id <- c(geno.id, 
                 pop$heifers[[4]][pop$heifers[[4]]@mother %in% dams@id]@id)
  }
  
  # Update year's records
  generation = generation + 1
  Records <- recording(Records, mtfile, 
                       c(pop$youngBulls[[3]], pop$heifers[[4]]))
  Covars <- runCovars(Covars, Records, 
                      c(pop$youngBulls[[3]], pop$heifers[[4]]))
  
  # Pre-test: select young males to enter progeny testing
  pop$youngBulls[[1]]@ebv <- as.matrix(with(Records,
                                            if(selection=="extreme"){tEbv[match(pop$youngBulls[[1]]@id, 
                                                                                Records$IId)]}else{nEbv[match(pop$youngBulls[[1]]@id, 
                                                                                                              Records$IId)]}))
  # assign breeding values (based on selection option: nEBV or tEBV) to be able to select
  
  category = "young_bulls"
  Accuracy <- record.accuracy(Accuracy, 
                              pop$youngBulls[[1]], Records, model)
  youngbulls <- pop$youngBulls[[1]]
  # estimate accuracy
  
  pop$waitingBulls <- c(pop$waitingBulls,
                        selectInd(pop$youngBulls[[1]], nWaitingBulls, use="ebv"))
  # select top candidates to enter progeny test
  
  pop$youngBulls[[1]] <- NULL
  # all non-selected candidates are culled
  
  # End of lactation
  Records$Lactation[Records$Lactation == 5] = NA
  Records$Lactation[Records$Lactation == 4] = 5
  Records$Lactation[Records$Lactation == 3] = 4
  Records$Lactation[Records$Lactation == 2] = 3
  Records$Lactation[Records$Lactation == 1] = 2
  Records$Lactation[Records$IId %in% pop$heifers[[1]]@id] = 1
  
  # run evaluation
  preparePAR(paste0(program, model))
  runRENUM(Records)
  varComp(model)
  Records = runBLUP(Records)
  
  pop$heifers[[1]]@ebv <- as.matrix(with(Records,
                                         if(selection=="baseline"){nEbv[match(pop$heifers[[1]]@id, 
                                                                              Records$IId)]}else{tEbv[match(pop$heifers[[1]]@id, 
                                                                                                            Records$IId)]}))
  
  
  # assign breeding values (based on selection option: nEBV or tEBV) to be able to select
  
  if(program=="PED"){
    category <- "heifers"
    Accuracy <- record.accuracy(Accuracy,
                                pop$heifers[[2]], Records, model)
    Bias <- record.bias(Bias,
                        pop$heifers[[2]], Records, model)
    category <- "1st_lact"
    Accuracy <- record.accuracy(Accuracy,
                                pop$heifers[[1]], Records, model)
    category <- "cows"
    Accuracy <- record.accuracy(Accuracy, 
                                c(dams, multiplier), Records, model)
    
  }else{
    if(sum(pop$heifers[[2]]@id %in% geno.id) != 0){
    category <- "heifers_geno"
    Accuracy <- record.accuracy(Accuracy,
                                pop$heifers[[2]][pop$heifers[[2]]@id %in% geno.id],
                                Records, model)
    }
      
    category <- "heifers_no_geno"
    Accuracy <- record.accuracy(Accuracy,
                                  pop$heifers[[2]][!(pop$heifers[[2]]@id %in% geno.id)],
                                  Records, model)
    
    if(sum(pop$heifers[[1]]@id %in% geno.id) != 0){
      category <- "1st_lact_geno"
      Accuracy <- record.accuracy(Accuracy,
                                    pop$heifers[[1]][pop$heifers[[1]]@id %in% geno.id],
                                    Records, model)
    }
      
    category <- "1st_lact_no_geno"
    Accuracy <- record.accuracy(Accuracy,
                                  pop$heifers[[1]][!(pop$heifers[[1]]@id %in% geno.id)],
                                  Records, model)
      
    cows <- c(dams, multiplier)
    if(sum(cows@id %in% geno.id) != 0){
      category <- "cows_geno"
      Accuracy <- record.accuracy(Accuracy,
                                cows[cows@id %in% geno.id],
                                Records, model)
  }
    
  category <- "cows_no_geno"
  Accuracy <- record.accuracy(Accuracy,
                              cows[!(cows@id %in% geno.id)],
                              Records, model)
  rm(cows)
  }  
  
  category <- "heifers"
  Bias <- record.bias(Bias, 
                      pop$heifers[[2]],
                      Records, model)
  category <- "cows"
  Bias <- record.bias(Bias, 
                      c(pop$heifers[[1]], dams, multiplier),
                      Records, model)
  # estimate accuracies
  
  # moving categories
  for(j in 2:4){
    pop$eliteDams[[j]]@ebv <- as.matrix(with(Records,
                                             if(selection=="baseline"){nEbv[match(pop$eliteDams[[j]]@id, 
                                                                                  Records$IId)]}else{tEbv[match(pop$eliteDams[[j]]@id, 
                                                                                                                Records$IId)]}))
    pop$commercial[[j]]@ebv <- as.matrix(with(Records,
                                              if(selection=="baseline"){nEbv[match(pop$commercial[[j]]@id, 
                                                                                   Records$IId)]}else{tEbv[match(pop$commercial[[j]]@id, 
                                                                                                                 Records$IId)]}))
    # assign ebvs (based on selection option: nEBV or tEBV)
    
    pop$eliteDams[[j]] <- selectInd(pop$eliteDams[[j]],
                                     nInd = (1-cullRate)*pop$eliteDams[[j]]@nInd,
                                     use="ebv")
    pop$commercial[[j]] <- selectInd(pop$commercial[[j]],
                                     nInd = (1-cullRate)*pop$commercial[[j]]@nInd, 
                                     use="ebv")
    # select best females
  }
  # females move up one age group and 30% are culled every year
  
  pop$eliteDams <- c(pop$eliteDams,
                      selectInd(pop$heifers[[1]],
                                nInd = nEliteDams,
                                use = "ebv"))
  # best 1st lactation cows are kept as elite dams
  
  tmp <- pop$heifers[[1]][setdiff(pop$heifers[[1]]@id,
                                  pop$eliteDams[[5]]@id)]
  # 1st lactation not selected as elite dams
  
  pop$commercial <- c(pop$commercial,
                      selectInd(tmp, nInd = nCommercial,
                                use = "ebv"))
  # 45% of non-selected 1st lactation cows are kept as commercial
  
  rm(tmp); pop$heifers[[1]] <- NULL; pop$eliteDamss[[1]] <- NULL; pop$commercial[[1]] <- NULL
  # 30% non-selected 1st lactation cows are culled. All cows that end 5th lactation are culled.
  
  active <- c(do.call(c, unlist(pop$eliteDams, recursive = FALSE))@id,
              do.call(c, unlist(pop$commercial, recursive = FALSE))@id,
              do.call(c, unlist(pop$heifers, recursive = FALSE))@id)
  # gather id of all active cows
  Records$Lactation[!(Records$IId %in% active)] = NA
  rm(active)
  # update data file: culled cows are "removed"
  
  # End of progeny-testing: best bulls are selected on EBV
  pop$waitingBulls[[1]]@ebv <- as.matrix(with(Records,
                                              if(selection=="extreme"){tEbv[match(pop$waitingBulls[[1]]@id, 
                                                                                  Records$IId)]}else{nEbv[match(pop$waitingBulls[[1]]@id, 
                                                                                                                Records$IId)]}))
  # assign ebv (based on selection option: nEBV or tEBV) to be able to select
  
  category <- "proven_bulls"
  Accuracy <- record.accuracy(Accuracy,
                              pop$waitingBulls[[1]], Records, model)
  category <- "males"
  Bias <- record.bias(Bias, 
                      c(youngbulls,
                        pop$waitingBulls[[1]]), Records, model)
  # estimate accuracy
  
  pop$eliteSires <- c(pop$eliteSires,
                      selectInd(pop$waitingBulls[[1]], nInd = nEliteSires, 
                                use="ebv"))
  # select best performing bulls to become elite sires
  
  pop$waitingBulls[[1]] <- NULL ; pop$eliteSires[[5]] <- NULL
  # non-selected candidates are culled. Sires used for 5 years are culled. 
  
}
# ------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------ Correct scales for bias estimation ------------------------------------------ 

s <- c(sd(Records$nTbv[Records$Generation == (generation - nBreeding)]),
       sd(Records$mTbv[Records$Generation == (generation - nBreeding)]),
       sd(Records$tTbv[Records$Generation == (generation - nBreeding)]))
# estimate sd for first observation for each genome
Bias <- Bias %>% mutate(b0_n = ifelse(Generation >= 32, b0_n/s[1], b0_n),
                        b0_m = ifelse(Generation >= 32, b0_m/s[2], b0_m),
                        b0_t = ifelse(Generation >= 32, b0_t/s[3], b0_t))
# divide b0 results by sd
