# mtDNA Project: Should mitochondrial variation be accounted for on BVE
# Dairy Cattle Breeding Simulations
# Roslin Institute
# from July 2020
# Dairy Cattle Single Poligenic Trait with multiple observations

for(i in (generation+1):((generation+1) + nBreedingGen)){
  cat("Working on Generation:", i, "Program:", program, "Model:", model, "...\n")
    
  # define mating groups for each category:
  dams <- do.call(c, unlist(breeding$eliteDam, recursive=FALSE))
  sires <- do.call(c, unlist(breeding$eliteSire, recursive=FALSE))
  testing <- do.call(c, unlist(breeding$progenyTest, recursive=FALSE))
  multiplier <- do.call(c, unlist(breeding$Cows, recursive=FALSE))
    
  # Reproduction happens every year
  offs0 <- randCross2(dams, sires, nCrosses=dams@nInd, nProgeny=2)
  offs1 <- randCross2(multiplier, testing, nCrosses=multiplier@nInd, nProgeny=1)
  offs2 <- randCross2(breeding$youngFemales[[1]], testing, nCrosses=nFM, nProgeny=1)
    
  offs1@sex[] <- rep("F", offs1@nInd); offs2@sex[] <- rep("F", offs2@nInd)
  # Setting multiplier offspring to only females (males won't enter the system)
  
  # update to current year newborns
  breeding$maleCalfs <- selectInd(offs0, sex="M", (offs0@nInd*0.97)/2 , use="rand")
  # male calves origin only from damxsire mating --- considering 3% will not make to this stage
  breeding$femCalfs <- selectInd(c(offs0, offs1, offs2), sex="F", (c(offs0, offs1, offs2)@nInd - 2500)*0.97, use="rand")
  # female calves origin from all matings --- considering 3% will not make to this stage
    
  # every year: update category
  breeding$youngBulls <- c(breeding$youngBulls, selectInd(breeding$maleCalfs, nYB, use="rand"))
  breeding$youngFemales <- c(breeding$youngFemales, selectInd(breeding$femCalfs, nFM, use="rand"))
  # move calves to young bulls and young females category
  
  generation = generation + 1
  Records <- recording(Records, c(breeding$youngBulls[[3]], breeding$youngFemales[[4]]), mtfile)
  Records$gv_corr[Records$Generation==generation] <- Records$nTbv[Records$Generation==generation] - gv 
  # record selected offspring (new generation)
  if(program=="GEN"){
      snpData(c(breeding$youngBulls[[3]], breeding$youngFemales[[4]][dams@id %in% breeding$youngFemales[[4]]@mother]))
  }
  # if genomic scheme uppdate snp file
  
  # SELECTION FOR PROGENY TEST: every year, candidates from 2 GA
  # young bulls reach puberty and enter progeny test
  breeding$youngBulls[[1]]@ebv <- as.matrix(with(Records, tEbv[match(breeding$youngBulls[[1]]@id, Records$IId)]))
  # estimate EBV for young bulls
    
  category <- "young_bulls"
  Accuracy <- record.accuracy(Accuracy, breeding$youngBulls[[1]], Records, model)
    
  breeding$progenyTest <- c(breeding$progenyTest, selectInd(breeding$youngBulls[[1]], nPTB, use="ebv"))
  # select young bulls to enter progeny test based on EBV
    
  breeding$youngBulls[[1]] <- NULL
  # not selected are culled
    
  # LACTATION END (2GA)
  Records$Lactation[Records$Lactation == 5] = NA
  Records$Lactation[Records$Lactation == 4] = 5
  Records$Lactation[Records$Lactation == 3] = 4
  Records$Lactation[Records$Lactation == 2] = 3
  Records$Lactation[Records$Lactation == 1] = 2
  Records$Lactation[Records$IId %in% breeding$youngFemales[[1]]@id] = 1
  
  Records = runBLUP(Records)
  
  breeding$youngFemales[[1]]@ebv <- as.matrix(with(Records, tEbv[match(breeding$youngFemales[[1]]@id, Records$IId)]))
  breeding$youngFemales[[2]]@ebv <- as.matrix(with(Records, tEbv[match(breeding$youngFemales[[2]]@id, Records$IId)]))
  dams@ebv <- as.matrix(with(Records, tEbv[match(dams@id, Records$IId)]))
  multiplier@ebv <- as.matrix(with(Records, tEbv[match(multiplier@id, Records$IId)]))
  # set EBV for accuracy check
  
  #-----Accuracy of three female satges: compare efficiency of selecting cows----------------------------------------#
  category <- "heifers"
  Accuracy <- record.accuracy(Accuracy, breeding$youngFemales[[2]], Records, model)
  category <- "1st_lact"
  Accuracy <- record.accuracy(Accuracy, breeding$youngFemales[[1]], Records, model)
  category <- "cows"
  Accuracy <- record.accuracy(Accuracy, c(dams, multiplier), Records, model)
  #-----------------------------------------------------------------------------------------------------------------#

  # moving categories
  for(j in 1:4){
    breeding$eliteDam[[j]] <- selectInd(breeding$eliteDam[[j]], nInd = 0.7*breeding$eliteDam[[j]]@nInd, use="rand")
    breeding$Cows[[j]] <- selectInd(breeding$Cows[[j]], nInd = 0.7*breeding$Cows[[j]]@nInd, use="rand")
  }
  # 30% of females from all categories will be culled
  
  breeding$youngFemales[[1]]@pheno <- as.matrix(with(Records, PhenoL1[match(breeding$youngFemales[[1]]@id, Records$IId)]))
  # set heifers phenotypes to be able to select
  
  breeding$eliteDam <- c(breeding$eliteDam, selectInd(breeding$youngFemales[[1]], nED, use="pheno"))
  # best young females move to elite Dams category
  
  tmp <- breeding$youngFemales[[1]][setdiff(breeding$youngFemales[[1]]@id, breeding$eliteDam[[5]]@id)]
  breeding$Cows <- c(breeding$Cows, selectInd(tmp, nCW, use="pheno"))
  rm(tmp)
  # females that don't make to elite go to multiplier
    
  breeding$eliteDam[[1]] <- NULL; breeding$Cows[[1]] <- NULL; breeding$youngFemales[[1]] <- NULL
  # Elite Dams and Multipliers that end 5th lac are culled; FM[[1]] category is excluded
  
  active <- c(do.call(c, unlist(breeding$eliteDam, recursive=FALSE))@id, 
              do.call(c, unlist(breeding$Cows, recursive=FALSE))@id,
              do.call(c, unlist(breeding$youngFemales, recursive=FALSE))@id)
  Records$Lactation[!(Records$IId %in% active)] = NA
  rm(active)
  
  # every 4 yrs: END OF PROGENY TEST
  # Generating EBV for progeny test bulls
  # Selecting best progeny tested bulls based on EBV
  breeding$progenyTest[[1]]@ebv <- as.matrix(with(Records, tEbv[match(breeding$progenyTest[[1]]@id, Records$IId)]))
  # set EBV for tested bulls (higher accuracies)
  
  category <- "proven_bulls"
  Accuracy <- record.accuracy(Accuracy, breeding$progenyTest[[1]], Records, model)
    
  breeding$eliteSire <- c(breeding$eliteSire, selectInd(breeding$progenyTest[[1]], nES, use="ebv"))
  # Elite sire selected from tested bulls (4GA)
  breeding$progenyTest[[1]] <- NULL; breeding$eliteSire[[1]] <- NULL
  # Non-selected PTB and Elite sire @ age 10: Involuntary culling
  
  Covars <- runCovars(Covars, Records)
}


# MOVE TO RUNNING SCRIP
rm(dams, sires, multiplier, testing, breeding, offs0, offs1, offs2)
s <- c(sd(Records$nTbv[Records$Generation == (generation - nBreedingGen)]),
        sd(Records$mTbv[Records$Generation == (generation - nBreedingGen)]),
        sd(Records$tTbv[Records$Generation == (generation - nBreedingGen)]))

Accuracy <- Accuracy %>% mutate(b0_n = ifelse(Generation >= 32, b0_n/s[1], b0_n),
                                b0_m = ifelse(Generation >= 32, b0_m/s[2], b0_m),
                                b0_t = ifelse(Generation >= 32, b0_t/s[3], b0_t))


