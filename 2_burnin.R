# mtDNA Project: Should mitochondrial variation be accounted for on BVE
# Dairy Cattle Breeding Simulations
# Roslin Institute
# from July 2020
# Burn-in: Fill-in breeding population with its categories and sub-categories 
#            and run progeny test schemefor 20 generations


# save founder population records to be able to trace mitogenome
generation = 0
Records <- NULL
Covars <- NULL
program = "burnin"
model = "std"

# start records and save results for (co)variance
Records <- recording(Records, founders, mtfile)
Covars <- runCovars(Covars, Records)

# Estabilish breeding population with categories and sub-categories
breeding <- structure(list("youngBulls" = vector("list", 2), "progenyTest" = vector("list", 4), "eliteSire" = vector("list", 4),
                           "youngFemales" = vector("list", 3), "Cows" = vector("list", 4), "eliteDam" = vector("list", 4)),
                      class = 'Population')

# 10 generations required to fill-in all categories
# Because new levels are added to the end of the list, smaller numbers represent older groups (1-oldest , 4-newest)
generation = 1
breeding$eliteSire[[1]] <- randCross(founders, nES); breeding$eliteSire[[1]]@sex[] <- "M"
Records <- recording(Records, breeding$eliteSire[[1]], mtfile)
Covars <- runCovars(Covars, Records)

generation = 2
breeding$eliteSire[[2]] <- randCross(founders, nES); breeding$eliteSire[[2]]@sex[] <- "M"
Records <- recording(Records, breeding$eliteSire[[2]], mtfile)
Covars <- runCovars(Covars, Records)

generation = 3
breeding$eliteSire[[3]] <- randCross(founders, nES); breeding$eliteSire[[3]]@sex[] <- "M"
Records <- recording(Records, breeding$eliteSire[[3]], mtfile)
Covars <- runCovars(Covars, Records)

generation = 4
breeding$eliteSire[[4]] <- randCross(founders, nES); breeding$eliteSire[[4]]@sex[] <- "M"
breeding$eliteDam[[1]] <- randCross(founders, nED*0.34); breeding$eliteDam[[1]]@sex[] <- "F" #nED*0.34
breeding$Cows[[1]] <- randCross(founders, nCW*0.34); breeding$Cows[[1]]@sex[] <- "F" #nCW*0.34
Records <- recording(Records, c(breeding$eliteSire[[4]],
                                breeding$eliteDam[[1]],
                                breeding$Cows[[1]]), mtfile)
Covars <- runCovars(Covars, Records)
Records$Lactation[Records$Generation == generation & Records$Gender == "F"] <- 5

generation = 5
breeding$progenyTest[[1]] <- randCross(founders, nPTB); breeding$progenyTest[[1]]@sex[] <- "M"
breeding$eliteDam[[2]] <- randCross(founders, nED*0.49); breeding$eliteDam[[2]]@sex[] <- "F" #nED*0.49
breeding$Cows[[2]] <- randCross(founders, nCW*0.49); breeding$Cows[[2]]@sex[] <- "F" #nCW*0.49
Records <- recording(Records, c(breeding$progenyTest[[1]],
                                breeding$eliteDam[[2]],
                                breeding$Cows[[2]]), mtfile)
Covars <- runCovars(Covars, Records)
Records$Lactation[Records$Generation == generation & Records$Gender == "F"] <- 4

generation = 6
breeding$progenyTest[[2]] <- randCross(founders, nPTB); breeding$progenyTest[[2]]@sex[] <- "M"
breeding$eliteDam[[3]] <- randCross(founders, nED*0.70); breeding$eliteDam[[3]]@sex[] <- "F" #nED*0.70
breeding$Cows[[3]] <- randCross(founders, nCW*0.70); breeding$Cows[[3]]@sex[] <- "F" #nCW*0.70
Records <- recording(Records, c(breeding$progenyTest[[2]],
                                breeding$eliteDam[[3]],
                                breeding$Cows[[3]]), mtfile)
Covars <- runCovars(Covars, Records)
Records$Lactation[Records$Generation == generation & Records$Gender == "F"] <- 3

generation = 7
breeding$progenyTest[[3]] <- randCross(founders, nPTB); breeding$progenyTest[[3]]@sex[] <- "M"
breeding$eliteDam[[4]] <- randCross(founders, nED); breeding$eliteDam[[4]]@sex[] <- "F" #nED
breeding$Cows[[4]] <- randCross(founders, nCW); breeding$Cows[[4]]@sex[] <- "F" #nCW
Records <- recording(Records, c(breeding$progenyTest[[3]],
                                breeding$eliteDam[[4]],
                                breeding$Cows[[4]]), mtfile)
Covars <- runCovars(Covars, Records)
Records$Lactation[Records$Generation == generation & Records$Gender == "F"] <- 2

generation = 8
breeding$progenyTest[[4]] <- randCross(founders, nPTB); breeding$progenyTest[[4]]@sex[] <- "M"
breeding$youngFemales[[1]] <- randCross(founders, nFM); breeding$youngFemales[[1]]@sex[] <- "F"
Records <- recording(Records, c(breeding$progenyTest[[4]],
                                breeding$youngFemales[[1]]), mtfile)
Covars <- runCovars(Covars, Records)
Records$Lactation[Records$Generation == generation & Records$Gender == "F"] <- 1

generation = 9
breeding$youngBulls[[1]] <- randCross(founders, nYB); breeding$youngBulls[[1]]@sex[] <- "M"
breeding$youngFemales[[2]] <- randCross(founders, nFM); breeding$youngFemales[[2]]@sex[] <- "F"
Records <- recording(Records, c(breeding$youngBulls[[1]],
                                breeding$youngFemales[[2]]), mtfile)
Covars <- runCovars(Covars, Records)

generation = 10
breeding$youngBulls[[2]] <- randCross(founders, nYB); breeding$youngBulls[[2]]@sex[] <- "M"
breeding$youngFemales[[3]] <- randCross(founders, nFM); breeding$youngFemales[[3]]@sex[] <- "F"
Records <- recording(Records, c(breeding$youngBulls[[2]],
                                breeding$youngFemales[[3]]), mtfile)
Covars <- runCovars(Covars, Records)

# Estimate breeding values to start evaluation scenario
preparePAR("PED")
Records = runBLUP(Records)
file.remove("Blupf901.dat")

nBurninGen = 20
program = "burnin"

#------------ CHANGE TO CALL BREEDING SCHEME SCRIPT ------------#

Accuracy = NULL
for(i in (generation+1):((generation+1) + nBurninGen)){
  cat("Working on Generation:", i, "...\n")
  
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
  
  breeding$maleCalfs <- selectInd(offs0, sex="M", (offs0@nInd*0.97)/2 , use="rand")
  breeding$femCalfs <- selectInd(c(offs0, offs1, offs2), sex="F", (c(offs0, offs1, offs2)@nInd - 2500)*0.97, use="rand")
  # update to current year newborns (3% death rate of calfs)
  
  # every year: update category
  breeding$youngBulls <- c(breeding$youngBulls, selectInd(breeding$maleCalfs, nYB, use="rand"))
  breeding$youngFemales <- c(breeding$youngFemales, selectInd(breeding$femCalfs, nFM, use="rand"))
  # move calfs to young bulls and young females category
  
  generation = generation + 1
  Records <- recording(Records, c(breeding$youngBulls[[3]], breeding$youngFemales[[4]]), mtfile)
  # record selected offspring (new generation)
  
  # SELECTION FOR PROGENY TEST: every year, candidates from 2 GA
  # young bulls reach puberty and enter progeny test
  breeding$youngBulls[[1]]@ebv <- as.matrix(with(Records, tEbv[match(breeding$youngBulls[[1]]@id, Records$IId)]))
  # estimate EBV for young bulls (without offspring - low accuracy)
  
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
rm(dams, sires, multiplier, testing, offs0, offs1, offs2, FOUNDERPOP, founders)

# correcting scales for bias estimation
s <- c(sd(Records$nTbv[Records$Generation == (generation - nBurninGen)]),
       sd(Records$mTbv[Records$Generation == (generation - nBurninGen)]),
       sd(Records$tTbv[Records$Generation == (generation - nBurninGen)]))
# estimate sd for first observation for each genome
Accuracy <- Accuracy %>% mutate(b0_n = b0_n/s[1],
                                b0_m = b0_m/s[2],
                                b0_t = b0_t/s[3])
# divide b0 results by sd


# save data
Records <- Records %>% filter(Generation >= generation-10)
Records <- Records %>% mutate(gv_corr = nTbv - mean(nTbv))
gv <- mean(Records$nTbv)

file.remove("Blupf901.dat")
save.image("burnin.RData")
