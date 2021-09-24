# Evaluation of the impact of accounting for mitochondrial variation on breeding values estimation for dairy cattle
# Gabriela Mafra Fortuna
# Highlander Lab
# The Roslin Institute 
# July 2020 - updated Aug 2021
# ------------------------------------------------------------------------------------------------------------------------
# ----------------------------------- create mt population depending on trait scenario -----------------------------------
cat("Generating mitochondrial pop ", traitScen, "\n")
mtDNA <- newPop(mtDNA, simParam = if(traitScen=="maxQTL"){SP2}else{SP3})

# Store haplotypes, define maternal lineages and haplogroups
mtfile <- as.matrix(pullSegSiteHaplo(mtDNA, simParam = if(traitScen=="maxQTL"){SP2}else{SP3}))
mtfile <- as_tibble(mtfile)
mtfile <- mtfile %>% unite(haplotype, 1:length(mtfile), sep="")
mtfile <- mtfile %>% add_column(ML = mtDNA@id,
                                mTbv = mtDNA@gv[,1]) %>% 
  filter(!duplicated(haplotype))

mtDNAx <- mtDNA
# save original pop
mtDNA <- selectInd(mtDNA, nInd = nrow(mtfile), simParam = if(traitScen=="maxQTL"){SP2}else{SP3},
                   candidates = mtfile$ML, use="rand")

# Generate mtDNA inverse matrix based on trait scenario
cat("Generating mtDNA inverse matrix...\n")
mt_ref = NULL
mtdnaGinv(mtDNA, if(traitScen=="maxQTL"){SP2}else{SP3})

# ------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------- Fill in animal categories ----------------------------------------------
generation = 0
program = "PED"
model = "std"
selection = "baseline"

# start recording
Records <- recording(Records, mtfile, founders)

# 10 generations required to fill-in all categories
# Because new levels are added to the end of the list, smaller numbers represent older groups (1-oldest , 4-newest)
generation = 1
pop$eliteSires[[1]] <- randCross(founders, nEliteSires); pop$eliteSires[[1]]@sex[] <- "M"
Records <- recording(Records, mtfile, pop$eliteSires[[1]])

generation = 2
pop$eliteSires[[2]] <- randCross(founders, nEliteSires); pop$eliteSires[[2]]@sex[] <- "M"
Records <- recording(Records, mtfile, pop$eliteSires[[2]])

generation = 3
pop$eliteSires[[3]] <- randCross(founders, nEliteSires); pop$eliteSires[[3]]@sex[] <- "M"
Records <- recording(Records, mtfile, pop$eliteSires[[3]])

generation = 4
pop$eliteSires[[4]] <- randCross(founders, nEliteSires); pop$eliteSires[[4]]@sex[] <- "M"
pop$eliteDams[[1]] <- randCross(founders, 0.34*nEliteDams); pop$eliteDams[[1]]@sex[] <- "F" #nEliteDams*0.34
pop$commercial[[1]] <- randCross(founders, 0.34*nCommercial); pop$commercial[[1]]@sex[] <- "F" #nCommercial*0.34
Records <- recording(Records, mtfile, 
                     c(pop$eliteSires[[4]],
                                pop$eliteDams[[1]],
                                pop$commercial[[1]]))
Records$Lactation[Records$Generation == generation & Records$Sex == "F"] <- 5

generation = 5
pop$waitingBulls[[1]] <- randCross(founders, nWaitingBulls); pop$waitingBulls[[1]]@sex[] <- "M"
pop$eliteDams[[2]] <- randCross(founders, 0.49*nEliteDams); pop$eliteDams[[2]]@sex[] <- "F" #nEliteDams*0.49
pop$commercial[[2]] <- randCross(founders, 0.49*nCommercial); pop$commercial[[2]]@sex[] <- "F" #nCommercial*0.49
Records <- recording(Records, mtfile,
                     c(pop$waitingBulls[[1]],
                                pop$eliteDams[[2]],
                                pop$commercial[[2]]))
Records$Lactation[Records$Generation == generation & Records$Sex == "F"] <- 4

generation = 6
pop$waitingBulls[[2]] <- randCross(founders, nWaitingBulls); pop$waitingBulls[[2]]@sex[] <- "M"
pop$eliteDams[[3]] <- randCross(founders, (1-cullRate)*nEliteDams); pop$eliteDams[[3]]@sex[] <- "F" #nEliteDams*0.70
pop$commercial[[3]] <- randCross(founders, (1-cullRate)*nCommercial); pop$commercial[[3]]@sex[] <- "F" #nCommercial*0.70
Records <- recording(Records, mtfile, 
                     c(pop$waitingBulls[[2]],
                                pop$eliteDams[[3]],
                                pop$commercial[[3]]))
Records$Lactation[Records$Generation == generation & Records$Sex == "F"] <- 3

generation = 7
pop$waitingBulls[[3]] <- randCross(founders, nWaitingBulls); pop$waitingBulls[[3]]@sex[] <- "M"
pop$eliteDams[[4]] <- randCross(founders, nEliteDams); pop$eliteDams[[4]]@sex[] <- "F" #nEliteDams
pop$commercial[[4]] <- randCross(founders, nCommercial); pop$commercial[[4]]@sex[] <- "F" #nCommercial
Records <- recording(Records, mtfile,
                     c(pop$waitingBulls[[3]],
                                pop$eliteDams[[4]],
                                pop$commercial[[4]]))
Records$Lactation[Records$Generation == generation & Records$Sex == "F"] <- 2

generation = 8
pop$waitingBulls[[4]] <- randCross(founders, nWaitingBulls); pop$waitingBulls[[4]]@sex[] <- "M"
pop$heifers[[1]] <- randCross(founders, nHeifers); pop$heifers[[1]]@sex[] <- "F"
Records <- recording(Records, mtfile, 
                     c(pop$waitingBulls[[4]],
                                pop$heifers[[1]]))
Records$Lactation[Records$Generation == generation & Records$Sex == "F"] <- 1

generation = 9
pop$youngBulls[[1]] <- randCross(founders, nYoungBulls); pop$youngBulls[[1]]@sex[] <- "M"
pop$heifers[[2]] <- randCross(founders, nHeifers); pop$heifers[[2]]@sex[] <- "F"
Records <- recording(Records, mtfile, 
                     c(pop$youngBulls[[1]],
                                pop$heifers[[2]]))

generation = 10
pop$youngBulls[[2]] <- randCross(founders, nYoungBulls); pop$youngBulls[[2]]@sex[] <- "M"
pop$heifers[[3]] <- randCross(founders, nHeifers); pop$heifers[[3]]@sex[] <- "F"
Records <- recording(Records, mtfile, 
                     c(pop$youngBulls[[2]],
                                pop$heifers[[3]]))


# Estimate breeding values to start evaluation scenario
preparePAR(paste0(program, model))
runRENUM(Records, mt_ref, program, model)
Records = runBLUP(Records)
#file.remove("Blupf901.dat")

# ------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------ Run Burn-in generations -----------------------------------------------
for(year in 1:nBreeding){
  cat("Burn-in year:", year, "generation:", generation, "...\n")
  
  # define matting for each category:
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
  
  # Update year's records
  generation = generation + 1
  Records <- recording(Records, mtfile, 
                       c(pop$youngBulls[[3]], pop$heifers[[4]]))
  Covars <- runCovars(Covars, Records, 
                      c(pop$youngBulls[[3]], pop$heifers[[4]]))
  
  # Pre-test: select young males to enter progeny testing
  pop$youngBulls[[1]]@ebv <- as.matrix(with(Records, 
                                            nEbv[match(pop$youngBulls[[1]]@id, 
                                                       Records$IId)]))
  # assign breeding values to be able to select
  
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
  runRENUM(Records, mt_ref, program, model)
  Records = runBLUP(Records)
  
  pop$heifers[[1]]@ebv <- as.matrix(with(Records, 
                                         nEbv[match(pop$heifers[[1]]@id, 
                                                    Records$IId)]))
  
  # assign breeding values to be able to select
  
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
  Bias <- record.bias(Bias,
                      c(pop$heifers[[1]], dams, multiplier), Records, model)
  # estimate accuracies
  
  # moving categories
  for(j in 2:4){
    pop$eliteDams[[j]]@ebv <- as.matrix(with(Records, 
                                             nEbv[match(pop$eliteDams[[j]]@id, 
                                                        Records$IId)]))
    pop$commercial[[j]]@ebv <- as.matrix(with(Records, 
                                              nEbv[match(pop$commercial[[j]]@id, 
                                                         Records$IId)]))
    # assign ebvs
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
  
  rm(tmp); pop$heifers[[1]] <- NULL; pop$eliteDams[[1]] <- NULL; pop$commercial[[1]] <- NULL
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
                                              nEbv[match(pop$waitingBulls[[1]]@id, 
                                                         Records$IId)]))
  # assign ebv to be able to select
  
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
  
  pop$waitingBulls[[1]] <- NULL ; pop$eliteSires[[1]] <- NULL
  # non-selected candidates are culled. Sires used for 5 years are culled. 
  
}
# ------------------------------------------------------------------------------------------------------------------------
# ---------------------------------------- Reference population for Genetic models ---------------------------------------
# combine active population per category for genotyping
dams <- do.call(c, unlist(pop$eliteDams))
# All elite dams are genotyped (3165)
cows <- do.call(c, unlist(pop$commercial))
cows <- cows[cows@mother %in% dams@id]
# all elite dam daughters (phenotyped) are genotyped (55)
heifers <- c(pop$heifers[[1]], pop$heifers[[2]])
heifers <- heifers[heifers@mother %in% dams@id]

ref.pop <- c(dams, cows)
# reference pop: phenotyped + genotyped (3220)

# save genotyped femlaes ids
geno.id <- c(ref.pop@id, heifers@id) 

# all young bulls are genotyped
bulls <- c(pop$eliteSires[[4]], pop$eliteSires[[3]],
           do.call(c, unlist(pop$waitingBulls)),
           do.call(c, unlist(pop$youngBulls)))

# combine population to generate snp file
active <- c(ref.pop, bulls, heifers)
rm(ref.pop, bulls)

# ------------------------------------------------------------------------------------------------------------------------
# ----------------------------------- Correct scales for bias estimation & save data -------------------------------------

s <- c(sd(Records$nTbv[Records$Generation == (generation - nBreeding)]),
       sd(Records$mTbv[Records$Generation == (generation - nBreeding)]),
       sd(Records$tTbv[Records$Generation == (generation - nBreeding)]))
# estimate sd for first observation for each genome
Bias <- Bias %>% mutate(b0_n = b0_n/s[1],
                        b0_m = b0_m/s[2],
                        b0_t = b0_t/s[3])
# divide b0 results by sd

# save data
Records <- Records %>% filter(Generation >= generation-10)
Records <- Records %>% mutate(gv_corr = nTbv - mean(nTbv))
gv <- mean(Records$nTbv)














