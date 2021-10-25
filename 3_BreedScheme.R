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

for(year in (year+1):(year+nBreeding)){
  cat("Breeding year:", year, "generation:", generation, "...\n")
  
  # define matting groups for each category:
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
    cat("Generating updated snp file...\n")
    snpData(c(pop$youngBulls[[3]], pop$heifers[[4]][pop$heifers[[4]]@mother %in% dams@id]))
    
    geno.id <- c(geno.id, 
                 pop$heifers[[4]][pop$heifers[[4]]@mother %in% dams@id]@id)
  }
  
  # Update year's records
  generation = generation + 1
  Records <- recording(Records, mtfile, 
                       c(pop$youngBulls[[3]], pop$heifers[[4]]))
  
  # Pre-test: select young males to enter progeny testing
  pop$youngBulls[[1]]@ebv <- add.ebv(pop$youngBulls[[1]], Records, selection)
  # assign breeding values to be able to select
  
  category = "young_bulls"
  catSummary <- summarise.category(catSummary, pop$youngBulls[[1]], Records, mtDNAx)
  youngbulls <- pop$youngBulls[[1]]
  # record summary for category
  
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
  Records = runBLUP(Records, if(model == "mt"){mtdna_ids= mt_ref})
  
  genSummary <- summarise.generation(genSummary, c(pop$youngBulls[[2]], pop$heifers[[4]]), Records, mtDNAx)
  
  pop$heifers[[1]]@ebv <- add.ebv(pop$heifers[[1]], Records, selection)
  # assign breeding values to be able to select
  
  categories <- if(program == "GEN"){
    c("heifers_geno", "heifers_no_geno", "1st_lact_geno", "1st_lact_no_geno", "cows_geno", "cows_no_geno", "proven_bulls")
  }else{c("heifers", "1st_lact", "cows", "proven_bulls")}
  
  agroups <- if(program == "GEN"){
    list(pop$heifers[[2]][pop$heifers[[2]]@id %in% geno.id], pop$heifers[[2]][!(pop$heifers[[2]]@id %in% geno.id)],
         pop$heifers[[1]][pop$heifers[[1]]@id %in% geno.id], pop$heifers[[1]][!(pop$heifers[[1]]@id %in% geno.id)],
         c(dams, multiplier)[c(dams, multiplier)@id %in% geno.id], c(dams, multiplier)[!(c(dams, multiplier)@id %in% geno.id)],
         pop$waitingBulls[[1]])
  }else{list(pop$heifers[[2]], pop$heifers[[1]], c(dams, multiplier), pop$waitingBulls[[1]])}
  
  for(i in 1:length(agroups)){
    if(agroups[[i]]@nInd !=0){
      category = categories[i]
      catSummary <- summarise.category(catSummary, agroups[[i]], Records, mtDNAx)
    }
  }
  # record summary for category
  
  # moving categories
  for(j in 2:4){
    # assign ebvs
    pop$eliteDams[[j]]@ebv <- add.ebv(pop$eliteDams[[j]], Records, selection)
    pop$commercial[[j]]@ebv <- add.ebv(pop$commercial[[j]], Records, selection)
    
    # select best females
    pop$eliteDams[[j]] <- selectInd(pop$eliteDams[[j]],
                                    nInd = (1-cullRate)*pop$eliteDams[[j]]@nInd,
                                    use="ebv")
    pop$commercial[[j]] <- selectInd(pop$commercial[[j]],
                                     nInd = (1-cullRate)*pop$commercial[[j]]@nInd, 
                                     use="ebv")
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
  pop$waitingBulls[[1]]@ebv <- add.ebv(pop$waitingBulls[[1]], Records, selection)
  # assign ebv to be able to select
  
  pop$eliteSires <- c(pop$eliteSires,
                      selectInd(pop$waitingBulls[[1]], nInd = nEliteSires, 
                                use="ebv"))
  # select best performing bulls to become elite sires
  
  pop$waitingBulls[[1]] <- NULL ; pop$eliteSires[[1]] <- NULL
  # non-selected candidates are culled. Sires used for 5 years are culled. 
  
}
# ------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------ Correct scales for bias estimation ------------------------------------------ 

s <- c(sd(Records$nTbv[Records$Generation == (generation - nBreeding)]),
       sd(Records$mTbv[Records$Generation == (generation - nBreeding)]),
       sd(Records$tTbv[Records$Generation == (generation - nBreeding)]))
# estimate sd for first observation for each genome
#Bias <- Bias %>% mutate(b0_n = ifelse(Generation >= 32, b0_n/s[1], b0_n),
#                        b0_m = ifelse(Generation >= 32, b0_m/s[2], b0_m),
#                        b0_t = ifelse(Generation >= 32, b0_t/s[3], b0_t))
# divide b0 results by sd
