# Evaluation of the impact of accounting for mitochondrial variation on breeding values estimation for dairy cattle
# Gabriela Mafra Fortuna
# Highlander Lab
# The Roslin Institute 
# July 2020 - updated Aug 2021

#library(AlphaSimR)
#library(tidyverse)
#library(Matrix)

# Simulation Parameters ---------------
# initiate datafiles
Records <- NULL
Covars <- NULL
Accuracy <- NULL
Bias <- NULL
geno.id <- NULL

# Global Parameters ---------------
# nuclear
nFem = 1000 
nQtl = 10^3
nSnp = 10^3
nChr = 10

# mito
mtInd = 5000
mNe = 1000
mu = 10*(2.5e-8)
histMtNe = c(1500, 2000, 2500, 3500, 7000, 10000, 17000, 62000)
histMtGen = c(25,  155,  455,  655, 1755,  2355,  3355, 33155)

# Trait Parameters ----------------
sd = 1890 # Kg
varP = sd^2
varA = varP*0.25
varM = varP*0.05
varPe = varP*0.1 
varE = varP*0.6
# lactation means (Kg)
mLac1 = 6733; mLac2 = 7440; mLac3 = 7344; mLac4 = 7482; mLac5 = 7168

# Population Parameters -----------
pop = list("heifers" = vector("list", 3), 
           "eliteDams" = vector("list", 4), 
           "commercial" = vector("list", 4),
           "youngBulls" = vector("list", 2), 
           "waitingBulls" = vector("list", 4), 
           "eliteSires" = vector("list", 4))

# selection rates
cullRate = 0.3 
concepRateH = 0.625
concepRateC = 0.35
#lnoConcep = 3

iYoungBulls = 0.4
iWaitingBulls = 0.1
iEliteDams = 0.25

# males
nEliteSires = 5
nWaitingBulls = nEliteSires/iWaitingBulls # 500
nYoungBulls = nWaitingBulls/iYoungBulls

# females
nHeifers = 25*(nWaitingBulls*4)
nReplacement = 790 # take all females from elite crossing
# initial number of heifers + the amount of females culled yearly
nEliteDams = iEliteDams*nHeifers
nCommercial = (1 - (iEliteDams+cullRate))*nHeifers

# Herd pars
#sHerd = 100
#nHerd = nFem/sHerd

# Functions ---------------

# Data Recording
recording <- function(datafile, mtdna, pop){
  
  if(generation == 0){
    id = sample(mtdna$ML, nFem, replace = TRUE)
    mt = with(mtdna, mTbv[match(id, mtdna$ML)])
    id = c(sapply(id, function(id) c(rep(NA, 1), id)))
    mt = c(sapply(mt, function(mt) c(rep(0, 1), mt)))
    }
  
  
  datafile = rbind(datafile,
                   tibble(Generation = generation,
                          Program    = program,
                          Model      = model,
                          Selection  = selection,
                          IId        = pop@id,
                          Sex        = pop@sex,
                          FId        = pop@father,
                          MId        = pop@mother,
                          ML         = if(generation==0){id}else{with(datafile, ML[match(pop@mother, datafile$IId)])},
                          nTbv       = pop@gv[,1],
                          mTbv       = if(generation==0){mt}else{with(datafile, mTbv[match(pop@mother, datafile$IId)])},
                          tTbv       = sum(nTbv, mTbv),
                          nEbv       = 0,
                          mEbv       = 0,
                          tEbv       = 0,
                          pe         = ifelse(pop@sex=="F", rnorm(pop@nInd)*sqrt(varPe), NA),
                          Lactation  = NA,
                          Pheno1     = ifelse(pop@sex=="F",
                                              mLac1 + nTbv + mTbv + 
                                                pe + rnorm(pop@nInd)*sqrt(varE), NA),
                          Pheno2     = ifelse(pop@sex=="F",
                                              mLac2 + nTbv + mTbv + 
                                                pe + rnorm(pop@nInd)*sqrt(varE), NA),
                          Pheno3     = ifelse(pop@sex=="F",
                                              mLac3 + nTbv + mTbv + 
                                                pe + rnorm(pop@nInd)*sqrt(varE), NA),
                          Pheno4     = ifelse(pop@sex=="F",
                                              mLac4 + nTbv + mTbv + 
                                                pe + rnorm(pop@nInd)*sqrt(varE), NA),
                          Pheno5     = ifelse(pop@sex=="F",
                                              mLac5 + nTbv + mTbv + 
                                                pe + rnorm(pop@nInd)*sqrt(varE), NA),  

                          gv_corr    = ifelse(Generation < 32, nTbv, 0)))
}

record.accuracy = function(datafile, pop, recfile, model){
  recfile <- recfile %>% filter(IId %in% pop@id)
  
  # creates mtDNA pop for calculating genic var
  y <- recfile %>% select(IId, ML)
  mtDNApop=NULL
  while(nrow(y) > 0){
    p <- distinct(y, ML)
    mtDNApop = c(mtDNApop, selectInd(mtDNAx, nInd=nrow(p), simParam = if(traitScen=="maxQTL"){SP2}else{SP3}, 
                                     use="rand", candidates = p$ML))
    y <- y %>% group_by(ML) %>% filter(!(row_number() == 1))
  }
  mtDNApop <- mergePops(mtDNApop)
  
  
  datafile <- rbind(datafile,
                    tibble(Generation = generation,
                           Breeding_year = year,
                           SelGroup   = category,  
                           Program    = program,
                           Model      = model,
                           Selop      = selection,
                           nInd       = pop@nInd,
                           ids        = list(pop@id),
                           # estimated bv
                           nEbv       = list(with(recfile, nEbv[match(pop@id, recfile$IId)])),
                           mEbv       = list(with(recfile, mEbv[match(pop@id, recfile$IId)])),
                           tEbv       = list(with(recfile, tEbv[match(pop@id, recfile$IId)])),
                           # true bv
                           nTbv       = list(with(recfile, nTbv[match(pop@id, recfile$IId)])),
                           mTbv       = list(with(recfile, mTbv[match(pop@id, recfile$IId)])),
                           tTbv       = list(with(recfile, tTbv[match(pop@id, recfile$IId)])),
                           # gen & var trends
                           nVar        = varG(pop), # genetic var
                           genicVarN   = genicVarG(pop, SP), # genic var
                           mVar        = var(recfile$mTbv)*((pop@nInd-1)/pop@nInd), 
                           genicVarM   = genicVarG(mtDNApop, simParam = if(traitScen=="maxQTL"){SP2}else{SP3}), 
                           tVar        = var(recfile$tTbv)*((pop@nInd-1)/pop@nInd),
                           genicVarT   = sum(genicVarN, genicVarM),
                           nMean       = mean(recfile$nTbv),
                           mMean       = mean(recfile$mTbv),
                           tMean       = mean(recfile$tTbv)
                    ))
}

record.bias = function(datafile, pop, recfile, model){
  recfile <- recfile %>% filter(IId %in% pop@id)
  datafile <- rbind(datafile,
                    tibble(Generation = generation,
                           Breeding_year = year,
                           SelGroup   = category,
                           nInd       = pop@nInd,
                           Program    = program,
                           Model      = model,
                           Selop      = selection,
                           # 1: nTbv vs nEbv - all models
                           b0_n       = lm(gv_corr~nEbv, recfile)$coefficients[[1]],
                           b1_n       = lm(gv_corr~nEbv, recfile)$coefficients[[2]],
                           # 2: mTbv vs mEbv - mt models
                           b0_m       = ifelse(Model == "mt",
                                               lm(mTbv~mEbv, recfile)$coefficients[[1]], NA),
                           b1_m       = ifelse(Model == "mt",
                                               lm(mTbv~mEbv, recfile)$coefficients[[2]], NA),
                           # 3: tTbv vs tEbv - all models
                           b0_t       = lm((mTbv + gv_corr)~tEbv, recfile)$coefficients[[1]],
                           b1_t       = lm((mTbv + gv_corr)~tEbv, recfile)$coefficients[[2]],
                    ))
}

runCovars <- function(datafile, recfile, pop){
  recfile <- recfile %>% filter(IId %in% pop@id)
  
  # creates mtDNA pop for calculating genic var
  y <- recfile %>% select(IId, ML)
  mtDNApop=NULL
  while(nrow(y) > 0){
    p <- distinct(y, ML)
    mtDNApop = c(mtDNApop, selectInd(mtDNAx, nInd=nrow(p), simParam = if(traitScen=="maxQTL"){SP2}else{SP3},
                                     use="rand", candidates = p$ML))
    y <- y %>% group_by(ML) %>% filter(!(row_number() == 1))
  }
  mtDNApop <- mergePops(mtDNApop)
  
  datafile <- rbind(datafile,
                    tibble(Generation  = generation,
                           Breeding_year = year,
                           nInd       = pop@nInd,
                           Program    = program,
                           Model      = model,
                           Selop      = selection,
                           result      = cor(recfile$nTbv, recfile$mTbv),
                           nVar        = var(recfile$nTbv)*((nrow(recfile)-1)/nrow(recfile)),
                           genicVarN   = genicVarG(pop, SP),
                           mVar        = var(recfile$mTbv)*((nrow(recfile)-1)/nrow(recfile)),
                           genicVarM   = genicVarG(mtDNApop, simParam = if(traitScen=="maxQTL"){SP2}else{SP3}),
                           tVar        = var(recfile$tTbv)*((nrow(recfile)-1)/nrow(recfile)),
                           genicVarT   = sum(genicVarN, genicVarM),
                           nMean       = mean(recfile$nTbv),
                           mMean       = mean(recfile$mTbv),
                           tMean       = mean(recfile$tTbv),
                           ML          = sum(!duplicated(recfile$ML) & !is.na(recfile$ML))
                    ))
} 

# RENUMF90 function
runRENUM = function(datafile, mtdna_ids, program, model){
  # create pedigree file
  pedFile <- datafile[, c("IId", "FId", "MId")]
  write.table(pedFile, "Blupf90.ped", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ", na="0")
  rm(pedFile)
  
  # create phenotype file
  phenoFile = datafile %>%
    filter(Sex == "F", !(is.na(Lactation)))
  # only females with at least one lactation
  phenoFile$Pheno <- ifelse(phenoFile$Lactation == 1, phenoFile$Pheno1,
                            ifelse(phenoFile$Lactation == 2, phenoFile$Pheno2,
                                   ifelse(phenoFile$Lactation == 3, phenoFile$Pheno3,
                                          ifelse(phenoFile$Lactation == 4, phenoFile$Pheno4, phenoFile$Pheno5))))
  
  phenoFile <- phenoFile[, c("IId", "Pheno", "Lactation", "ML")]
  
  # Join mtDNA cross-ref file with data file
  # Check if inner_join is correct and if it deals with repetitions, etc...
  pheno_export = inner_join(phenoFile, mtdna_ids, by = c("ML" = "mt_id"))
 
  write.table(pheno_export, "Blupf901.dat", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ", na = "0", append = TRUE)
  rm(phenoFile, pheno_export)
  
  # insert variance components into parameter file
  system(paste0('sed "s/varE/', varE, '/g" renumf900.par > renumf901.par'))
  file.remove("renumf900.par")
  system(paste0('sed "s/varA/', varA, '/g" renumf901.par > renumf902.par'))
  file.remove("renumf901.par")
  system(paste0('sed "s/varPe/', varPe, '/g" renumf902.par > renumf90.par'))
  file.remove("renumf902.par")
  
  # call RENUMF90
  system(command = "echo renumf90.par | $HOME/bin/renumf90 | tee renum.log")
  
  if(model=="mt"){
    # On line 16 - add new effect and levels
    nML <- as.numeric(nrow(mtdna_ids))
    system(command="awk 'NR==16 { print \" 6        nML cross\";}1' renf90.par > tmp1.par")
    system(paste0('sed "s/nML/', nML,'/g" tmp1.par > tmp.par '))
    file.remove("tmp1.par", "renf90.par")
    
    if(program=="GEN"){system(command = "sh mtdnarenf90.sh")}else{system(command="sh pedmtdnarenf90.sh")}
    }
}

# AIREMLF90 function
varComp <- function(model){
  # replace variances on renf90.par (only standard models)
  system(command = "sed -i '' '17s/'.*'/10000/' renf90.par")
  system(command = "sed -i '' '25s/'.*'/10000/' renf90.par")
  system(command = "sed -i '' '33s/'.*'/10000/' renf90.par")
  
  system(command = "echo renf90.par | $HOME/bin/airemlf90 | tee airemlf90.log")
  
  # return new variances
  # varA
  system(command = "sed -n '6p' airemlf90.log | tee vars.txt")
  varA <<- read_table2("vars.txt", col_names = FALSE)[[1]]
  # varPe
  system(command = "sed -n '8p' airemlf90.log | tee vars.txt")
  varPe <<- read_table2("vars.txt", col_names = FALSE)[[1]]
  
  if(model == "mt"){
    # varM
    system(command = "sed -n '10p' airemlf90.log | tee vars.txt")
    varM <<- read.table2("vars.txt", col_names = FALSE)[[1]]
    
    # varE
    system(command = "sed -n '12p' airemlf90.log | tee vars.txt")
    varE <<- read_table2("vars.txt", col_names = FALSE)[[1]]
  }else{
    # varE
    system(command = "sed -n '10p' airemlf90.log | tee vars.txt")
    varE <<- read_table2("vars.txt", col_names = FALSE)[[1]]
  }
  
  # replace variances on renf90.par (only standard models)
  system(paste0('sed "17s/.*/',varE, '/g" renf90.par > renf901.par'))
  file.remove("renf90.par")
  system(paste0('sed "25s/.*/',varA, '/g" renf901.par > renf902.par'))
  file.remove("renf901.par")
  system(paste0('sed "33s/.*/',varPe, '/g" renf902.par > renf90.par'))
  file.remove("renf902.par")
}
 
# run BLUPF90 function
runBLUP = function(datafile, mtdna_ids = NULL){
  system(command = "echo renf90.par | $HOME/bin/blupf90 | tee blup.log")
  # run BLUPF90
  
  # recover solutions
  sol  <- read_table2("solutions", col_names=FALSE, skip = 1,
                      col_types = cols(.default = col_double(),
                                       X1 = col_double(),
                                       X2 = col_double(),
                                       X3 = col_double(),
                                       X4 = col_double()))
  colnames(sol) = c("Trait", "Effect", "Level", "Solution")
  
  # Extract nEBV from file
  nebv = sol %>%
    filter(Trait == 1 & Effect == 2) %>% # renadd03 -> 3
    select("Level", "Solution")
  
  renadd <- read.table("renadd02.ped")
  renadd <- renadd[order(renadd$V1), c(1,10)]
  nebv <- merge(renadd, nebv, by.x="V1", by.y="Level")
  
  # Update nuclear breeding values in database:
  nebv <- as_tibble(with(nebv, Solution[match(datafile$IId, nebv$V10)]))
  nebv <- nebv %>% replace(is.na(.), 0)
  datafile$nEbv <- nebv$value 
  
  # Extract mEBV from file (only mt models)
  if(is.null(mtdna_ids)){
    datafile <- datafile %>% mutate(tEbv = nEbv)
    }else{
    mebv = sol %>% 
      filter(Trait == 1 & Effect == 4) %>%
      select("Level", "Solution")
    
    # "mebv" has results for 100 mtdnas - RETRIEVE CORRECT IDS TO MATCH RESULTS
    
    # retrieve level coding - solutions with new and original ids
    mtresult <- inner_join(mebv, mtdna_ids, by = c("Level" = "mt_new_id"))
    
    # Update mito breeding values in database: compare col1 with ML (X1 = original ids)
    mebv <- as_tibble(with(mebv, Solution[match(datafile$ML, mtresult$mt_id)]))
    mebv <- mebv %>% replace(is.na(.), 0)
    datafile$mEbv <- mebv$value 
    
    datafile <- datafile %>% mutate(tEbv = nEbv + mEbv)
  }

  #  datafile$tEbv <- sum(datafile$nEbv, datafile$mEbv)
  
  return(datafile)
}

# generate snp file
snpData <- function(pop){
  cat("Generating SNP file...\n")
  id <- pop@id
  pr <- pullSnpGeno(pop)
  prt <- apply(pr[,1:ncol(pr)], MARGIN = 1, FUN = paste, sep = "", collapse= "")
  snp.file <- cbind(id, prt)
  write.table(snp.file, "mrk.tmp", col.names=FALSE, row.names=FALSE, quote=FALSE, sep=" ", append = TRUE)
  rm(snp.file, id, pr, prt)
  system(command = "sh prepgeno.sh")
  #system(command = "awk '{printf(\"%-8s %s\\n\", $1,i$2)}' mrk.tmp > snp.dat")
}

sink("prepgeno.sh", type="output")
writeLines("#!/bin/bash
# Prepare Marker Genotypes file for blupf90 format
# And automatically cut old markers for blupf90 limit (25k)

# Blupf90 fixed format
awk '{printf(\"%-8s %s\\n\", $1,i$2)}' mrk.tmp > mrk2.tmp

gen=$(< mrk2.tmp wc -l)
# hard limit = 25000
limit=10200 # value does not make sense here but increasing is a problem to run
# *change when running on server!!

if [ \"$gen\" -gt \"$limit\" ]
 then
   cold=(`expr $gen - $limit + 1`)
   tail -n +$cold mrk2.tmp > snp.dat
 else
   cp mrk2.tmp snp.dat
fi
")
sink()

# modify renf90.par to include mitochondrial effect
sink("mtdnarenf90.sh", type="output")
writeLines("#!/bin/bash
# Correct renf90.par for mtdna
# first - overwrite line 7 (3 --> 4)           
awk 'NR==7 {$0=\"    4\"}1' tmp.par > tmp1.par
rm tmp.par

# Add block text to line 35
awk 'NR==35 { print \" RANDOM_GROUP|END|4|END|RANDOM_TYPE|END|user_file|END|FILE|END|mtdnaGinv.txt|END|(CO)VARIANCES|END|178605 \" ;}1' tmp1.par > tmp2.par
rm tmp1.par

# Break into different lines
sed 's/|END|/\\n/g' tmp2.par > renf90.par
rm tmp2.par
")
sink()

sink("pedmtdnarenf90.sh", type="output")
writeLines("#!/bin/bash
# Correct renf90.par for mtdna
# first - overwrite line 7 (3 --> 4)           
awk 'NR==7 {$0=\"    4\"}1' tmp.par > tmp1.par
rm tmp.par

# Add block text to line 35
awk 'NR==35 { print \" RANDOM_GROUP|END|4|END|RANDOM_TYPE|END|diagonal|END|FILE|END|            |END|(CO)VARIANCES|END|178605 \" ;}1' tmp1.par > tmp2.par
rm tmp1.par

# Break into different lines
sed 's/|END|/\\n/g' tmp2.par > renf90.par
rm tmp2.par
")
sink()

# Create Initial mtDNA Ginv matrix
mtdnaGinv <- function(mtdnaPop, simParam){
  M = pullSnpGeno(pop=mtdnaPop, simParam=simParam)
  p = colMeans(M)
  
  P = matrix(p, nrow = nrow(M), ncol = ncol(M), byrow = TRUE)
  Z = M-P
  k = (sum(p*(1-p)))
  
  ZZ = Z%*%t(Z)
  G = ZZ/k
  
  nid = nrow(G)
  Gstar = (0.99*G) + (diag(0.01, nid))
  
  Ginv = solve(Gstar)
  
  # Save inverted matrix in BLUPF90 readable format 
  G31 = as(forceSymmetric(Ginv), "dsTMatrix")
  #G31 = as(Ginv, "dsTMatrix")
  G32 = summary(G31)
  write.table(G32, file = "mtdnaGinv.txt", row.names = FALSE, col.names = FALSE)
  
  # Create cross-ref file based on the exported rownames
  # this will deal with the correction of ML_ids during BLUP estimations
  mt_ref <<- tibble(mt_id = rownames(Z), mt_new_id = seq(1:length(mt_id)))
}


# prepare BLUPF90 parameters according to model running
preparePAR <- function(model = c("PEDstd", "PEDmt", "GENstd", "GENmt")){
  if(model == "PEDstd" | model == "PEDmt"){
    sink("renumf900.par", type="output")
    writeLines("#renumf90 parametar file
DATAFILE
Blupf901.dat
TRAITS
2
FIELDS_PASSED TO OUTPUT
1 4 5
#original_id mtdna_original_id mtdna_order_covariances
WEIGHT(S)
  
RESIDUAL_VARIANCE
varE
EFFECT
3 cross alpha
#lactation order
EFFECT
1 cross alpha
#animal
RANDOM
animal
OPTIONAL
pe
FILE
Blupf90.ped
PED_DEPTH
0
INBREEDING
pedigree
(CO)VARIANCES
varA
(CO)VARIANCES_PE
varPe
OPTION use_yams
")
    sink()
  }else{
    sink("renumf900.par", type="output")
    writeLines("#renumf90 parametar file
DATAFILE
Blupf901.dat
TRAITS
2
FIELDS_PASSED TO OUTPUT
1 4 5
#original_id mtdna_original_id mtdna_order_covariances
WEIGHT(S)

RESIDUAL_VARIANCE
varE
EFFECT
3 cross alpha
#lactation order
EFFECT
1 cross alpha
#animal
RANDOM
animal
OPTIONAL
pe
FILE
Blupf90.ped
SNP_FILE
snp.dat
PED_DEPTH
0
INBREEDING
pedigree
(CO)VARIANCES
varA
(CO)VARIANCES_PE
varPe
OPTION use_yams
OPTION no_quality_control
OPTION thrStopCorAG 0
")
    sink()
  }
}
