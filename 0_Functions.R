# mtDNA Project: Should mitochondrial variation be accounted for on BVE
# Dairy Cattle Breeding Simulations
# Roslin Institute
# from July 2020
# Definition of genetic and populational parameters and functions to be used thrughout simulations

# NOTE: Change code from macos to linux when running on server (line 158/9)

# Global Parameters
# Genetic parameters
# nuclear genome parameters
nInd = 1000
nQtl = 10^3
nSnp = 10^3
nChr = 10

# mitogenome parameters
mtInd = 5000
mNe = 1000
mu = 10*(2.5e-8)
histMtNe = c(1500, 2000, 2500, 3500, 7000, 10000, 17000, 62000)
histMtGen = c(25,  155,  455,  655, 1755,  2355,  3355, 33155)

# Population parameters
# population sizes
nYB = 418; nPTB = 50; nES = 5 
# number of young Bulls, waiting-bulls, elite Sires
nFM = 7110; nCW = 4802; nED = 250
# number of young Females, cows, elite Dams

# lactation means
mLac1 = 6733; mLac2 = 7440; mLac3 = 7344; mLac4 = 7482; mLac5 = 7168
# sd for MY from SPEHAR et al., 2017 (1890) -> varP = sd^2 = 3572100
varP = 3572100
varA = varP*0.25
varM = varP*0.05
varE = varP*0.6
varPe = varP*0.1
#---------------------------------------------------------------------------------------
# Required functions

# Recording function
recording <- function(datafile, pop, mtdna){
  
  if(generation == 0){id = sample(mtdna$ML, nInd, replace = TRUE)
  mt = with(mtdna, mTbv[match(id, mtdna$ML)])
  id = c(sapply(id, function(id) c(rep(NA, 1), id)))
  mt = c(sapply(mt, function(mt) c(rep(0, 1), mt)))}
  # deals with assigning ML only to females
  
  datafile = rbind(datafile, 
                   tibble(Generation = generation,
                          Program    = program,
                          Model      = model,
                          IId        = pop@id,
                          FId        = pop@father,
                          MId        = pop@mother,
                          Gender     = pop@sex,
                          nTbv       = pop@gv[,1],
                          nEbv       = 0, # comes form BLUP estimation (all models)
                          ML         = if(generation==0){id}else{with(datafile, ML[match(pop@mother, datafile$IId)])},
                          mTbv       = if(generation==0){mt}else{with(datafile, mTbv[match(pop@mother, datafile$IId)])},
                          mEbv       = 0, # comes from BLUP estimation (only mt models)
                          tTbv       = nTbv + mTbv, # (all models)
                          tEbv       = 0, # = nEBV | nEBV + mEBV
                          pe         = ifelse(pop@sex=="F",
                                              rnorm(pop@nInd)*sqrt(varPe), NA),
                          Lactation  = NA,
                          PhenoL1    = ifelse(pop@sex=="F",
                                              mLac1 + nTbv + mTbv + pe + rnorm(pop@nInd)*sqrt(varE), NA),
                          PhenoL2    = ifelse(pop@sex=="F",
                                              mLac2 + nTbv + mTbv + pe + rnorm(pop@nInd)*sqrt(varE), NA),
                          PhenoL3    = ifelse(pop@sex=="F",
                                              mLac3 + nTbv + mTbv + pe + rnorm(pop@nInd)*sqrt(varE), NA),
                          PhenoL4    = ifelse(pop@sex=="F",
                                              mLac4 + nTbv + mTbv + pe + rnorm(pop@nInd)*sqrt(varE), NA),
                          PhenoL5    = ifelse(pop@sex=="F",
                                              mLac5 + nTbv + mTbv + pe + rnorm(pop@nInd)*sqrt(varE), NA),
                          gv_corr    = ifelse(Generation < 32, nTbv, 0)))
}

record.accuracy = function(datafile, pop, recfile, model){
  recfile <- recfile %>% filter(IId %in% pop@id)
#  if(generation > 31){recfile <- recfile %>% mutate(nTbv = nTbv-gv)}
  
  datafile <- rbind(datafile,
                   tibble(Generation = generation,
                          SelGroup   = category,
                          Program    = program,
                          Model      = model,
                          # 1: nTbv vs nEbv - all models
                          acc_n      = cor(recfile$nTbv, recfile$nEbv),
                          r2_n       = acc_n^2,
                          b0_n       = lm(gv_corr~nEbv, recfile)$coefficients[[1]],
                          b1_n       = lm(gv_corr~nEbv, recfile)$coefficients[[2]],
                          # 2: mTbv vs mEbv - mt models
                          acc_m      = ifelse(Model == "mt",
                                              cor(recfile$mTbv, recfile$mEbv), NA),
                          r2_m       = ifelse(Model == "mt",
                                              acc_m^2, NA),
                          b0_m       = ifelse(Model == "mt",
                                              lm(mTbv~mEbv, recfile)$coefficients[[1]], NA),
                          b1_m       = ifelse(Model == "mt",
                                              lm(mTbv~mEbv, recfile)$coefficients[[2]], NA),
                          # 3: tTbv vs tEbv - all models
                          acc_t      = cor(recfile$tTbv, recfile$tEbv),
                          r2_t       = acc_t^2,
                          b0_t       = lm((mTbv + gv_corr)~tEbv, recfile)$coefficients[[1]],
                          b1_t       = lm(tTbv~tEbv, recfile)$coefficients[[2]]))
                          

}

runCovars <- function(datafile, recfile){
  recfile <- recfile %>% filter(Generation == generation)
  
  datafile <- rbind(datafile,
                    tibble(Generation  = generation,
                           result      = cov(recfile$nTbv, recfile$mTbv),
                           mtVar       = var(recfile$mTbv),
                           nVar        = var(recfile$nTbv),
                           nMean       = mean(recfile$nTbv),
                           mMean       = mean(recfile$mTbv),
                           ML          = sum(!duplicated(recfile$ML) & !is.na(recfile$ML))
                           ))
} 
# REMEMBER: divide vars by first observation
  

# run BLUPF90 function
runBLUP = function(datafile){
  # create pedigree file
  pedFile <- datafile[, c("IId", "FId", "MId")]
  write.table(pedFile, "Blupf90.ped", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ", na="0")
  rm(pedFile)

  # create phenotype file
  phenoFile = datafile %>%
    filter(Gender == "F", !(is.na(Lactation)))
    # only females with at least one lactation
  phenoFile$Pheno <- ifelse(phenoFile$Lactation == 1, phenoFile$PhenoL1,
                                   ifelse(phenoFile$Lactation == 2, phenoFile$PhenoL2,
                                          ifelse(phenoFile$Lactation == 3, phenoFile$PhenoL3,
                                                 ifelse(phenoFile$Lactation == 4, phenoFile$PhenoL4, phenoFile$PhenoL5))))
  
  phenoFile <- phenoFile[, c("IId", "Pheno", "Lactation", "IId", "ML")]
  write.table(phenoFile, "Blupf901.dat", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ", na = "0", append = TRUE)
  rm(phenoFile)
  
  # estimations
  system(command = "echo renumf90.par | $HOME/bin/renumf90 | tee renum.log")
  # run RENUMF90
  if("GEN" %in% datafile$Program & "mt" %in% datafile$Model){
    #system(command = "sed -i '' '40i\\'$'\\n''mtdnaGinv.txt' renf90.par") # macos
    system(command = "sed -i '40i''mtdnaGinv.txt' renf90.par") # linux
    system(command = "sed -i '41d' renf90.par") # remove extra line  
    
    # insert Ginv file name into correct place
    system(command="awk '/diagonal/{count++;if(count==2){sub(\"diagonal\",\"user_file\")}}1' renf90.par > renf901.par")
    # replace second occurrence of diagonal to user_file
    file.remove("renf90.par")
    file.rename("renf901.par", "renf90.par")
  }
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
  if(4 %in% sol$Effect){
    mebv = sol %>% 
      filter(Trait == 1 & Effect == 4) %>%
      select("Level", "Solution")
    
    # retrieve level coding
    system(command = "awk '/Effect group 4/ {matched = 1} matched' renf90.tables > newfile")
    renum <- read_table2("newfile", col_names=FALSE, skip = 2,
                         col_types = cols(.default = col_double(),
                                          X1 = col_double(),
                                          X2 = col_double(),
                                          X3 = col_double()))
    
    mebv <- merge(renum, mebv, by.x="X3", by.y="Level")
    # level = renumbered values (ids)

    # Update mito breeding values in database: compare col1 with ML (X1 = original ids)
    mebv <- as_tibble(with(mebv, Solution[match(datafile$ML, mebv$X1)]))
    mebv <- mebv %>% replace(is.na(.), 0)
    datafile$mEbv <- mebv$value
  }
  
  datafile <- datafile %>% mutate(tEbv = nEbv + mEbv)
  #datafile$tEbv["mt" %in% datafile$Model] = datafile$nEbv + datafile$mEbv
  #datafile$tEbv["std" %in% datafile$Model] = datafile$nEbv
  
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
}


# prepare BLUPF90 parameters according to model running
preparePAR <- function(model = c("PED", "PEDmt", "GEN", "GENmt")){
  if(model == "PED"){
    sink("renumf90.par", type="output")
    writeLines("#renumf90 parametar file
DATAFILE
Blupf901.dat
TRAITS
2
FIELDS_PASSED TO OUTPUT
1 3 4
#id lactation pe
WEIGHT(S)
  
RESIDUAL_VARIANCE
0.6
EFFECT
3 cross alpha
#lactation order
EFFECT
1 cross alpha
#animal
RANDOM
animal
FILE
Blupf90.ped
PED_DEPTH
0
INBREEDING
pedigree
(CO)VARIANCES
0.25
EFFECT
4 cross alpha
#pe
RANDOM
diagonal
(CO)VARIANCES
0.10
")
sink()
  }else if(model=="PEDmt"){
    sink("renumf90.par", type="output")
    writeLines("#renumf90 parametar file
DATAFILE
Blupf901.dat
TRAITS
2
FIELDS_PASSED TO OUTPUT
1 3 4 5
#id lactation pe mt
WEIGHT(S)

RESIDUAL_VARIANCE
0.6
EFFECT
3 cross alpha
#lactation order
EFFECT
1 cross alpha
#animal
RANDOM
animal
FILE
Blupf90.ped
PED_DEPTH
0
INBREEDING
pedigree
(CO)VARIANCES
0.25
EFFECT
4 cross alpha
#pe
RANDOM
diagonal
(CO)VARIANCES
0.10
EFFECT
5 cross alpha
#ml
RANDOM
diagonal
(CO)VARIANCES
0.05
")
sink()
  }else if(model=="GEN"){
    sink("renumf90.par", type="output")
    writeLines("#renumf90 parametar file
DATAFILE
Blupf901.dat
TRAITS
2
FIELDS_PASSED TO OUTPUT
1 3 4
#id lactation pe
WEIGHT(S)

RESIDUAL_VARIANCE
0.6
EFFECT
3 cross alpha
#lactation order
EFFECT
1 cross alpha
#animal
RANDOM
animal
FILE
Blupf90.ped
SNP_FILE
snp.dat
PED_DEPTH
0
INBREEDING
pedigree
(CO)VARIANCES
0.25
EFFECT
4 cross alpha
#pe
RANDOM
diagonal
(CO)VARIANCES
0.10
OPTION no_quality_control
OPTION thrStopCorAG 0
")
sink()
  }else{
    sink("renumf90.par", type="output")
    writeLines("#renumf90 parametar file
DATAFILE
Blupf901.dat
TRAITS
2
FIELDS_PASSED TO OUTPUT
1 3 4 5
#id lactation pe ml
WEIGHT(S)

RESIDUAL_VARIANCE
0.6
EFFECT
3 cross alpha
#lactation order
EFFECT
1 cross alpha
#animal
RANDOM
animal
FILE
Blupf90.ped
SNP_FILE
snp.dat
PED_DEPTH
0
INBREEDING
pedigree
(CO)VARIANCES
0.25
EFFECT
4 cross alpha
#pe
RANDOM
diagonal
(CO)VARIANCES
0.10
EFFECT
5 cross alpha
#ml
RANDOM
diagonal
(CO)VARIANCES
0.05
OPTION no_quality_control
OPTION thrStopCorAG 0
")
sink()
 }
  
}







