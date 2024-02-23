# Accounting for nuclear - and mito-genome in dairy cattle breeding: a simulation study
# The Roslin Institute
# Since July 2020
# Execution script: run burn-in plus 4 evaluation scenarios with replicates

# Simulating four breeding scenarios, two trait scenarios (mt causal loci)
# - breeding scenarios: PBLUP | mtPBLUP | GBLUP | mtGBLUP
# - trait scenarios: all loci are causal | 1 locous is causal
# We are using the same breeding scheme but selecting on pedigree-based or genome-based
# estimated breeding values 

rm(list=ls())
# load packages
library(AlphaSimR)
library(tidyverse)
library(Matrix)

# Define replicates and trait scenario: traitScen
#     maxQTL -> all segSites are causal
#     minQTL -> 1 segSite is causal

scriptsDir = "/home/v1gfortu/scripts/"

rep = 10

for(run in 1:rep){
  for(ts in 1:2){
    if(ts == 1){
      
      # load functions
      source(paste0(scriptsDir, "0_Functions.R"))
      
      cat("All segregating sites are causal \n")
      traitScen = "maxQTL"
      
      source(paste0(scriptsDir, "1_Founders_mt.R"))
      
    }else{
      load("founders.RData")
      cat("One segregating site is causal \n")
      
      traitScen = "minQTL"
      mQtl = 1
      mSnp = i - mQtl
      
      mtDNA <- mtDNA_sv
      rm(mtDNA_sv, SP2)
      
      SP2 = SimParam$new(mtDNA)$
        restrSegSites(overlap=TRUE)$
        addTraitA(nQtlPerChr = mQtl, mean = 0, var=varM)$
        addSnpChip(nSnpPerChr = mSnp)
      mtDNA <- newPop(mtDNA, simParam=SP2)
      
      # store haplotyopes, define maternal lineages and haploGroups
      mtfile <- as.matrix(pullSegSiteHaplo(mtDNA, simParam = SP2))
      mtfile <- as.data.frame(mtfile)
      mtfile <- mtfile %>% unite(haplotype, 1:length(mtfile), sep = "")
      
      #sum(!duplicated(mtfile))
      y <- unique(mtfile)
      
      y <- rownames_to_column(y, "haploGroup")
      
      mtfile <- mtfile %>% mutate(HaploGroup = with(y, haploGroup[match(mtfile$haplotype,  y$haplotype)]), 
                                  ML = mtDNA@id, mTbv = mtDNA@gv)
      
      mtfile <- mtfile %>% filter(!duplicated(HaploGroup)) %>% select(-HaploGroup)
      
      mtDNA <- selectInd(mtDNA, nInd = nrow(mtfile), simParam = SP2, candidates = mtfile$ML, use="rand")
      rm(y)
      
    }
    
    # Generate base population
    source(paste0(scriptsDir, "2_burnin.R"))

    # Breeding Scenarios
    # (1) Standard Progeny-Testing:
    rm(list = ls())
    load("burnin.RData")
    file.remove("Blupf901.dat", "renumf90.par")
    preparePAR(model="PED")
    nBreedingGen = 20
    program = "PED"
    model = "std"
    
    source(paste0(scriptsDir, "3_BreedScheme.R"))
    # save data
    rm(list=setdiff(ls(), c("Accuracy", "Covars", "model", "program", "traitScen", "run", "ts")))
    assign(paste0("acc_",model, program, traitScen, run), Accuracy); rm(Accuracy)
    assign(paste0("cov_",model, program, traitScen, run), Covars); rm(Covars)
    
    save.image(file=paste0(model, program, traitScen, run, ".RData"))
    
    # (2) Mitochondria Progeny-Testing: 
    rm(list = ls())
    load("burnin.RData")
    file.remove("Blupf901.dat", "renumf90.par")
    preparePAR(model="PEDmt")
    nBreedingGen = 20
    program = "PED"
    model = "mt"
    
    source(paste0(scriptsDir, "3_BreedScheme.R"))
    # save data
    rm(list=setdiff(ls(), c("Accuracy", "Covars", "model", "program", "traitScen", "run", "ts")))
    assign(paste0("acc_",model, program, traitScen, run), Accuracy); rm(Accuracy)
    assign(paste0("cov_",model, program, traitScen, run), Covars); rm(Covars)
    
    save.image(file=paste0(model, program, traitScen, run, ".RData"))
    
    # (3) Standard Genomic-Testing:
    rm(list = ls())
    load("burnin.RData")
    file.remove("Blupf901.dat", "renumf90.par")
    preparePAR(model="GEN")
    nBreedingGen = 20
    program = "GEN"
    model = "std"
    
    active<- c(do.call(c, unlist(breeding$eliteSire, recursive = FALSE)),
               do.call(c, unlist(breeding$progenyTest, recursive = FALSE)),
               do.call(c, unlist(breeding$youngBulls, recursive = FALSE)),
               selectInd(c(do.call(c, unlist(breeding$eliteDam, recursive=FALSE)),
                           do.call(c, unlist(breeding$youngFemales, recursive=FALSE)), 
                           do.call(c, unlist(breeding$Cows, recursive = FALSE))), 8944, use="rand"))
    # create training population considering all males and random females
    snpData(active)
    
    source(paste0(scriptsDir, "3_BreedScheme.R"))
    # save data
    rm(list=setdiff(ls(), c("Accuracy", "Covars", "model", "program", "traitScen", "run", "ts")))
    assign(paste0("acc_",model, program, traitScen, run), Accuracy); rm(Accuracy)
    assign(paste0("cov_",model, program, traitScen, run), Covars); rm(Covars)
    
    save.image(file=paste0(model, program, traitScen, run, ".RData"))
    
    # (4) Mitochondria Genomic-Testing:
    rm(list = ls())
    load("burnin.RData")
    file.remove("Blupf901.dat", "renumf90.par")
    preparePAR(model="GENmt")
    mtdnaGinv(mtDNA, SP2)
    nBreedingGen = 20
    program = "GEN"
    model = "mt"
    
    active<- c(do.call(c, unlist(breeding$eliteSire, recursive = FALSE)),
               do.call(c, unlist(breeding$progenyTest, recursive = FALSE)),
               do.call(c, unlist(breeding$youngBulls, recursive = FALSE)),
               selectInd(c(do.call(c, unlist(breeding$eliteDam, recursive=FALSE)),
                           do.call(c, unlist(breeding$youngFemales, recursive=FALSE)), 
                           do.call(c, unlist(breeding$Cows, recursive = FALSE))), 8944, use="rand"))
    # create training population considering all males and random females
    snpData(active)
    
    source(paste0(scriptsDir, "3_BreedScheme.R"))
    # save data
    rm(list=setdiff(ls(), c("Accuracy", "Covars", "model", "program", "traitScen", "run", "ts")))
    assign(paste0("acc_",model, program, traitScen, run), Accuracy); rm(Accuracy)
    assign(paste0("cov_",model, program, traitScen, run), Covars); rm(Covars)
    
    save.image(file=paste0(model, program, traitScen, run, ".RData"))
    
    #clear old files
    file.remove("Blupf901.dat","Blupf90.ped","blup.log","burnin.RData","freqdata.count","Gen_call_rate","mrk2.tmp", 
                "mrk.tmp","mtdnaGinv.txt","newfile","renadd02.ped","renf90.dat", 
                "renf90.fields","renf90.inb","renf90.par","renf90.tables","renumf90.par", 
                "renum.log","solutions","sum2pq","snp.dat","snp.dat_XrefID")
  }
}
