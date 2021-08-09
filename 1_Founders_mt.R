# mtDNA Project: Should mitochondrial variation be accounted for on BVE
# Dairy Cattle Breeding Simulations
# Roslin Institute
# from July 2020
# Create Founder Haplotypes (nuclear and mitochondrial)

#library(AlphaSimR)
#library(tidyverse)
#library(Matrix)

# Simulating nuclear genome
cat("Simulating nuclear genome...\n")
FOUNDERPOP <- runMacs(nInd = nInd*2,
                      nChr = nChr,
                      segSites = nQtl+nSnp,
                      species = "CATTLE",
                      ploidy = 2L)

SP = SimParam$new(FOUNDERPOP)$
   restrSegSites(nQtl, nSnp, overlap=FALSE)$
   addTraitA(nQtlPerChr=nQtl, mean = 0, var = varA)$
   setSexes("yes_sys")$
   addSnpChip(nSnp)
founders <- newPop(FOUNDERPOP)


# Simulating mtDNA 
cat("Simulating mitogenome...\n")
i = 0
while (i < 400 ) {
   mtDNA = runMacs2(nInd=mtInd, nChr=1, 
                    Ne= mNe,
                    bp=16202,
                    genLen = 0,
                    mutRate = mu, 
                    histNe = histMtNe,
                    histGen = histMtGen,
                    ploidy=1)
   i = mtDNA@nLoci
}

mtDNA_sv = mtDNA
# solution to run different scenarios

mQtl = i
mSnp = mQtl
      
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

save.image("founders.RData")
