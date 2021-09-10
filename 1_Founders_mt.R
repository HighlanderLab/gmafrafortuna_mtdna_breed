# Evaluation of the impact of accounting for mitochondrial variation on breeding values estimation for dairy cattle
# Gabriela Mafra Fortuna
# Highlander Lab
# The Roslin Institute 
# July 2020 - updated Aug 2021
# ------------------------------------------------------------------------------------------------------------------------
# ---------------------------------- Generate initial haplotypes & nuclear founder pop ----------------------------------- 
# ---------------------------------- Nuclear ----------------------------------- 
cat("Simulating nuclear haplotypes...\n")
FOUNDERPOP <- runMacs(nInd = nFem*2,
                      nChr = nChr,
                      segSites = nQtl+nSnp,
                      species = "CATTLE",
                      ploidy = 2L)

SP = SimParam$new(FOUNDERPOP)$
  restrSegSites(nQtl, nSnp, overlap=FALSE)$
  addTraitA(nQtlPerChr=nQtl, mean=0, var=varA)$
  setSexes("yes_sys")$
  addSnpChip(nSnp)

cat("Simulating founder population...\n")
founders <- newPop(FOUNDERPOP)
rm(FOUNDERPOP)

# ---------------------------------- Mitochondrial -----------------------------------
cat("Simulating mitochondrial haplotypes...\n")
i = 0
while(i < 400){
  mtDNA = runMacs2(nInd = mtInd, 
                   nChr = 1,
                   Ne = mNe,
                   bp = 16202,
                   genLen = 0,
                   mutRate = mu,
                   histNe = histMtNe,
                   histGen = histMtGen,
                   ploidy = 1)
  i = mtDNA@nLoci
}

# Qtl and Snp par depending on scenario
mQtlmax = i
mSnpmax = mQtlmax

mQtlmin = 1
mSnpmin = i

SP2 = SimParam$new(mtDNA)$
  restrSegSites(overlap=TRUE)$
  addTraitA(nQtlPerChr=mQtlmax, mean = 0, var = varM)$
  addSnpChip(nSnpPerChr = mSnpmax)

SP3 = SimParam$new(mtDNA)$
  restrSegSites(overlap=TRUE)$
  addTraitA(nQtlPerChr=mQtlmin, mean = 0, var = varM)$
  addSnpChip(nSnpPerChr = mSnpmin)

cat("Saving data as 'founders.RData'...\n")
save.image("founders.RData")



















