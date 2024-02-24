#
# mtDNA results 100 rep
# Accuracy by animal category
#

library(tidyverse)

# holder dataset
df = NULL

# list directories in folder (replicates)
folders = list.dirs()
# remove ./
folders = folders[-1]

# iterate over replicates folders
for(folder in folders){
  
  # list all files in folder for the maxQTL scenario
  files = list.files(path=folder,pattern = "*maxQTL_\\d*.RData",full.names = "TRUE")
  
  # iterate over files
  for(file in files){
    
    # load file and summarise category results
    load(file)
    
    # filter the last year of breeding cycle,
    # select col Category, Scenario, IIds and nEbv
    # add replicate col (REP)
    # bring nTbv from Records table 
    tmp = catSummary %>% 
      filter(Breeding_year == 40) %>%
      select(Category, Scenario, IIds, nEbv) %>% 
      mutate(REP = rep) %>% 
      unnest_longer(c(IIds, nEbv))%>% 
      mutate(nTbv = Records$gv_corr[match(IIds, Records$IId)]) %>%
      select(-IIds)
    
    # combine to holer dataset
    df = rbind(df, tmp)
    
  }
}

# save results
write_csv(df, "category_sum.csv")

# select only STANDARD and BASELINE scenarios
# calculate accuray by category/scenario/rep
df2 = df %>% 
  filter(Scenario %in% c("GENmtbaseline", 
                         "GENstdbaseline", 
                         "PEDmtbaseline", 
                         "PEDstdbaseline")) %>% 
  group_by(Category, Scenario, REP) %>% 
  summarise(acc = cor(nEbv, nTbv))
#write_csv(df2, "acc_cat_100.csv")

#df = read_csv("~/Desktop/eddie_mount/data/acc_cat_100.csv")
#head(df)

# Split the Scenario variable into PROGRAM and SCENARIO
# rename proven_bulls and young_bulls from GS Program
df2 = df %>% mutate(
  Program = ifelse(Scenario %in% c("GENmtbaseline", "GENstdbaseline"), 
                   "GS", "PT"),
  Scenario = ifelse(Scenario %in% c("GENstdbaseline", "PEDstdbaseline"), 
                    "STANDARD", "BASELINE"),
  Category = ifelse(Category %in% c("proven_bulls") & Program %in% c("GS"),
                    "proven_bulls_geno",
                    ifelse(Category %in% c("young_bulls") & Program %in% c("GS"),
                           "young_bulls_geno", Category)))

# calculate mean/sd for replicates and diff in means
df3 = df2 %>% 
  group_by(Category, Program, REP) %>% 
  mutate(diff = lag(acc, default = first(acc)) - acc) %>%
  ungroup() 

df3 = df3 %>%
  group_by(Category, Scenario, Program) %>%
  summarise(mean_acc  = mean(acc),
            sd_acc    = sd(acc),
            mean_diff = mean(diff),
            sd_diff   = sd(diff),
            .groups = "keep")

# write final table to file
write_csv(df3,
          "~/Desktop/mDNA_jds_submission/review/final_acc_100r.csv")
