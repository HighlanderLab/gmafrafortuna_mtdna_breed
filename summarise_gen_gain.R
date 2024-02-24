# module load roslin/R/4.2.1
library(tidyverse)
df = NULL
# list directories in folder (replicates)
folders = list.dirs()
# remove ./
folders = folders[-1]

for(folder in folders){
  
  # list files in folder (scenarios) that match pattern (only maxQTL scenarios)
  files = list.files(path=folder, 
                     pattern = "*maxQTL_\\d*.RData", full.names = "TRUE")
  
  # iterate over files
  for(file in files){
    
    load(file)
    
    tmp = Records %>% filter(Generation == 30 | Generation == 50) %>%
      select(Generation, Program, Model, Selection, mTbv, nTbv) %>%
      unite("Scenario", Program:Selection, remove = TRUE) %>%
      mutate(tTbv = nTbv + mTbv)
    
    sigma_t0 = sd(tmp %>% filter(Generation == min(Generation)) %>% pull(tTbv))
    
    means = tmp %>% group_by(Generation, Scenario) %>% 
      summarise_at(vars(c(mTbv, nTbv, tTbv)), mean)
    
    df = rbind(df,
               tibble(Scenario = means$Scenario[2],
                      REP = rep,
                      nTbv = (means$nTbv[means$Generation == max(means$Generation)] - 
                                means$nTbv[means$Generation == min(means$Generation)])/sigma_t0,
                      mTbv = (means$mTbv[means$Generation == max(means$Generation)] - 
                                means$mTbv[means$Generation == min(means$Generation)])/sigma_t0,
                      tTbv = (means$tTbv[means$Generation == max(means$Generation)] - 
                                means$tTbv[means$Generation == min(means$Generation)])/sigma_t0,
               ))
  }
}

# write table to csv
write_csv(df,file = "gen_gain_std_final.csv", quote = "none")

df = read_csv("~/Desktop/eddie_mount/data/gen_gain_std_final.csv")

mn = df %>% group_by(Scenario) %>% 
  summarise_at(vars(nTbv, mTbv, tTbv), mean) %>%
  pivot_longer(
    cols = c(nTbv, tTbv, mTbv), 
    names_to = "Genome", 
    values_to = "Gain"
  )
head(mn)

df = df %>% group_by(Scenario) %>% 
  summarise_at(vars(nTbv, mTbv, tTbv), sd) %>%
  pivot_longer(
    cols = c(nTbv, tTbv, mTbv), 
    names_to = "Genome", 
    values_to = "sdd"
  ) %>% left_join(mn)

df = df %>% mutate(Program = ifelse(Scenario == "GEN_mt_baseline" |
                              Scenario == "GEN_mt_extreme" |
                              Scenario == "GEN_mt_optimum" |
                              Scenario == "GEN_std_baseline", 
                            "Genomic Selection",
                            "Progeny Test"),
           Scenario = ifelse(Scenario == "GEN_mt_baseline" | 
                               Scenario == "PED_mt_baseline"  , "BASELINE",
                             ifelse(Scenario == "GEN_mt_extreme" |
                                      Scenario == "PED_mt_extreme", 
                                    "EXTREME",
                                    ifelse(Scenario == "GEN_mt_optimum" |
                                             Scenario == "PED_mt_optimum", 
                                           "OPTIMUM", 
                                           "STANDARD"))))

target = c("STANDARD", "BASELINE", "OPTIMUM", "EXTREME")
target2 = c("Progeny Test", "Genomic Selection")
df$Scenario <- factor(as.factor(df$Scenario),
                   levels = target)
df$Program = factor(as.factor(df$Program), levels = target2)

pd <- position_dodge(0.001)
q = ggplot(df, aes(x = Scenario, y = Gain,colour = Program, shape = Genome, fill = Program))+
  geom_errorbar(aes(ymin = Gain - sdd, ymax = Gain + sdd),
                width=0.6,
                size=0.5,
                position = position_dodge(width = .9))+
  geom_point( size = 3, position = position_dodge(width = .9))+
  labs(y="Genetic gain (sd)", x="Scenario") 
q

library(extrafont)
library(ggbreak)
font_import()
loadfonts(device="pdf")  
cl <- c(
  NA,
  "black"
)
sp <- c(1, 24)

p = q + scale_shape_manual(values = c(21, 24, 22)) +
  scale_color_manual(values = c("black", "black"))+
  scale_fill_manual(values = c("white", "black"))+
  theme_bw() +
  theme(text=element_text(family="Times New Roman", face="bold", size=22))+
  theme(axis.text.x = element_text(size=22, angle = 0, vjust = 1, hjust=0.5))+
  theme(legend.text = element_text(size=26),
        legend.position="top",
        axis.text = element_text(size=22),
        axis.title = element_text(size=26),
        plot.margin = unit(c(1,1,1,1), "cm"))+
  guides(fill = guide_legend(override.aes = list(shape=21)))+ 
  scale_y_break(c(0.45, 2))

p

PaperSize=15
ggsave(plot = p, filename = "~/Desktop/mDNA_jds_submission/review/Genetic_gain_std_paper_100r.png",
       height = PaperSize, width = PaperSize * 2.5, unit = "cm")

