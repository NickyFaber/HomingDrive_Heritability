rm(list = ls())

########################################
########## Environment set-up ##########
########################################

setwd("/Users/s2018147/Library/CloudStorage/OneDrive-WageningenUniversity&Research/PhD/03_Genetic_diversity")

# Set this to either get nice plots or correct statistics (!). Rerun data processing before switching between tasks.
doin_stats = TRUE

library(tidyverse)
library(glmmTMB)
library(bbmle)
library(DHARMa)
library(janitor)
library(multcompView)
library(car)
library(lme4)
library(merTools)
library(sjPlot)
library(insight)
library(glmmTMB)
library(patchwork)
library(parameters)
library(emmeans)

# Figure aesthetics
PaperTheme <- theme_bw(base_size = 11, base_family = "sans") + 
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        plot.title=element_text(size=12, hjust=0.5), 
        legend.title=element_text(size=12),
        legend.position = "bottom", 
        legend.justification = "center",
        axis.title=element_text(size=12),
        plot.background = element_rect(fill = "transparent", color=NA), # region outside plot
        panel.background = element_rect(fill = "white"), # region inside plot
        legend.background = element_rect(colour = 'transparent', fill='transparent'), # transparent legends
        legend.box.background = element_rect(fill = "white", color="white")) # white legend box

#########################################################
########## Loading, merging, and cleaning data ##########
#########################################################

##### Read data files
results <- read.csv2(file = "scored_Nicky.csv", na.strings = "")
resultz <- results[-1:-2,]
colnames(resultz) <- results[2,]

extra_results <- read.csv2(file = "scored_Fleur_Jimmie.csv", na.strings = "")
extra_resultz <- extra_results[-1:-2,]
colnames(extra_resultz) <- extra_results[2,]

rm(results, extra_results)

##### Checks for duplicates and removing them all
duplicates1 <- get_dupes(resultz, population, line, mother, cross, genotype) # Ten duplicate labels
duplicates2 <- get_dupes(extra_resultz, population, line, mother, cross, genotype) # Two duplicate labels

resultz <- resultz %>% 
  group_by(population, line, mother, cross, genotype) %>% 
  filter(n() == 1)
extra_resultz <- extra_resultz %>% 
  group_by(population, line, mother, cross, genotype) %>% 
  filter(n() == 1)

##### Merging data files
orphan_vials <- anti_join(extra_resultz, resultz, by = c("population", "line", "mother", "cross", "genotype"))

merged <- full_join(resultz, extra_resultz, by = c("population", "line", "mother", "cross", "genotype")) %>% 
  mutate(fliesScoredBy = coalesce(fliesScoredBy.x, fliesScoredBy.y),
         
         males_whiteeyes_RG = coalesce(males_whiteeyes_RG.x, males_whiteeyes_RG.y),
         males_whiteeyes_G = coalesce(males_whiteeyes_G.x, males_whiteeyes_G.y),
         males_whiteeyes_R = coalesce(males_whiteeyes_R.x, males_whiteeyes_R.y),
         males_whiteeyes_N = coalesce(males_whiteeyes_N.x, males_whiteeyes_N.y),
         
         females_redeyes_R = coalesce(females_redeyes_R.x, females_redeyes_R.y),
         females_redeyes_NR = coalesce(females_redeyes_NR.x, females_redeyes_NR.y),
         females_redeyes_unclear = coalesce(females_redeyes_unclear.x, females_redeyes_unclear.y),
         
         males_redeyes_unclear = coalesce(males_redeyes_unclear.x, males_redeyes_unclear.y),
         females_whiteeyes_NR = coalesce(females_whiteeyes_NR.x, females_whiteeyes_NR.y),
         offspring_unclear = coalesce(offspring_unclear.x, offspring_unclear.y),
         
         unhatched = sum(as.numeric(unhatched.x), as.numeric(unhatched.y), na.rm = TRUE),
         dead = sum(as.numeric(dead.x), as.numeric(dead.y), na.rm = TRUE)) %>% 
  dplyr::select(-males_whiteeyes_RG.x,-males_whiteeyes_G.x,-males_whiteeyes_R.x,-males_whiteeyes_N.x,-females_redeyes_R.x,
         -females_redeyes_NR.x,-females_redeyes_unclear.x,-males_redeyes_unclear.x,-females_whiteeyes_NR.x,
         -offspring_unclear.x, -unhatched.x, -dead.x, -fliesScoredBy.x,
         -males_whiteeyes_RG.y,-males_whiteeyes_G.y,-males_whiteeyes_R.y,-males_whiteeyes_N.y,-females_redeyes_R.y,
         -females_redeyes_NR.y,-females_redeyes_unclear.y,-males_redeyes_unclear.y,-females_whiteeyes_NR.y,
         -offspring_unclear.y, -unhatched.y, -dead.y, -fliesScoredBy.y) %>%
  dplyr::select(-vialNr.y) %>%
  rename(vialNr = vialNr.x,
         notes.vials = notes.x,
         notes.flies = notes.y) %>%
  relocate(c(vialNr,vialsScoredBy,fliesScoredBy), .before = population) %>%
  relocate(c(mediumLooks, notes.vials, notes.flies), .after = dead)

no_vial_info <- merged %>%
  filter(is.na(vialsScoredBy))

no_fly_info <- merged %>%
  filter(is.na(fliesScoredBy))

weird_single_replicate <- merged %>%
  filter(population == "AK" & line == "6")

unreadable_label <- merged %>%
  filter(mother == "?")

failed_crosses <- merged %>%
  # Filter out any crosses with red-eyed males or white-eyes females (indicating that virgin flies were in fact not virgin)
  filter(!(males_redeyes_unclear == "0" & females_whiteeyes_NR == "0"))

results_complete <- merged %>%
  filter(!is.na(vialsScoredBy) & !is.na(fliesScoredBy)) %>%
  filter(!(population == "AK" & line == "6")) %>%
  filter(males_redeyes_unclear == "0" & females_whiteeyes_NR == "0") %>%
  filter(!(mother == "?"))

##### Recode line numbers and mother numbers so it looks nice in plots
results_complete$line[results_complete$population=="KA" & results_complete$line==7] <- 5
results_complete$line[results_complete$population=="MA" & results_complete$line==6] <- 2
results_complete$line[results_complete$population=="MA" & results_complete$line==7] <- 3
results_complete$line[results_complete$population=="MU" & results_complete$line==6] <- 2
results_complete$line[results_complete$population=="RE" & results_complete$line==6] <- 5
results_complete$line[results_complete$population=="UM" & results_complete$line==6] <- 2

results_complete$mother[results_complete$population=="AK" & results_complete$line==1 & results_complete$mother==6] <- 1
results_complete$mother[results_complete$population=="AK" & results_complete$line==2 & results_complete$mother==6] <- 3
results_complete$mother[results_complete$population=="AK" & results_complete$line==4 & results_complete$mother==6] <- 2
results_complete$mother[results_complete$population=="AK" & results_complete$line==5 & results_complete$mother==6] <- 5
results_complete$mother[results_complete$population=="KA" & results_complete$line==2 & results_complete$mother==6] <- 2
results_complete$mother[results_complete$population=="KA" & results_complete$line==4 & results_complete$mother==6] <- 2
results_complete$mother[results_complete$population=="MA" & results_complete$line==2 & results_complete$mother==6] <- 2
results_complete$mother[results_complete$population=="MA" & results_complete$line==2 & results_complete$mother==5] <- 4
results_complete$mother[results_complete$population=="MA" & results_complete$line==3 & results_complete$mother==6] <- 3
results_complete$mother[results_complete$population=="MA" & results_complete$line==4 & results_complete$mother==6] <- 2
results_complete$mother[results_complete$population=="MA" & results_complete$line==5 & results_complete$mother==6] <- 2
results_complete$mother[results_complete$population=="MU" & results_complete$line==1 & results_complete$mother==6] <- 1
results_complete$mother[results_complete$population=="MU" & results_complete$line==2 & results_complete$mother==6] <- 3
results_complete$mother[results_complete$population=="MU" & results_complete$line==3 & results_complete$mother==6] <- 2
results_complete$mother[results_complete$population=="MU" & results_complete$line==4 & results_complete$mother==6] <- 3
results_complete$mother[results_complete$population=="MU" & results_complete$line==5 & results_complete$mother==6] <- 1
results_complete$mother[results_complete$population=="RE" & results_complete$line==1 & results_complete$mother==6] <- 2
results_complete$mother[results_complete$population=="RE" & results_complete$line==2 & results_complete$mother==4] <- 2
results_complete$mother[results_complete$population=="RE" & results_complete$line==5 & results_complete$mother==6] <- 1
results_complete$mother[results_complete$population=="UM" & results_complete$line==1 & results_complete$mother==6] <- 3
results_complete$mother[results_complete$population=="UM" & results_complete$line==4 & results_complete$mother==6] <- 1

##### Putting everything in the right format for either plotting or statistics
if (doin_stats) {
  results_complete$population <- factor(results_complete$population, 
                                        levels = c("AK","KA","MA","MU","RE","UM"),
                                        labels = c("AK","KA","MA","MU","RE","UM"))
  
  results_complete$line <- paste0(results_complete$population, "_", 
                                  results_complete$line)
  results_complete$mother <- paste0(results_complete$line, "_", 
                                    results_complete$mother, "_",
                                    results_complete$genotype)
  results_complete$cross <- paste0(results_complete$mother, "_", 
                                   results_complete$cross)
} else {
  results_complete$population <- factor(results_complete$population, 
                              levels = c("AK","KA","MA","MU","RE","UM"),
                              labels = c("Akaa,\nFinland","Karensminde,\nDenmark","Mauternbach,\nAustria",
                                           "Munich,\nGermany","Recarei,\nPortugal","Uman,\nUkraine"))
}

results_complete$vialNr <- factor(results_complete$vialNr)
results_complete$line <- factor(results_complete$line)
results_complete$mother <- factor(results_complete$mother)
results_complete$cross <- factor(results_complete$cross)
results_complete$genotype <- factor(results_complete$genotype, 
                                    levels = c("CO", "GD"),
                                    labels = c("Control", "Gene drive"))
results_complete$vialsScoredBy <- factor(results_complete$vialsScoredBy, 
                                         levels = c("Nicky", "Gabriella","Fleur", "Jimmie"),
                                         labels = c("Nicky", "Gabriella","Fleur", "Jimmie"))
results_complete$fliesScoredBy <- factor(results_complete$fliesScoredBy, 
                                         levels = c("Nicky", "Gabriella","Fleur", "Jimmie"),
                                         labels = c("Nicky", "Gabriella","Fleur", "Jimmie"))

results_complete$males_whiteeyes_RG <- as.numeric(results_complete$males_whiteeyes_RG)
results_complete$males_whiteeyes_R <- as.numeric(results_complete$males_whiteeyes_R)
results_complete$males_whiteeyes_G <- as.numeric(results_complete$males_whiteeyes_G)
results_complete$males_whiteeyes_N <- as.numeric(results_complete$males_whiteeyes_N)
results_complete$females_redeyes_R <- as.numeric(results_complete$females_redeyes_R)
results_complete$females_redeyes_NR <- as.numeric(results_complete$females_redeyes_NR)
results_complete$females_redeyes_unclear <- as.numeric(results_complete$females_redeyes_unclear)
results_complete$males_redeyes_unclear <- as.numeric(results_complete$males_redeyes_unclear)
results_complete$females_whiteeyes_NR <- as.numeric(results_complete$females_whiteeyes_NR)
results_complete$offspring_unclear <- as.numeric(results_complete$offspring_unclear)
results_complete$unhatched <- as.numeric(results_complete$unhatched)
results_complete$dead <- as.numeric(results_complete$dead)

###########################################
########## Calculate inheritance ##########
###########################################

results_calculated <- results_complete %>%
# Calculate offspring numbers and inheritance rates
  mutate(totalOffspring = females_redeyes_R + females_redeyes_NR + females_redeyes_unclear + 
           males_whiteeyes_RG + males_whiteeyes_R + males_whiteeyes_G + males_whiteeyes_N +
           unhatched + dead + offspring_unclear,
         totalAliveOffspring = females_redeyes_R + females_redeyes_NR + females_redeyes_unclear + 
           males_whiteeyes_RG + males_whiteeyes_R + males_whiteeyes_G + males_whiteeyes_N + offspring_unclear, 
         totalMaleOffspring = males_whiteeyes_RG + males_whiteeyes_R + males_whiteeyes_G + males_whiteeyes_N,
         
         cas9Inheritance = (males_whiteeyes_RG + males_whiteeyes_G) / totalMaleOffspring,
         GDInheritance = (males_whiteeyes_RG + males_whiteeyes_R) / totalMaleOffspring,
         sexRatio = totalMaleOffspring / totalSexedOffspring)

# Remove incorrect crosses (no GD or no Cas9)
results_calculated <- results_calculated %>%
  filter(!(genotype == "Gene drive" & totalMaleOffspring > 0 & (GDInheritance == 0 | cas9Inheritance == 0)))

# Number of correct crosses
nrow(results_calculated)
nrow(results_calculated[results_calculated$genotype=="Gene drive",])
nrow(results_calculated[results_calculated$genotype=="Control",])

# Gene drive crosses with at least one male offspring
results_calculated_drive <- results_calculated %>%
  filter(genotype == "Gene drive" & totalMaleOffspring > 0)
nrow(results_calculated_drive)

################################
########## Plot plots ##########
################################

#------------------------------------------####
#--------- Gene drive inheritance ---------####
#------------------------------------------####

pdf(file = "GD.pdf", useDingbats = TRUE)

# By population

p1 <- ggplot() +
  geom_hline(yintercept = 0.5, linetype = 3, colour = "grey") +
  geom_violin(data = results_calculated_drive, aes(x = population, y = GDInheritance), scale = "count", bw = "nrd0") + 
  geom_boxplot(data = results_calculated_drive, aes(x = population, y = GDInheritance), 
               position=position_dodge(width=1), linewidth = 0.4, width = 0.25, outliers = FALSE) +
  geom_point(data = results_calculated_drive, aes(x = population, y = GDInheritance, 
                                                                colour = population, size = totalMaleOffspring), 
             alpha = 0.5, position = position_jitter(width = 0.25)) +
  scale_color_manual(values = c(viridisLite::viridis(6)[1:5],"#ffd416")) +
  scale_size(range = c(0.5, 3)) +
  ylim(c(0,NA)) +
  xlab("Population") +
  ylab("Gene drive inheritance") +
  PaperTheme + guides(colour = "none", size = guide_legend("Total scored offspring")) + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
p1

ggsave(plot = p1, filename = "graphs_paper/inheritance_pop.pdf", height = 15, width = 20.7, unit = "cm")

# By population and line

p2 <- ggplot() +
  facet_grid(population ~ .) +
  geom_hline(yintercept = 0.5, linetype = 3, colour = "grey") +
  geom_violin(data = results_calculated_drive, aes(x = line, y = GDInheritance), scale = "count", bw = "nrd0") + 
  geom_boxplot(data = results_calculated_drive, aes(x = line, y = GDInheritance), 
               position=position_dodge(width=1), linewidth = 0.4, width = 0.25, outliers = FALSE) +
  geom_point(data = results_calculated_drive, aes(x = line, y = GDInheritance, colour = population, size = totalMaleOffspring), 
             alpha = 0.5, position = position_jitter(width = 0.25)) +
  scale_color_manual(values = c(viridisLite::viridis(6)[1:5],"#ffd416")) +
  scale_size(range = c(0.5, 3)) +
  ylim(c(0,NA)) +
  xlab("Isofemale line number") +
  ylab("Drive inheritance") +
  PaperTheme + guides(colour = "none", size = guide_legend("Total scored offspring"))
p2

ggsave(plot = p2, filename = "graphs_paper/inheritance_pop_line.pdf", height = 25, width = 20.7, unit = "cm")

# By genotype and population and line and mother

p3 <- ggplot() +
  facet_grid(population ~ line) +
  geom_boxplot(data = results_calculated_drive, aes(x = mother, y = GDInheritance, group = interaction(line, mother)), 
               position=position_dodge(width=1), linewidth = 0.2, width = 0.5, outliers = FALSE) +
  geom_point(data = results_calculated_drive, aes(x = mother, y = GDInheritance, group = interaction(line, mother), 
                                                                colour = population, size = totalMaleOffspring), 
             alpha = 0.5, position = position_jitter(width = 0.25)) +
  geom_hline(yintercept = 0.5, colour = "grey", linetype = 3) +
  scale_color_manual(values = c(viridisLite::viridis(6)[1:5],"#ffd416")) +
  scale_size(range = c(0.5, 4)) +
  ylim(c(0,NA)) +
  xlab("Sibling batch") +
  ylab("Gene drive inheritance") +
  ggtitle("Isofemale line") +
  PaperTheme + guides(colour = "none", size = guide_legend("Total scored offspring"))
p3

ggsave(plot = p3, filename = "graphs_paper/inheritance_alllevels.pdf", height = 20, width = 20.7, unit = "cm")

# Correlation between gene drive inheritance and total number of offspring

p4 <- ggplot() +
  geom_hline(yintercept = 0.5, linetype = 3, colour = "grey") +
  geom_smooth(data = results_calculated_drive, aes(x = totalOffspring, y = GDInheritance), colour = "black", method = "lm", se = TRUE) +
  geom_point(data = results_calculated_drive, aes(x = totalOffspring, y = GDInheritance)) +
  scale_color_viridis_d(name = "Population") +
  ylim(c(0,1)) +
  xlab("Total number of offspring") +
  ylab("Drive inheritance") +
  PaperTheme
p4

ggsave(plot = p4, filename = "graphs_paper/inheritance_offspring_correlation.pdf", height = 15, width = 20.7, unit = "cm")

dev.off()

#------------------------------------####
#--------- Cas9 inheritance ---------####
#------------------------------------####

pdf(file = "Cas9.pdf", useDingbats = TRUE)

# By population

p1 <- ggplot() +
  geom_hline(yintercept = 0.5, linetype = 3, colour = "grey") +
  geom_violin(data = results_calculated_drive, aes(x = population, y = cas9Inheritance), scale = "count", bw = "nrd0") + 
  geom_boxplot(data = results_calculated_drive, aes(x = population, y = cas9Inheritance), 
               position=position_dodge(width=1), linewidth = 0.4, width = 0.25, outliers = FALSE) +
  geom_point(data = results_calculated_drive, aes(x = population, y = cas9Inheritance, 
                                                  colour = population, size = totalMaleOffspring), 
             alpha = 0.5, position = position_jitter(width = 0.25)) +
  scale_color_manual(values = c(viridisLite::viridis(6)[1:5],"#ffd416")) +
  scale_size(range = c(0.5, 3)) +
  ylim(c(0,NA)) +
  xlab("Population") +
  ylab("Cas9 inheritance") +
  PaperTheme + guides(colour = "none", size = guide_legend("Total scored offspring")) + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
p1

# By population and line

p2 <- ggplot() +
  facet_grid(population ~ .) +
  geom_hline(yintercept = 0.5, linetype = 3, colour = "grey") +
  geom_violin(data = results_calculated_drive, aes(x = line, y = cas9Inheritance), scale = "count", bw = "nrd0") + 
  geom_boxplot(data = results_calculated_drive, aes(x = line, y = cas9Inheritance), 
               position=position_dodge(width=1), linewidth = 0.4, width = 0.25, outliers = FALSE) +
  geom_point(data = results_calculated_drive, aes(x = line, y = cas9Inheritance, colour = population, size = totalMaleOffspring), 
             alpha = 0.5, position = position_jitter(width = 0.25)) +
  scale_color_manual(values = c(viridisLite::viridis(6)[1:5],"#ffd416")) +
  scale_size(range = c(0.5, 3)) +
  ylim(c(0,NA)) +
  xlab("Isofemale line number") +
  ylab("Cas9 inheritance") +
  PaperTheme + guides(colour = "none", size = guide_legend("Total scored offspring"))
p2

# By genotype and population and line and mother

p3 <- ggplot() +
  facet_grid(population ~ line) +
  geom_boxplot(data = results_calculated_drive, aes(x = mother, y = cas9Inheritance, group = interaction(line, mother)), 
               position=position_dodge(width=1), linewidth = 0.2, width = 0.5, outliers = FALSE) +
  geom_point(data = results_calculated_drive, aes(x = mother, y = cas9Inheritance, group = interaction(line, mother), 
                                                  colour = population, size = totalMaleOffspring), 
             alpha = 0.5, position = position_jitter(width = 0.25)) +
  geom_hline(yintercept = 0.5, colour = "grey", linetype = 3) +
  scale_color_manual(values = c(viridisLite::viridis(6)[1:5],"#ffd416")) +
  scale_size(range = c(0.5, 4)) +
  ylim(c(0,NA)) +
  xlab("Sibling batch") +
  ylab("Cas9 inheritance") +
  ggtitle("Isofemale line") +
  PaperTheme + guides(colour = "none", size = guide_legend("Total scored offspring"))
p3

ggsave(plot = p3, filename = "graphs_paper/cas9_inheritance_alllevels.pdf", height = 20, width = 20.7, unit = "cm")

# Correlation between gene drive inheritance and total number of offspring

p4 <- ggplot() +
  geom_hline(yintercept = 0.5, linetype = 3, colour = "grey") +
  geom_smooth(data = results_calculated_drive, aes(x = totalOffspring, y = cas9Inheritance), colour = "black", method = "lm", se = TRUE) +
  geom_point(data = results_calculated_drive, aes(x = totalOffspring, y = cas9Inheritance)) +
  scale_color_viridis_d(name = "Population") +
  ylim(c(0,1)) +
  xlab("Total number of offspring") +
  ylab("Cas9 inheritance") +
  PaperTheme
p4

ggsave(plot = p4, filename = "graphs_paper/cas9_inheritance_offspring_correlation.pdf", height = 15, width = 20.7, unit = "cm")

dev.off()

#---------------------------------------####
#--------- Number of offspring ---------####
#---------------------------------------####

pdf(file = "Offspring.pdf", useDingbats = TRUE)

# By genotype

p1 <- ggplot() +
  geom_violin(data = results_calculated, aes(x = genotype, y = totalOffspring), scale = "count", bw = "nrd0") + 
  geom_boxplot(data = results_calculated, aes(x = genotype, y = totalOffspring), 
               position=position_dodge(width=1), linewidth = 0.4, width = 0.25, outliers = FALSE) +
  geom_jitter(data = results_calculated, aes(x = genotype, y = totalOffspring, colour = genotype), 
             alpha = 0.5, width = 0.25, stroke = NA) +
  scale_color_manual(values = c(viridisLite::viridis(1),"#ffd416")) +
  xlab("Genotype") +
  ylab("Number of offspring") +
  PaperTheme + guides(colour = "none")
p1

ggsave(plot = p1, filename = "graphs_paper/offspring.pdf", height = 10, width = 10.7, unit = "cm")

# By genotype and population

p2 <- ggplot() +
  facet_grid(. ~ population) +
  geom_violin(data = results_calculated, aes(x = genotype, y = totalOffspring), scale = "count", bw = "nrd0") + 
  geom_boxplot(data = results_calculated, aes(x = genotype, y = totalOffspring), 
               position=position_dodge(width=1), linewidth = 0.4, width = 0.25, outliers = FALSE) +
  geom_jitter(data = results_calculated, aes(x = genotype, y = totalOffspring, colour = genotype), 
              alpha = 0.5, width = 0.25, stroke = NA) +
  scale_color_manual(values = c(viridisLite::viridis(1),"#ffd416")) +
  xlab("Genotype") +
  ylab("Number of offspring") +
  PaperTheme + guides(colour = "none")
p2

ggsave(plot = p2, filename = "graphs_paper/offspring_pop.pdf", height = 8, width = 20.7, unit = "cm")

# By genotype and population and line

p3 <- ggplot() +
  facet_grid(population ~ line) +
  geom_violin(data = results_calculated, aes(x = genotype, y = totalOffspring), scale = "count", bw = "nrd0") + 
  geom_boxplot(data = results_calculated, aes(x = genotype, y = totalOffspring), 
               position=position_dodge(width=1), linewidth = 0.4, width = 0.25, outliers = FALSE) +
  geom_jitter(data = results_calculated, aes(x = genotype, y = totalOffspring, colour = genotype), 
              alpha = 0.5, width = 0.25, stroke = NA, size = 2) +
  scale_color_manual(values = c(viridisLite::viridis(1),"#ffd416")) +
  xlab("Genotype") +
  ylab("Number of offspring") +
  ggtitle("Isofemale line") +
  PaperTheme + guides(colour = "none")
p3

ggsave(plot = p3, filename = "graphs_paper/offspring_pop_line.pdf", height = 20, width = 20.7, unit = "cm")

# By genotype and population and line and mother

p4 <- ggplot() +
  facet_grid(population ~ line) +
  geom_boxplot(data = results_calculated, aes(x = mother, y = totalOffspring, group = interaction(genotype, mother)), 
               position=position_dodge(width=1), linewidth = 0.2, width = 0.5, outliers = FALSE) +
  geom_point(data = results_calculated, aes(x = mother, y = totalOffspring, group = interaction(genotype, mother), colour = genotype), 
             position = position_jitterdodge(dodge.width = 1), alpha = 0.5, stroke = NA, size = 2) +
  scale_color_manual(values = c(viridisLite::viridis(1),"#ffd416"),name = "Genotype") +
  xlab("Sibling batch") +
  ylab("Number of offspring") +
  ggtitle("Isofemale line") +
  PaperTheme + guides(colour = guide_legend(nrow = 1))
p4

ggsave(plot = p4, filename = "graphs_paper/offspring_alllevels.pdf", height = 20, width = 20.7, unit = "cm")

dev.off()

####################################
########## Fitting models ##########
####################################

#---------------------------------------####
#--------- Number of offspring ---------####
#---------------------------------------####

# Analysing variance components --> with population, with line, with mother, no cross
P_full <- glmmTMB(totalOffspring ~ 1 + genotype + (1 | population*genotype) + (1 | population:line*genotype) + (1 | population:line:mother),
                  family = poisson, 
                  ziformula = ~ 1,
                  data = results_calculated)

summary(P_full)
plot(simulateResiduals(P_full)) # Poisson distribution does not fit

testCategorical(simulateResiduals(P_full), catPred = results_calculated$genotype)
testUniformity(P_full)
testDispersion(P_full)
testOutliers(P_full, type = 'bootstrap')

##### What does this overdispersion look like
sims <- 100
simo=simulate(P_full, nsim = sims)
results_calculated$type = "Observed" 
results_calculated$number = 1
Simdat=results_calculated
for (sim in 1:sims) {
  print(sim)
  SimdatExtra=results_calculated
  SimdatExtra$totalOffspring=simo[[sim]]
  SimdatExtra=transform(SimdatExtra,  
                        type="Simulated",
                        number = sim)
  Simdat=rbind(Simdat, SimdatExtra) 
}
Simdat$type <- factor(Simdat$type, 
                      levels = c("Observed", "Simulated"))
Simdat$number <- factor(Simdat$number)

SimdatOutline=results_calculated
SimdatOutline$type = "Simulated"

p1 = ggplot() + 
  geom_density(data = Simdat, aes(x=totalOffspring, colour=type, group = interaction(number,type)), fill = NA, linewidth = 0.1, alpha = 0.01, trim=TRUE) + 
  geom_density(data = Simdat[Simdat$type=="Observed",], aes(x=totalOffspring, color = type, group = interaction(number,type)), fill=NA, alpha = 1, linewidth = 2, trim=TRUE) + 
  scale_colour_manual(values = c("#69b3a2", "#404080")) +
  scale_fill_manual(values = c("#69b3a2", "#404080")) +
  ylab("Density") + 
  xlab("Total number of offspring") + 
  PaperTheme + guides(colour = "none", fill = "none") +
  ggtitle("Poisson distribution with zero-inflation")
p1

# Analysing variance components --> with population, with line, with mother, no cross
NB_full <- glmmTMB(totalOffspring ~ 1 + genotype + (1 | population*genotype) + (1 | population:line*genotype) + (1 | population:line:mother),
                  family = nbinom2, 
                  ziformula = ~ 1,
                  data = results_calculated)

summary(NB_full)
plot(simulateResiduals(NB_full)) # Poisson distribution does not fit

testUniformity(NB_full)
testDispersion(NB_full)
testOutliers(NB_full, type = 'bootstrap')

##### What does this overdispersion look like
sims <- 100
simo=simulate(NB_full, nsim = sims)
results_calculated$type = "Observed" 
results_calculated$number = 1
Simdat=results_calculated
for (sim in 1:sims) {
  print(sim)
  SimdatExtra=results_calculated
  SimdatExtra$totalOffspring=simo[[sim]]
  SimdatExtra=transform(SimdatExtra,  
                        type="Simulated",
                        number = sim)
  Simdat=rbind(Simdat, SimdatExtra) 
}
Simdat$type <- factor(Simdat$type, 
                      levels = c("Observed", "Simulated"))
Simdat$number <- factor(Simdat$number)

SimdatOutline=results_calculated
SimdatOutline$type = "Simulated"

p2 = ggplot() + 
  geom_density(data = Simdat, aes(x=totalOffspring, colour=type, group = interaction(number,type)), fill = NA, linewidth = 0.1, alpha = 0.01, trim=TRUE) + 
  geom_density(data = Simdat[Simdat$type=="Observed",], aes(x=totalOffspring, color=type, group = interaction(number,type)), fill = NA, linewidth = 2, alpha = 1, trim=TRUE) + 
  scale_colour_manual(values = c("#69b3a2", "#404080")) +
  scale_fill_manual(values = c("#69b3a2", "#404080")) +
  ylab("Density") + 
  xlab("Total number of offspring") + 
  ggtitle("Negative binomial distribution with zero-inflation") +
  PaperTheme + guides(colour = "none", fill = "none")
p2

# Analysing variance components --> with population, with line, with mother, no cross
G_full <- glmmTMB(totalOffspring ~ 1 + genotype + (1 | population*genotype) + (1 | population:line*genotype) + (1 | population:line:mother),
                  family = gaussian, 
                  ziformula = ~ 1,
                  data = results_calculated)

summary(G_full)
plot(simulateResiduals(G_full))

testCategorical(simulateResiduals(G_full), catPred = results_calculated$genotype)
testUniformity(G_full)
testDispersion(G_full)
testOutliers(G_full, type = 'bootstrap')

##### What does this overdispersion look like
sims <- 100
simo=simulate(G_full, nsim = sims)
results_calculated$type = "Observed" 
results_calculated$number = 1
Simdat=results_calculated
for (sim in 1:sims) {
  print(sim)
  SimdatExtra=results_calculated
  SimdatExtra$totalOffspring=simo[[sim]]
  SimdatExtra=transform(SimdatExtra,  
                        type="Simulated",
                        number = sim)
  Simdat=rbind(Simdat, SimdatExtra) 
}
Simdat$type <- factor(Simdat$type, 
                      levels = c("Observed", "Simulated"))
Simdat$number <- factor(Simdat$number)

SimdatOutline=results_calculated
SimdatOutline$type = "Simulated"

p3 = ggplot() + 
  geom_density(data = Simdat, aes(x=totalOffspring, colour=type, group = interaction(number,type)), fill = NA, linewidth = 0.1, alpha = 0.01, trim=TRUE) + 
  geom_density(data = Simdat[Simdat$type=="Observed",], aes(x=totalOffspring, color=type, group = interaction(number,type)), fill = NA, linewidth = 2, alpha = 1, trim=TRUE) + 
  scale_colour_manual(values = c("#69b3a2", "#404080")) +
  scale_fill_manual(values = c("#69b3a2", "#404080"), name = "Data") +
  ylab("Density") + 
  xlab("Total number of offspring") + 
  ggtitle("Gaussian distribution with zero-inflation") +
  PaperTheme
p3

p <- p1 / p2 / p3 +
  plot_annotation(tag_levels = 'A'); p

ggsave(plot = p, filename = "graphs_paper/distribution_offspring.pdf", height = 20, width = 15, unit = "cm")

# Checking which model fits best
stats::anova(P_full, NB_full, G_full)
# Negative binomial has convergence issues, so only testing the other two
stats::anova(P_full, G_full)

##### Reducing the best model to the simplest form, starting with the full model
G_genotype_population.genotype_line.genotype_mother <- glmmTMB(totalOffspring ~ 1 + genotype + (1 | population*genotype) + (1 | population:line*genotype) + (1 | population:line:mother),
                   family = gaussian, 
                   ziformula = ~ 1,
                   data = results_calculated)

summary(G_genotype_population.genotype_line.genotype_mother)

# Removing genotype fixed effect
G_population.genotype_line.genotype_mother <- glmmTMB(totalOffspring ~ 1 + (1 | population*genotype) + (1 | population:line*genotype) + (1 | population:line:mother),
                              family = gaussian, 
                              ziformula = ~ 1,
                              data = results_calculated)

summary(G_population.genotype_line.genotype_mother)

# Removing genotype population interaction
G_populationgenotype_line.genotype_mother <- glmmTMB(totalOffspring ~ 1 + (1 | population+genotype) + (1 | population:line*genotype) + (1 | population:line:mother),
                                                      family = gaussian, 
                                                      ziformula = ~ 1,
                                                      data = results_calculated)

summary(G_populationgenotype_line.genotype_mother)

# Removing genotype line interaction
G_populationgenotype_linegenotype_mother <- glmmTMB(totalOffspring ~ 1 + (1 | population+genotype) + (1 | population:line+genotype) + (1 | population:line:mother),
                                                     family = gaussian, 
                                                     ziformula = ~ 1,
                                                     data = results_calculated)

summary(G_populationgenotype_linegenotype_mother)

# Removing genotype from population
G_population_linegenotype_mother <- glmmTMB(totalOffspring ~ 1 + (1 | population) + (1 | population:line+genotype) + (1 | population:line:mother),
                                    family = gaussian, 
                                    ziformula = ~ 1,
                                    data = results_calculated)

summary(G_population_linegenotype_mother)

# Removing genotype altogether
G_population_line_mother <- glmmTMB(totalOffspring ~ 1 + (1 | population) + (1 | population:line) + (1 | population:line:mother),
                                             family = gaussian, 
                                             ziformula = ~ 1,
                                             data = results_calculated)

summary(G_population_line_mother)

# Removing mother
G_population_line <- glmmTMB(totalOffspring ~ 1 + (1 | population) + (1 | population:line),
                              family = gaussian, 
                              ziformula = ~ 1,
                              data = results_calculated)

summary(G_population_line)

# Removing population
G_line <- glmmTMB(totalOffspring ~ 1 + (1 | line),
                              family = gaussian, 
                              ziformula = ~ 1,
                              data = results_calculated)

summary(G_line)

# Removing line
G_null <- glmmTMB(totalOffspring ~ 1,
                         family = gaussian, 
                         ziformula = ~ 1,
                         data = results_calculated)

summary(G_null)

# Checking which model fits best
stats::anova(G_genotype_population.genotype_line.genotype_mother, G_population.genotype_line.genotype_mother,
             G_populationgenotype_line.genotype_mother, G_populationgenotype_linegenotype_mother, G_population_linegenotype_mother,
             G_population_line_mother, G_population_line, G_line, G_null)
# G_population_line_mother has a convergence problem, removing that
stats::anova(G_genotype_population.genotype_line.genotype_mother, G_population.genotype_line.genotype_mother,
             G_populationgenotype_line.genotype_mother, G_populationgenotype_linegenotype_mother, G_population_linegenotype_mother,
             G_population_line, G_line, G_null)

#-----------------------------------------------------------####
#--------- Number of offspring variance components ---------####
#-----------------------------------------------------------####

summary(G_genotype_population.genotype_line.genotype_mother)

variances <- data.frame(level = c("population","genotype","population.genotype","population.line","population.line.genotype","population.line.mother"),
                        variance = c(9.012e+00,1.304e-07,1.134e-03,5.396e+01,1.711e-06,5.277e+01)) # Values from the model fit above
variances$level <- factor(variances$level,
                          levels = c("population","genotype","population.genotype","population.line","population.line.genotype","population.line.mother"),
                          labels = c("Population","Genotype","Population*Genotype","Isofemale line","Isofemale line*Genotype","Sibling batch"))

p1 <- ggplot(data = variances) +
  geom_segment(aes(x=level, y=variance, xend=level, yend=0)) +
  geom_point(aes(x=level, y=variance), size=4, color="#f177ae") +
  geom_text(aes(x=level, y=variance + 3, label = sprintf("%.03e", variance))) +
  xlab("Variance component") + 
  ylab("Amount of variance explained (identity)") +
  PaperTheme
p1

ggsave(plot = p1, filename = "graphs_paper/variance_explained_offspring.pdf", height = 10, width = 20.7, unit = "cm")

#------------------------------------------####
#--------- Gene drive inheritance ---------####
#------------------------------------------####

# Analysing variance components --> with population, with line, with mother, no cross
B_full <- glmmTMB(GDInheritance ~ 1 + fliesScoredBy + (1 | population) + (1 | population:line) + (1 | population:line:mother),
                  family = binomial, 
                  data = results_calculated_drive,
                  weights = totalMaleOffspring)
summary(B_full)
plot(simulateResiduals(B_full)) # There is overdispersion
testUniformity(B_full)
testDispersion(B_full)
testOutliers(B_full, type = 'bootstrap')

##### What does this overdispersion look like
sims <- 100
simo=simulate(B_full, nsim = sims)
results_calculated_drive$type = "Observed" 
results_calculated_drive$number = 1
Simdat=results_calculated_drive
for (sim in 1:sims) {
  print(sim)
  SimdatExtra=results_calculated_drive
  SimdatExtra$GDInheritance=simo[[sim]][,1]/SimdatExtra$totalMaleOffspring
  SimdatExtra=transform(SimdatExtra,  
                        type="Simulated",
                        number = sim)
  Simdat=rbind(Simdat, SimdatExtra) 
}
Simdat$type <- factor(Simdat$type, 
                      levels = c("Observed", "Simulated"))
Simdat$number <- factor(Simdat$number)

SimdatOutline=results_calculated_drive
SimdatOutline$type = "Simulated"

p1 = ggplot() + 
  geom_density(data = Simdat, aes(x=GDInheritance, colour=type, group = interaction(number,type)), fill = NA, linewidth = 0.1, alpha = 0.01) + 
  geom_density(data = Simdat[Simdat$type=="Observed",], aes(x=GDInheritance, color=type, group = interaction(number,type)), fill = NA, linewidth = 2, alpha = 1) + 
  scale_colour_manual(values = c("#69b3a2", "#404080")) +
  scale_fill_manual(values = c("#69b3a2", "#404080")) +
  xlim(c(0,1)) + 
  ylim(c(0,6.3)) +
  ylab("Density") + 
  xlab("Gene drive inheritance") + 
  ggtitle("Binomial distribution") +
  PaperTheme + guides(colour = "none", fill = "none")
p1

# We could solve overdispersion in two ways:

# Using a beta-binomial distribution
BB_full <- glmmTMB(GDInheritance ~ 1 + fliesScoredBy + (1 | population) + (1 | population:line) + (1 | population:line:mother),
                   family = betabinomial, 
                   data = results_calculated_drive,
                   weights = totalMaleOffspring)
summary(BB_full)
plot(simulateResiduals(BB_full))


##### What does this overdispersion look like
sims <- 100
simo=simulate(BB_full, nsim = sims)
results_calculated_drive$type = "Observed" 
results_calculated_drive$number = 1
Simdat=results_calculated_drive
for (sim in 1:sims) {
  print(sim)
  SimdatExtra=results_calculated_drive
  SimdatExtra$GDInheritance=simo[[sim]][,1]/SimdatExtra$totalMaleOffspring
  SimdatExtra=transform(SimdatExtra,  
                        type="Simulated",
                        number = sim)
  Simdat=rbind(Simdat, SimdatExtra) 
}
Simdat$type <- factor(Simdat$type, 
                      levels = c("Observed", "Simulated"))
Simdat$number <- factor(Simdat$number)

SimdatOutline=results_calculated_drive
SimdatOutline$type = "Simulated"

p2 = ggplot() + 
  geom_density(data = Simdat, aes(x=GDInheritance, colour=type, group = interaction(number,type)), fill = NA, linewidth = 0.1, alpha = 0.01) + 
  geom_density(data = Simdat[Simdat$type=="Observed",], aes(x=GDInheritance, color=type, group = interaction(number,type)), fill = NA, linewidth = 2, alpha = 1) + 
  scale_colour_manual(values = c("#69b3a2", "#404080")) +
  scale_fill_manual(values = c("#69b3a2", "#404080")) +
  xlim(c(0,1)) + 
  ylim(c(0,6.3)) +
  ylab("Density") + 
  xlab("Gene drive inheritance") + 
  ggtitle("Beta-binomial distribution") +
  PaperTheme + guides(colour = "none", fill = "none")
p2

# Observation-level random effect
OLR_full <- glmmTMB(GDInheritance ~ 1 + fliesScoredBy + (1 | population) + (1 | population:line) + (1 | population:line:mother) + (1 | population:line:mother:cross),
                  family = binomial, 
                  data = results_calculated_drive,
                  weights = totalMaleOffspring)
summary(OLR_full)
plot(simulateResiduals(OLR_full))

##### What does this overdispersion look like
sims <- 100
simo=simulate(OLR_full, nsim = sims)
results_calculated_drive$type = "Observed" 
results_calculated_drive$number = 1
Simdat=results_calculated_drive
for (sim in 1:sims) {
  print(sim)
  SimdatExtra=results_calculated_drive
  SimdatExtra$GDInheritance=simo[[sim]][,1]/SimdatExtra$totalMaleOffspring
  SimdatExtra=transform(SimdatExtra,  
                        type="Simulated",
                        number = sim)
  Simdat=rbind(Simdat, SimdatExtra) 
}
Simdat$type <- factor(Simdat$type, 
                      levels = c("Observed", "Simulated"))
Simdat$number <- factor(Simdat$number)

SimdatOutline=results_calculated_drive
SimdatOutline$type = "Simulated"

p3 = ggplot() + 
  geom_density(data = Simdat, aes(x=GDInheritance, colour=type, group = interaction(number,type)), fill = NA, linewidth = 0.1, alpha = 0.01) + 
  geom_density(data = Simdat[Simdat$type=="Observed",], aes(x=GDInheritance, color=type, group = interaction(number,type)), fill = NA, linewidth = 2, alpha = 1) + 
  scale_colour_manual(values = c("#69b3a2", "#404080")) +
  scale_fill_manual(values = c("#69b3a2", "#404080"), name = "Data") +
  xlim(c(0,1)) + 
  ylim(c(0,6.3)) +
  ylab("Density") + 
  xlab("Gene drive inheritance") + 
  ggtitle("Binomial distribution with an observation-level random effect") +
  PaperTheme
p3

p <- p1 / p2 / p3 +
  plot_annotation(tag_levels = 'A'); p

ggsave(plot = p, filename = "graphs_paper/distribution_GD.pdf", height = 20, width = 15, unit = "cm")

# Checking which model fits best, OLR_full is best
stats::anova(B_full, BB_full, OLR_full)

##### Reducing the best model to the simplest form, starting with the fixed effect

# Observation-level random effect
OLR_full <- glmmTMB(GDInheritance ~ 1 + fliesScoredBy + (1 | population) + (1 | population:line) + (1 | population:line:mother) + (1 | population:line:mother:cross),
                    family = binomial, 
                    data = results_calculated_drive,
                    weights = totalMaleOffspring)
summary(OLR_full)
vars <- get_variance(OLR_full, component = "all", tolerance = 1e-09)

# Removing scoredBy
OLR_population_line_mother_cross <- glmmTMB(GDInheritance ~ 1 + (1 | population) + (1 | population:line) + (1 | population:line:mother) + (1 | population:line:mother:cross),
                  family = binomial, 
                  data = results_calculated_drive,
                  weights = totalMaleOffspring)
summary(OLR_population_line_mother_cross)
vars <- get_variance(OLR_population_line_mother_cross, component = "all", tolerance = 1e-09)

# Removing population
OLR_line_mother_cross <- glmmTMB(GDInheritance ~ 1 + (1 | line) + (1 | line:mother) + (1 | line:mother:cross),
                  family = binomial, 
                  data = results_calculated_drive,
                  weights = totalMaleOffspring)
summary(OLR_line_mother_cross)
vars <- get_variance(OLR_line_mother_cross, component = "all", tolerance = 1e-09)

# Removing line
OLR_mother_cross <- glmmTMB(GDInheritance ~ 1 + (1 | mother) + (1 | mother:cross),
                               family = binomial, 
                               data = results_calculated_drive,
                               weights = totalMaleOffspring)
summary(OLR_mother_cross)
vars <- get_variance(OLR_mother_cross, component = "all", tolerance = 1e-09)

# Removing mother
OLR_cross <- glmmTMB(GDInheritance ~ 1 + (1 | cross),
                          family = binomial, 
                          data = results_calculated_drive,
                          weights = totalMaleOffspring)
summary(OLR_cross)
vars <- get_variance(OLR_cross, component = "all", tolerance = 1e-09)

# Adding cross significantly improves the model fit, the rest does not
stats::anova(OLR_full, 
             OLR_population_line_mother_cross, OLR_line_mother_cross, OLR_mother_cross, OLR_cross)

#-------------------------------------------------------------------####
#--------- Gene drive variance components and heritability ---------####
#-------------------------------------------------------------------####

# Using the model with all random effects to calculate heritability
summary(OLR_population_line_mother_cross)
vars <- get_variance(OLR_population_line_mother_cross, component = "all", tolerance = 1e-09)

population <- as.numeric(vars$var.intercept[1])
line <- as.numeric(vars$var.intercept[2])
mother <- as.numeric(vars$var.intercept[3])
cross <- as.numeric(vars$var.intercept[4])
total <- as.numeric(sum(vars$var.intercept) + vars$var.residual + vars$var.random)

# Broad-sense heritability
H2 <- (population + line) / (total); H2

variances <- data.frame(level = c("population", "isofemale line", "mother", "cross"),
                        variance = vars$var.intercept)
variances$level <- factor(variances$level,
                          levels = c("population", "isofemale line", "mother", "cross"),
                          labels = c("Population", "Isofemale line", "Sibling batch", "Individual"))

p1 <- ggplot(data = variances) +
  geom_segment(aes(x=level, y=variance, xend=level, yend=0)) +
  geom_point(aes(x=level, y=variance), size=4, color="#f177ae") +
  geom_text(aes(x=level, y=variance + 0.03, label = sprintf("%.03e", variance))) +
  xlab("Experimental level") + 
  ylab("Amount of variance explained (logit-link)") +
  PaperTheme
p1

ggsave(plot = p1, filename = "graphs_paper/variance_explained_inheritance.pdf", height = 10, width = 10.7, unit = "cm")
ggsave(plot = p1, filename = "variance_explained_inheritance.png", height = 10, width = 10.7, unit = "cm", dpi = 600)


