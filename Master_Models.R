### Wisnionski et al. "Corticosterone Varies with Metabolic Rate in Vertebrates"
### Author: Alyson R. Wisnionski
### Incept Date: 6/2/2025

###---Analyses for all PGLS relationships, statistics, and figures---###

library(tidyverse)
library(ggplot2) 
library(gt) #stats table 
library(webshot2) #gtsave
library(rphylopic) #animal phylos for plots

library(ape) #build phylogenetic tree
library(rotl) #pull from OTL
library(nlme) #gls with Brownian motion
library(geiger) #name.check
library(rr2) #R2 function
library(phytools) #phylo signal - Pagel's lambda
library(patchwork) #combining ggplots 


# Data set up -------------------------------------------------------------
#load in tree from "VertCleanTree.R"
tree <- read.nexus("VertTree.nex")

#laundry: load and clean data
#reptiles 🦎🐍🐢🐊
repdata <- read.csv("Reptile_MASTER.csv", sep = ",", header = TRUE) 
repclean <- repdata %>%
  rename(MSMR = MSMR.measured,
         Assay = Assay.Method) %>%
  mutate(Vert = "Reptile") %>%
  mutate(Species.Name = recode(Species.Name, 
                               "Cnemidophorus sexlineatus" = "Aspidoscelis sexlineatus",
                               "Cnemidophorus uniparens" = "Aspidoscelis uniparens",
                               "Crotalus helleri" = "Crotalus viridis",
                               "Egernia whitii" = "Liopholis whitii",
                               "Lacerta vivipara" = "Zootoca vivipara",
                               "Podarcis sicula" = "Podarcis siculus",
                               "Thamnophis elegans vagrans" = "Thamnophis errans",
                               "Ticiqua scincoides" = "Tiliqua scincoides"))

repmerge <- repclean %>%
  select(Species.Name, Common.Name, Lifespan, Body.Mass, MSMR, Basal.CC, Environment, Vert, Assay)

#birds 🦉🦅🐦🕊
birddata <- read.csv("Bird_MASTER.csv", sep = ",", header = TRUE)
birdclean <- birddata %>%
  rename(Assay = Assay.Method) %>%
  mutate(Vert = "Bird") %>%
  mutate(Species.Name = recode(Species.Name, 
                               "Diomedea melanophris" = "Thalassarche melanophris",
                               "Thallasarche melanophris" = "Thalassarche melanophris",
                               "Carduelis flammea" = "Acanthis flammea",
                               "Gallus domesticus" = "Gallus gallus",
                               "Grus canadensis" = "Antigone canadensis",
                               "Melozone aberti" = "Kieneria aberti",
                               "Parus caeruleus" = "Cyanistes caeruleus",
                               "Peucaea carpalis" = "Aimophila carpalis",
                               "Rhynchophanes mccownii" = "Calcarius mccownii",
                               "Spizella arborea" = "Spizelloides arborea",
                               "Thryesphilus rufalbus" = "Thryophilus rufalbus",)) %>%
  mutate(MSMR = MSMR / 1000) %>% #convert mW/g to W/g
  group_by(Species.Name) %>%
  summarise(across(everything(), #combine duplicate species by averaging values
                   ~ if (is.numeric(.)) {
                      if (n_distinct(.) > 1) mean(., na.rm = TRUE) else first(.x) } 
                   else { first(.) } )) %>%
  ungroup()

birdmerge <- birdclean %>%
  select(Species.Name, Common.Name, Lifespan, Body.Mass, MSMR, Basal.CC, Environment, Vert, Assay)

#mammals 🐭🦁🐒🦧
mammaldata <- read.csv("Mammal_MASTER.csv", sep = ",", header = TRUE)
mammalclean <- mammaldata %>%
  rename(Assay = Assay.Method) %>%
  mutate(Vert = "Mammal") %>%
  mutate(MSMR = MSMR / 1000) %>%
  mutate(Species.Name = recode(Species.Name, 
                               "Macropus rufogriseus" = "Notamacropus rufogriseus",
                               "Perognathus baileyi" = "Chaetodipus baileyi"))

mammalmerge <- mammalclean %>%
  select(Species.Name, Common.Name, Lifespan, Body.Mass, MSMR, Basal.CC, Environment, Vert, Assay)

#merge all vertebrate classes together
vertdata <- rbind(repmerge, birdmerge, mammalmerge)
rownames(vertdata) = vertdata$Species.Name #allows for the tree to map to dataframe
name.check(tree, vertdata) #OK


# Animal phylos for figures -----------------------------------------------
reptileimg <- get_phylopic(uuid = get_uuid(name = "Egernia saxatilis", n = 1)) 
birdimg <- get_phylopic(uuid = get_uuid(name = "Charadrius vociferus", n = 1))
mammalimg <- get_phylopic(uuid = get_uuid(name = "Sminthopsis murina", n = 1))


# Baseline CCST vs. Body Mass ---------------------------------------------
#---REPTILES
#drop NAs
rep_BasalCCmass_data <- repclean %>%
  drop_na(c(Basal.CC, Body.Mass))
rownames(rep_BasalCCmass_data) = rep_BasalCCmass_data$Species.Name

#drop species in the tree that are not included in this relationship
if (sum(is.na(vertdata$Basal.CC)) > 0 | sum(is.na(vertdata$Body.Mass)) > 0) {
  rep_BasalCCmass_tree <- drop.tip(tree, name.check(tree, rep_BasalCCmass_data)$tree_not_data)
} else {
  rep_BasalCCmass_tree <- tree
}

#build pgls model 
rep_BasalCCmass_pgls <- gls(log(Basal.CC) ~ log(Body.Mass), 
                            data = rep_BasalCCmass_data, 
                            correlation = corPagel(value = 0.1, 
                                                   phy = rep_BasalCCmass_tree, 
                                                   form = ~Species.Name)) 

#R2 values: specify reduced model as null model (no relationship, no phylogeny)
rep_BasalCCmass_R2_reduced <- lm(log(Basal.CC) ~ 1, 
                                 data = rep_BasalCCmass_data) 
rep_BasalCCmass_R2 <- R2(rep_BasalCCmass_pgls, rep_BasalCCmass_R2_reduced)

#plot
rep_BasalCCmass_plot<- ggplot(data = rep_BasalCCmass_data, aes(x = log(Body.Mass), y = log(Basal.CC))) + 
  geom_point(size = 3, color = "green4") +
  geom_abline(intercept = coefficients(summary(rep_BasalCCmass_pgls))[1,1], 
              slope = coefficients(summary(rep_BasalCCmass_pgls))[2,1]) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14)) +
  xlim(1, 15.5) +
  ylim(-1, 6.5) +
  labs(x = "Body Mass (ln(g))",
       y = "Baseline CCST (ln(ng/ml)") +
  annotate("text", x = 1, y = 6.5, label = "A", size = 6, 
           fontface = "bold", hjust = 0, vjust = 1) +
  add_phylopic(x = 2.5, y = 6.3, img = reptileimg, alpha = 1, height = 1.1)

ggsave("rep_BasalCCmass.png",
       width = 7,
       height = 5) 

#---BIRDS
#drop NAs
bird_BasalCCmass_data <- birdclean %>% 
  drop_na(c(Basal.CC, Body.Mass)) 
rownames(bird_BasalCCmass_data) = bird_BasalCCmass_data$Species.Name

#drop species in the tree that are not included in this relationship
if (sum(is.na(vertdata$Basal.CC)) > 0 | sum(is.na(vertdata$Body.Mass)) > 0) {
  bird_BasalCCmass_tree <- drop.tip(tree, name.check(tree, bird_BasalCCmass_data)$tree_not_data)
} else {
  bird_BasalCCmass_tree <- tree
}

#build pgls model 
bird_BasalCCmass_pgls <- gls(log(Basal.CC) ~ log(Body.Mass), 
                             data = bird_BasalCCmass_data, 
                             correlation = corPagel(value = 0, 
                                                    phy = bird_BasalCCmass_tree, 
                                                    form = ~Species.Name,
                                                    fixed = TRUE)) #lambda < 0 so manually bound at 0

#R2 values: specify reduced model as null model (no relationship, no phylogeny)
bird_BasalCCmass_R2_reduced <- lm(log(Basal.CC) ~ 1, 
                                  data = bird_BasalCCmass_data) 
bird_BasalCCmass_R2 <- R2(bird_BasalCCmass_pgls, bird_BasalCCmass_R2_reduced)

#plot
bird_BasalCCmass_plot <- ggplot(data = bird_BasalCCmass_data, aes(x = log(Body.Mass), y = log(Basal.CC))) + 
  geom_point(size = 3, color = "royalblue4") +
  geom_abline(intercept = coefficients(summary(bird_BasalCCmass_pgls))[1,1], 
              slope = coefficients(summary(bird_BasalCCmass_pgls))[2,1]) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14)) +
  xlim(1, 15.5) +
  ylim(-1, 6.5) +
  labs(x = "Body Mass (ln(g))",
       y = "Baseline CCST (ln(ng/ml)") +
  annotate("text", x = 1, y = 6.5, label = "B", size = 6, 
           fontface = "bold", hjust = 0, vjust = 1) +
  add_phylopic(x = 3.75, y = 6.3, img = birdimg, alpha = 1, height = 1)

ggsave("bird_BasalCCmass.png",
       width = 7,
       height = 5) 

#---MAMMALS
#drop NAs
mammal_BasalCCmass_data <- mammalclean %>% 
  drop_na(c(Basal.CC, Body.Mass)) 
rownames(mammal_BasalCCmass_data) = mammal_BasalCCmass_data$Species.Name

#drop species in the tree that are not included in this relationship
if (sum(is.na(vertdata$Basal.CC)) > 0 | sum(is.na(vertdata$Body.Mass)) > 0) {
  mammal_BasalCCmass_tree <- drop.tip(tree, name.check(tree, mammal_BasalCCmass_data)$tree_not_data)
} else {
  mammal_BasalCCmass_tree <- tree
}

#build pgls model 
mammal_BasalCCmass_pgls <- gls(log(Basal.CC) ~ log(Body.Mass), 
                               data = mammal_BasalCCmass_data, 
                               correlation = corPagel(value = 0.1, 
                                                      phy = mammal_BasalCCmass_tree, 
                                                      form = ~Species.Name)) 

#R2 values: specify reduced model as null model (no relationship, no phylogeny)
mammal_BasalCCmass_R2_reduced <- lm(log(Basal.CC) ~ 1, 
                                    data = mammal_BasalCCmass_data) 
mammal_BasalCCmass_R2 <- R2(mammal_BasalCCmass_pgls, mammal_BasalCCmass_R2_reduced)

#plot
mammal_BasalCCmass_plot <- ggplot(data = mammal_BasalCCmass_data, aes(x = log(Body.Mass), y = log(Basal.CC), shape = Dominant)) + 
  geom_point(size = 3, color = "violetred4") +
  geom_abline(intercept = coefficients(summary(mammal_BasalCCmass_pgls))[1,1], 
              slope = coefficients(summary(mammal_BasalCCmass_pgls))[2,1]) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14), 
        legend.position = "none") +
  xlim(1, 15.5) +
  ylim(-1, 6.6) +
  scale_shape_manual(values = c("D" = 1, "ND" = 16), labels = c("D" = "Dominant", "ND" = "Non-dominant")) +
  labs(x = "Body Mass (ln(g))",
       y = "Baseline CCST (ln(ng/ml)",
       shape = "Corticosterone") +
  annotate("text", x = 1, y = 6.5, label = "C", size = 6, 
           fontface = "bold", hjust = 0, vjust = 1) +
  add_phylopic(x = 6, y = 6.3, img = mammalimg, alpha = 1, height = 1)

ggsave("mammal_BasalCCmass.png",
       width = 7,
       height = 5) 


# Baseline CCST vs. MSMR --------------------------------------------------
#---REPTILES
#drop NAs
rep_CORTMSMR_data <- repclean %>% 
  drop_na(c(Basal.CC, MSMR)) 
rownames(rep_CORTMSMR_data) = rep_CORTMSMR_data$Species.Name

#drop species in the tree that are not included in this relationship
if (sum(is.na(vertdata$Basal.CC)) > 0 | sum(is.na(vertdata$MSMR)) > 0) {
  rep_CORTMSMR_tree <- drop.tip(tree, name.check(tree, rep_CORTMSMR_data)$tree_not_data)
} else {
  rep_CORTMSMR_tree <- tree
}

#build pgls model 
rep_CORTMSMR_pgls <- gls(log(Basal.CC) ~ log(MSMR), 
                         data = rep_CORTMSMR_data, 
                         correlation = corPagel(value = 0.1, 
                                                phy = rep_CORTMSMR_tree, 
                                                form = ~Species.Name)) 

#R2 values: specify reduced model as null model (no relationship, no phylogeny)
rep_CORTMSMR_R2_reduced <- lm(log(Basal.CC) ~ 1, 
                              data = rep_CORTMSMR_data) 
rep_CORTMSMR_R2 <- R2(rep_CORTMSMR_pgls, rep_CORTMSMR_R2_reduced)

#plot
rep_CORTMSMR_plot <- ggplot(data = rep_CORTMSMR_data, aes(x = log(MSMR), y = log(Basal.CC))) + 
  geom_point(size = 3, color = "green4") +
  geom_abline(intercept = coefficients(summary(rep_CORTMSMR_pgls))[1,1], 
              slope = coefficients(summary(rep_CORTMSMR_pgls))[2,1]) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14)) +
  xlim(-10, -3) +
  ylim(-2, 7) +
  labs(x = "MSMR (ln(W/g))",
       y = "Baseline CCST (ln[ng/ml])") +
  annotate("text", x = -10, y = 7, label = "A", size = 6, 
           fontface = "bold", hjust = 0, vjust = 1) +
  add_phylopic(x = -9, y = 6.8, img = reptileimg, alpha = 1, height = 1.1)

ggsave("rep_CORTMSMR.png",
       width = 7,
       height = 5) 

#---BIRDS
#drop NAs
bird_CORTMSMR_data <- birdclean %>% 
  drop_na(c(Basal.CC, MSMR)) 
rownames(bird_CORTMSMR_data) = bird_CORTMSMR_data$Species.Name

#drop species in the tree that are not included in this relationship
if (sum(is.na(vertdata$Basal.CC)) > 0 | sum(is.na(vertdata$MSMR)) > 0) {
  bird_CORTMSMR_tree <- drop.tip(tree, name.check(tree, bird_CORTMSMR_data)$tree_not_data)
} else {
  bird_CORTMSMR_tree <- tree
}

#build pgls model 
bird_CORTMSMR_pgls <- gls(log(Basal.CC) ~ log(MSMR), 
                          data = bird_CORTMSMR_data, 
                          correlation = corPagel(value = 0, 
                                                 phy = bird_CORTMSMR_tree, 
                                                 form = ~Species.Name,
                                                 fixed = TRUE)) #lambda < 0 so manually bound at 0

#R2 values: specify reduced model as null model (no relationship, no phylogeny)
bird_CORTMSMR_R2_reduced <- lm(log(Basal.CC) ~ 1, 
                               data = bird_CORTMSMR_data) 
bird_CORTMSMR_R2 <- R2(bird_CORTMSMR_pgls, bird_CORTMSMR_R2_reduced)

#plot
bird_CORTMSMR_plot <- ggplot(data = bird_CORTMSMR_data, aes(x = log(MSMR), y = log(Basal.CC))) + 
  geom_point(size = 3, color = "royalblue4") +
  geom_abline(intercept = coefficients(summary(bird_CORTMSMR_pgls))[1,1], 
              slope = coefficients(summary(bird_CORTMSMR_pgls))[2,1]) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14)) +
  xlim(-10, -3) +
  ylim(-2, 7) +
  labs(x = "MSMR (ln(W/g))",
       y = "Baseline CCST (ln[ng/ml])") +
  annotate("text", x = -10, y = 7, label = "B", size = 6, 
           fontface = "bold", hjust = 0, vjust = 1) +
  add_phylopic(x = -8.75, y = 6.8, img = birdimg, alpha = 1, height = 1)

ggsave("bird_CORTMSMR.png",
       width = 7,
       height = 5) 

#---MAMMALS
#drop NAs
mammal_CORTMSMR_data <- mammalclean %>% 
  drop_na(c(Basal.CC, MSMR)) 
rownames(mammal_CORTMSMR_data) = mammal_CORTMSMR_data$Species.Name

#drop species in the tree that are not included in this relationship
if (sum(is.na(vertdata$Basal.CC)) > 0 | sum(is.na(vertdata$MSMR)) > 0) {
  mammal_CORTMSMR_tree <- drop.tip(tree, name.check(tree, mammal_CORTMSMR_data)$tree_not_data)
} else {
  mammal_CORTMSMR_tree <- tree
}

#build pgls model 
mammal_CORTMSMR_pgls <- gls(log(Basal.CC) ~ log(MSMR), 
                            data = mammal_CORTMSMR_data, 
                            correlation = corPagel(value = 0.1,
                                                   phy = mammal_CORTMSMR_tree, 
                                                   form = ~Species.Name)) 

#R2 values: specify reduced model as null model (no relationship, no phylogeny)
mammal_CORTMSMR_R2_reduced <- lm(log(Basal.CC) ~ 1, 
                                 data = mammal_CORTMSMR_data) 
mammal_CORTMSMR_R2 <- R2(mammal_CORTMSMR_pgls, mammal_CORTMSMR_R2_reduced)

#plot
mammal_CORTMSMR_plot <- ggplot(data = mammal_CORTMSMR_data, aes(x = log(MSMR), y = log(Basal.CC), shape = Dominant)) + 
  geom_point(size = 3, color = "violetred4") +
  geom_abline(intercept = coefficients(summary(mammal_CORTMSMR_pgls))[1,1], 
              slope = coefficients(summary(mammal_CORTMSMR_pgls))[2,1]) +
  theme_classic() +
  scale_shape_manual(values = c("D" = 1, "ND" = 16), labels = c("D" = "Dominant", "ND" = "Non-dominant")) +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14), 
        legend.position = "none") +
  xlim(-10, -3) +
  ylim(-2, 7) +
  labs(x = "MSMR (ln(W/g))",
       y = "Baseline CCST (ln[ng/ml])",
       shape = "Corticosterone") +
  annotate("text", x = -10, y = 7, label = "C", size = 6, 
           fontface = "bold", hjust = 0, vjust = 1) +
  add_phylopic(x = -8.75, y = 6.8, img = mammalimg, alpha = 1, height = 1)

ggsave("mammal_CORTMSMR.png",
       width = 7,
       height = 5) 


# Lifespan Partial Regression ---------------------------------------------
#---REPTILES
#drop NAs
rep_LsCORT_data <- repclean %>% 
  drop_na(c(Basal.CC, Lifespan)) 
rownames(rep_LsCORT_data) = rep_LsCORT_data$Species.Name

#test metabolic rate residuals against stress concentrations
rep_LsCORT_lm <- lm(log(Lifespan) ~ log(MSMR), 
                    data = rep_LsCORT_data)

rep_LsCORT_residuals <- residuals(rep_LsCORT_lm)
rep_LsCORT_resmodel <- lm(rep_LsCORT_residuals ~ log(Basal.CC), 
                          data = rep_LsCORT_data %>%
                            drop_na(MSMR))
#---BIRDS
#drop NAs
bird_LsCORT_data <- birdclean %>% 
  drop_na(c(Basal.CC, Lifespan)) 
rownames(bird_LsCORT_data) = bird_LsCORT_data$Species.Name

#test metabolic rate residuals against stress concentrations
bird_LsCORT_lm <- lm(log(Lifespan) ~ log(MSMR), 
                     data = bird_LsCORT_data)

bird_LsCORT_residuals <- residuals(bird_LsCORT_lm)
bird_LsCORT_resmodel <- lm(bird_LsCORT_residuals ~ log(Basal.CC), 
                           data = bird_LsCORT_data %>%
                             drop_na(MSMR))
#---MAMMALS
#drop NAs
mammal_LsCORT_data <- mammalclean %>% 
  drop_na(c(Basal.CC, Lifespan)) 
rownames(mammal_LsCORT_data) = mammal_LsCORT_data$Species.Name

#test metabolic rate residuals against stress concentrations
mammal_LsCORT_lm <- lm(log(Lifespan) ~ log(MSMR), 
                       data = mammal_LsCORT_data)

mammal_LsCORT_residuals <- residuals(mammal_LsCORT_lm)
mammal_LsCORT_resmodel <- lm(mammal_LsCORT_residuals ~ log(Basal.CC), 
                             data = mammal_LsCORT_data %>%
                               drop_na(MSMR))


# Reptiles: Elevated vs. Baseline CCST ------------------------------------
#drop NAs
rep_ElevBasalCC_data <- repclean %>% 
  drop_na(c(Elevated.CC, Basal.CC))
rownames(rep_ElevBasalCC_data) = rep_ElevBasalCC_data$Species.Name

#drop species from tree that are not included in this relationship
if (sum(is.na(repclean$Elevated.CC)) > 0 | sum(is.na(repclean$Basal.CC)) > 0) {
  rep_ElevBasalCC_tree <- drop.tip(tree, name.check(tree, rep_ElevBasalCC_data)$tree_not_data)
} else {
  rep_ElevBasalCC_tree <- tree
}

#build pgls model 
rep_ElevBasalCC_pgls <- gls(log(Elevated.CC) ~ log(Basal.CC), 
                            data = rep_ElevBasalCC_data, 
                            correlation = corPagel(value = 0.1, 
                                                   phy = rep_ElevBasalCC_tree, 
                                                   form = ~Species.Name)) 
#R2 values
#specify reduced model as null (no relationship, no phylogeny)
rep_ElevBasalCC_R2_reduced <- lm(log(Elevated.CC) ~ 1, 
                                 data = rep_ElevBasalCC_data) 
rep_ElevBasalCC_R2 <- R2(rep_ElevBasalCC_pgls, rep_ElevBasalCC_R2_reduced)

#plot
rep_ElevBasalCC_plot <- ggplot(data = rep_ElevBasalCC_data, aes(x = log(Basal.CC), y = log(Elevated.CC))) + 
  geom_point(size = 2, color = "green4") +
  geom_abline(intercept = coefficients(summary(rep_ElevBasalCC_pgls))[1,1], 
              slope = coefficients(summary(rep_ElevBasalCC_pgls))[2,1]) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14)) +
  labs(x = "Baseline CCST (ln[ng/ml])",
       y = "Elevated CCST (ln[ng/ml])") +
  add_phylopic(x = 0.5, y = 6, img = reptileimg, alpha = 1, height = 1)

ggsave("rep_ElevBasalCC.png",
       width = 7,
       height = 5) 


# Mammals: CCST vs. Cortisol ----------------------------------------------
#drop NAs and CCSt-dominant species
mammal_GC_data <- mammalclean %>% 
  drop_na(c(Cortisol, Basal.CC)) %>%
  filter(Dominant == "ND")
rownames(mammal_GC_data) = mammal_GC_data$Species.Name

#drop species in the tree that are not included in this relationship
if (sum(is.na(mammalclean$Cortisol)) > 0 | sum(is.na(mammalclean$Basal.CC)) > 0) {
  mammal_GC_tree <- drop.tip(tree, name.check(tree, mammal_GC_data)$tree_not_data)
} else {
  mammal_GC_tree <- tree
}

#build pgls model 
mammal_GC_pgls <- gls(log(Basal.CC) ~ log(Cortisol), 
                            data = mammal_GC_data, 
                            correlation = corPagel(value = 0.1, 
                                                   phy = mammal_GC_tree, 
                                                   form = ~Species.Name)) 

#R2 values: specify reduced model as null model (no relationship, no phylogeny)
mammal_GC_R2_reduced <- lm(log(Basal.CC) ~ 1, 
                                 data = mammal_GC_data) 
mammal_GC_R2 <- R2(mammal_GC_pgls, mammal_GC_R2_reduced)

#plot
mammal_GC_plot <- ggplot(data = mammal_GC_data, aes(x = log(Cortisol), y = log(Basal.CC))) + 
  geom_point(size = 2, color = "violetred4") +
  geom_abline(intercept = coefficients(summary(mammal_GC_pgls))[1,1], 
              slope = coefficients(summary(mammal_GC_pgls))[2,1]) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14)) +
  xlim(0, 8) +
  ylim(-1, 5) +
  labs(x = "Baseline Cortisol (ln[ng/ml])",
       y = "Baseline CCST (ln[ng/ml])") +
  add_phylopic(x = 0.25, y = 4.9, img = mammalimg, alpha = 1, height = 0.75)

ggsave("mammal_GC.png",
       width = 7,
       height = 5) 

# Figures -----------------------------------------------------------------
#---Baseline CCST vs. Body Mass
rep_BasalCCmass_plot + bird_BasalCCmass_plot + mammal_BasalCCmass_plot + 
  plot_layout(axis_titles = "collect")

ggsave("BasalCCmass.png", 
       width = 13, 
       height = 6, 
       dpi = 200)

#---Baseline CCST vs. MSMR
rep_CORTMSMR_plot + bird_CORTMSMR_plot + mammal_CORTMSMR_plot + 
  plot_layout(axis_titles = "collect") 

ggsave("CORTMSMR.png", 
       width = 13, 
       height = 6,
       dpi = 200)

# Stats Table -------------------------------------------------------------
#spooky time 👻

#manually creating dataframe for all stats
#not sure if this is the best way to do this but we're doing it 
statstab <- data.frame(
  Model = c("Baseline CCST vs. Body Mass", "Baseline CCST vs. MSMR", "Elevated CCST vs. Baseline CCST", #reptiles
            "Baseline CCST vs. Body Mass", "Baseline CCST vs. MSMR", #birds
            "Baseline CCST vs. Body Mass", "Baseline CCST vs. MSMR"), #mammals
  n = c(nrow(rep_BasalCCmass_data), nrow(rep_CORTMSMR_data), nrow(rep_ElevBasalCC_data),
        nrow(bird_BasalCCmass_data), nrow(bird_CORTMSMR_data),
        nrow(mammal_BasalCCmass_data), nrow(mammal_CORTMSMR_data)),
  Slope = c(coefficients(rep_BasalCCmass_pgls)[2], coefficients(rep_CORTMSMR_pgls)[2], coefficients(rep_ElevBasalCC_pgls)[2],
            coefficients(bird_BasalCCmass_pgls)[2], coefficients(bird_CORTMSMR_pgls)[2],
            coefficients(mammal_BasalCCmass_pgls)[2], coefficients(mammal_CORTMSMR_pgls)[2]),
  Lower = c(intervals(rep_BasalCCmass_pgls)$coef[2,1], intervals(rep_CORTMSMR_pgls)$coef[2,1], intervals(rep_ElevBasalCC_pgls)$coef[2,1],
            intervals(bird_BasalCCmass_pgls)$coef[2,1], intervals(bird_CORTMSMR_pgls)$coef[2,1], 
            intervals(mammal_BasalCCmass_pgls)$coef[2,1], intervals(mammal_CORTMSMR_pgls)$coef[2,1]),
  Upper = c(intervals(rep_BasalCCmass_pgls)$coef[2,3], intervals(rep_CORTMSMR_pgls)$coef[2,3], intervals(rep_ElevBasalCC_pgls)$coef[2,3],
            intervals(bird_BasalCCmass_pgls)$coef[2,3], intervals(bird_CORTMSMR_pgls)$coef[2,3], 
            intervals(mammal_BasalCCmass_pgls)$coef[2,3], intervals(mammal_CORTMSMR_pgls)$coef[2,3]),
  Intercept = c(coefficients(rep_BasalCCmass_pgls)[1], coefficients(rep_CORTMSMR_pgls)[1], coefficients(rep_ElevBasalCC_pgls)[1],
                coefficients(bird_BasalCCmass_pgls)[1], coefficients(bird_CORTMSMR_pgls)[1],
                coefficients(mammal_BasalCCmass_pgls)[1], coefficients(mammal_CORTMSMR_pgls)[1]),
  R2.like = c(rep_BasalCCmass_R2[1], rep_CORTMSMR_R2[1], rep_ElevBasalCC_R2[1],
              bird_BasalCCmass_R2[1], bird_CORTMSMR_R2[1],
              mammal_BasalCCmass_R2[1], mammal_CORTMSMR_R2[1]),
  p = c(coefficients(summary(rep_BasalCCmass_pgls))[2,4], coefficients(summary(rep_CORTMSMR_pgls))[2,4], coefficients(summary(rep_ElevBasalCC_pgls))[2,4],
        coefficients(summary(bird_BasalCCmass_pgls))[2,4], coefficients(summary(bird_CORTMSMR_pgls))[2,4], 
        coefficients(summary(mammal_BasalCCmass_pgls))[2,4], coefficients(summary(mammal_CORTMSMR_pgls))[2,4]),
  Lambda = c(rep_BasalCCmass_pgls$modelStruct$corStruct[1], rep_CORTMSMR_pgls$modelStruct$corStruct[1], rep_ElevBasalCC_pgls$modelStruct$corStruct[1],
             bird_BasalCCmass_pgls$modelStruct$corStruct[1], bird_CORTMSMR_pgls$modelStruct$corStruct[1], 
             mammal_BasalCCmass_pgls$modelStruct$corStruct[1], mammal_CORTMSMR_pgls$modelStruct$corStruct[1]))

statstab <- statstab %>%
  mutate(n = as.character(as.integer(n))) %>%
  mutate(across(where(is.numeric), \(x) format(round(x, 2), nsmall = 2))) %>%
  mutate(p = case_when(p < 0.001 ~ '<0.001',
                       TRUE ~ as.character(p))) %>%
  mutate('Slope (95% CI)' = str_c(Slope, " (", Lower, ",", Upper, ")")) %>%
  select(Model, n, 'Slope (95% CI)', Intercept, R2.like, p, Lambda) 

#holy grail stats table 🙌
statstable <-
  gt(statstab, rowname_col = "Model") %>%
  cols_label(Lambda ~ "\U03BB",
             R2.like ~ "Likelihood R\U00B2") %>%
  tab_row_group(group = "Mammals", rows = 6:7) %>%
  tab_row_group(group = "Birds", rows = 4:5) %>%
  tab_row_group(group = "Reptiles", rows = 1:3) %>%
  tab_style(style = cell_text(weight = "bold"), 
            locations = cells_row_groups()) %>%
  tab_style(style = cell_text(align = "center", weight = "bold"), 
            locations = cells_column_labels(everything())) %>%
  tab_style(style = cell_text(align = "center"),
            locations = cells_body(everything())) %>%
  tab_style(style = cell_text(weight = "bold"),
            locations = cells_body(columns = p, rows = p <0.05)) %>%
  tab_options(table_body.border.bottom.width = px(2),
              data_row.padding.horizontal = px(25))

gtsave(statstable, filename = "Stats Table.png")




