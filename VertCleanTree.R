### Wisnionski et al. "Corticosterone Varies with Metabolic Rate in Vertebrates"
### Author: Alyson R. Wisnionski
### Incept Date: 6/2/2025

###---Building phylogenetic tree nexus file---###

library(tidyverse)

library(ape) #build tree
library(rotl) #pull from OTL
library(geiger) #name.check

#---data setup 
#laundry: load and clean data
repdata <- read.csv("Reptile_MASTER.csv", sep = ",", header = TRUE)
repclean <- repdata %>%
  rename(MSMR = MSMR.measured) %>%
  mutate(Vert = "Reptile") %>%
  select(Species.Name, Common.Name, Basal.CC)

birddata <- read.csv("Bird_MASTER.csv", sep = ",", header = TRUE)
birdclean <- birddata %>%
  mutate(Vert = "Bird") %>%
  select(Species.Name, Common.Name, Basal.CC) %>%
  filter(Species.Name != "Diomedea melanophris") #duplicate species

mammaldata <- read.csv("Mammal_MASTER.csv", sep = ",", header = TRUE)
mammalclean <- mammaldata %>%
  mutate(Vert = "Mammal") %>%
  select(Species.Name, Common.Name, Basal.CC)
  
#merge dataframes and rename species name changes in data to match tree
vertdata <- rbind(repclean, birdclean, mammalclean)
vertdata <- vertdata %>%
  mutate(Species.Name = recode(Species.Name, 
                               "Cnemidophorus sexlineatus" = "Aspidoscelis sexlineatus",
                               "Cnemidophorus uniparens" = "Aspidoscelis uniparens",
                               "Crotalus helleri" = "Crotalus viridis",
                               "Egernia whitii" = "Liopholis whitii",
                               "Lacerta vivipara" = "Zootoca vivipara",
                               "Podarcis sicula" = "Podarcis siculus",
                               "Thamnophis elegans vagrans" = "Thamnophis errans",
                               "Carduelis flammea" = "Acanthis flammea",
                               "Gallus domesticus" = "Gallus gallus",
                               "Grus canadensis" = "Antigone canadensis",
                               "Melozone aberti" = "Kieneria aberti",
                               "Parus caeruleus" = "Cyanistes caeruleus",
                               "Peucaea carpalis" = "Aimophila carpalis",
                               "Rhynchophanes mccownii" = "Calcarius mccownii",
                               "Spizella arborea" = "Spizelloides arborea",
                               "Thryesphilus rufalbus" = "Thryophilus rufalbus", 
                               "Thallasarche melanophris" = "Thalassarche melanophris",
                               "Ticiqua scincoides" = "Tiliqua scincoides",
                               "Macropus rufogriseus" = "Notamacropus rufogriseus",
                               "Perognathus baileyi" = "Chaetodipus baileyi"))

#---build tree
taxa <- tnrs_match_names(unique(vertdata$Species.Name)) #matching data species with OTL 
tree <- tol_induced_subtree(ott_ids = taxa$ott_id) #pull the subtree of the matched names 
tree <- compute.brlen(tree, method = "Grafen", power = 1) #compute branch lengths using Grafen method
tree$tip.label <- strip_ott_ids(tree$tip.label, remove_underscores = T) #strip ott ids

#rename tree labels
tree$tip.label[tree$tip.label == "Crocodylus johnsoni"] <- "Crocodylus johnstoni"
tree$tip.label[tree$tip.label == "Thalassarche melanophrys"] <- "Thalassarche melanophris"

cbind(sort(tree$tip.label), sort(unique(vertdata$Species.Name))) #view dataframe species and tip labels lines up

rownames(vertdata) = vertdata$Species.Name #check that tree and species from dataframe match
name.check(tree, vertdata) #OK

#saving tree 
write.nexus(tree, file = "VertTree.nex")

png("TreePic.png",
    height = 10000,
    width = 6000,
    res = 300)
plot(tree)
dev.off()
