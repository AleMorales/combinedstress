# Load packages and data ----------------------------------------------------
library(tidyverse)
library(readxl)


# Read the traits data -----------------------------------------------------

# Load the morphological traits
HarvestPSD = read_excel("Data/MorphologyTraitsFLD.xlsx", sheet = "Harvest", na = c("NA",""))
PlantPSD = read_excel("Data/MorphologyTraitsFLD.xlsx",   sheet = "Plant", na = c("NA","")) %>%
  mutate(FW = ifelse(FW < DW, FW*1000, FW))
OrganPSD = read_excel("Data/MorphologyTraitsFLD.xlsx",   sheet = "Organ", na = c("NA",""))
HarvestHTD = read_excel("Data/MorphologyTraitsHTD.xlsx", sheet = "Harvest", na = c("NA",""))
PlantHTD = read_excel("Data/MorphologyTraitsHTD.xlsx",   sheet = "Plant", na = c("NA","")) %>%
  mutate(FW = ifelse(FW < DW, FW*1000, FW))
OrganHTD = read_excel("Data/MorphologyTraitsHTD.xlsx",   sheet = "Organ", na = c("NA",""))



# Merge different experiments ---------------------------------------------------

# Use better genotype labels for figures
gencodes = c(Col = "Col-0", An = "An-1", Bur = "Bur-0", Bay = "Bay-0", LP2 = "Lp2-6")
HarvestPSD = mutate(HarvestPSD, Genotype = gencodes[Genotype])
HarvestHTD = mutate(HarvestHTD, Genotype = gencodes[Genotype])

# Add label to each table to facilitate better merging
HarvestPSD = mutate(HarvestPSD, XP = "PSD")
PlantPSD = mutate(PlantPSD, XP = "PSD")
OrganPSD = mutate(OrganPSD, XP = "PSD")
HarvestHTD = mutate(HarvestHTD, XP = "HTD")
PlantHTD = mutate(PlantHTD, XP = "HTD")
OrganHTD = mutate(OrganHTD, XP = "HTD")

# Merge the two experiments
Harvest = rbind(HarvestPSD, HarvestHTD)
Plant = rbind(PlantPSD, PlantHTD)
Organ = rbind(OrganPSD, OrganHTD)

# Cleanup
rm(HarvestHTD, HarvestPSD, PlantHTD, PlantPSD, OrganHTD, OrganPSD)


# Remove pilot experiments ------------------------------------------------

# HTD pilot experiments
removedID = filter(Harvest, FileID == "SaskiaBio1")$HarvestID
Harvest = filter(Harvest, FileID != "SaskiaBio1")
Plant = filter(Plant, !(HarvestID %in% removedID))
Organ = filter(Organ, !(HarvestID %in% removedID))

# PSD pilot experiments
removedID = filter(Harvest, FileID == "Kevin10")$HarvestID
Harvest = filter(Harvest, FileID != "Kevin10")
Plant = filter(Plant, !(HarvestID %in% removedID))
Organ = filter(Organ, !(HarvestID %in% removedID))


# Derive additional traits ------------------------------------------------

Organ_Plant = group_by(Organ, XP, HarvestID, Replicate) %>%
  summarise(LeafNumber = sum(BladeWidth > 0, na.rm = TRUE),
            PlantArea = sum(LeafArea, na.rm = TRUE),
            PetioleLength = mean(PetioleLength, na.rm = TRUE),
            BladeLength = mean(BladeLength, na.rm = TRUE),
            PetioleRatio = mean(PetioleLength/(BladeLength + PetioleLength), na.rm = TRUE),
            ShapeBlade = mean(BladeWidth/BladeLength, na.rm = TRUE)) %>%
  dplyr::select(XP, HarvestID, Replicate, LeafNumber, PlantArea, PetioleLength, BladeLength, PetioleRatio, ShapeBlade)

Plant = mutate(Plant, RWC = (FW - DW)/FW)

MorphTraits = full_join(Plant, Organ_Plant) %>% left_join(Harvest)
rm(Organ_Plant, Organ, Harvest, Plant)

MorphTraits = mutate(MorphTraits, 
                     Phyllochron = (DAS - 10)/LeafNumber,
                     LeafSize = PlantArea/LeafNumber,
                     LeafWeight = DW/LeafNumber,
                     SLA = PlantArea/DW*10, # cm2/g
                     RelativeAngle = Angle/90,
                     N = NA,
                     Treatment = paste0(Temperature, "_", Water),
                     Treatment = c(`21_C` = "C", `27_C` = "HT", `21_D` = "D",
                                   `27_D` = "HTD", `21_S` = "S", 
                                   `21_PSD` = "PSD")[Treatment])


# Add extra HTD experiment ------------------------------------------------

ExtraHTD = read_csv("Data/ChenYunHTD.csv")
ExtraHTD = mutate(ExtraHTD, 
                  RootLength = NA,
                  HypocotylLength = NA,
                  Angle = NA,
                  Bolting = NA,
                  XP = "HTD",
                  PetioleLength = NA,
                  BladeLength = NA,
                  PetioleRatio = NA,
                  ShapeBlade = NA,
                  Timepoint = 3,
                  FileID = NA,
                  Phyllochron = (DAS - 10)/LeafNumber,
                  LeafSize = PlantArea/LeafNumber,
                  LeafWeight = DW/LeafNumber,
                  RelativeAngle = NA,
                  RWC = (FW - DW)/FW,
                  HarvestID = NA,
                  SLA = PlantArea/DW*10, # cm2/g
                  Treatment = paste0(Temperature, "_", Water),
                  Treatment = c(`21_C` = "C", `27_C` = "HT", `21_D` = "D",
                                `27_D` = "HTD", `21_S` = "S", 
                                `21_PSD` = "PSD")[Treatment])


MorphTraits = rbind(MorphTraits, ExtraHTD)


# Add extra PSD experiment ------------------------------------------------

ExtraPSD = read_csv("Data/FloodingCN.csv")
ExtraPSD = mutate(ExtraPSD, 
                  RootLength = NA,
                  HypocotylLength = NA,
                  Angle = NA,
                  Bolting = NA,
                  XP = "PSD",
                  PetioleLength = NA,
                  BladeLength = NA,
                  PetioleRatio = NA,
                  ShapeBlade = NA,
                  FileID = NA,
                  Phyllochron = NA,
                  LeafSize = NA,
                  LeafWeight = NA,
                  RelativeAngle = NA,
                  DAS = NA,
                  LeafNumber = NA,
                  Temperature = 21,
                  HarvestID = NA,
                  Water = Treatment) %>% 
  rename(PlantArea = Area)

MorphTraits = rbind(MorphTraits, ExtraPSD)



# Further small adjustments -----------------------------------------------

MorphTraits = mutate(MorphTraits, Nmol = (N/100)/(SLA/1e4)/14, # mol N/m2
                     PlantArea = PlantArea/1e2, # cm2
                     LeafSize = LeafSize/1e2) # cm2

MorphTraits = filter(MorphTraits, Genotype != "Bur-0")

saveRDS(MorphTraits, file = "Intermediate/Morphology/AllTraits.rds")

