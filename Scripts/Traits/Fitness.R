# 01 Load packages and data -----------------------------------------------
library(tidyverse)
library(readxl)


# Fitness data ------------------------------------------------------------

# Load the fitness traits
Fitness = read_excel("Data/Fitness.xlsx") %>% 
  mutate(Treatment = c(C = "C", D1 = "D1", D2 = "D2", HtD = "HTD", Ht = "HT", 
                       FlD = "PSD", Flww = "S")[Treatment],
         TotalLeafRate = LN/(FT - 10), Phyllochron = 1/TotalLeafRate)

# Load the seed weights
SeedWeights = read_excel("Data/SeedWeights.xlsx", sheet = "RevisedSeedWeight") %>% 
  mutate(Treatment = c(C = "C", D1 = "D1", D2 = "D2", HtD = "HTD", Ht = "HT", 
                       FlD = "PSD", Flww = "S")[Treatment]) %>%
  filter(Rep == 1) %>% rename(Block = Experiment)

# Merge the two datasets
Fitness = full_join(Fitness, dplyr::select(SeedWeights, Genotype, Treatment, Block, ID, Yield))

# Reformat the treatments and do not distinguish between droughts
FitnessD = filter(Fitness, Treatment %in% c("D1", "D2")) %>% 
              mutate(Treatment = "D")
Fitness = rbind(Fitness, FitnessD)

# Save processed dataset --------------------------------------------------

saveRDS(Fitness, file = "Intermediate/Fitness/AllTraits.rds")



