
# Load packages and data --------------------------------------------------

library(tidyverse)

# Load the different traits
morph_traits = readRDS("Intermediate/Morphology/AllTraits.rds")

fitness_traits = readRDS("Intermediate/Fitness/AllTraits.rds")

photosynthesis_traits = readRDS("Intermediate/Photosynthesis/AllTraits.rds")

# Identify batches of data for random effect
morph_traits   = mutate(morph_traits, Batch = paste0(Genotype, XP, Block))

fitness_traits = mutate(fitness_traits, Batch = paste0(Genotype, Block))

# Select the specific timepoints and genotypes
morph_traits   = filter(morph_traits, (XP == "PSD" & Timepoint %in% c(3,5)) | 
                                      (XP == "HTD" & Timepoint == 3),
                                      Genotype != "Bur-0")

fitness_traits = filter(fitness_traits, Genotype != "Bur-0")

# Define the treatments as required by the linear models
morph_traits   = mutate(morph_traits,
                        Drought = ifelse(Treatment %in% c("D", "HTD", "PSD"), TRUE, FALSE),
                        MainTreatment = ifelse(Treatment %in% c("HT", "HTD"), "HT",
                                ifelse(Treatment %in% c("S","PSD"),
                                       ifelse(Timepoint == 3, "S", "PS"), "C")))

fitness_traits = filter(fitness_traits, !(Treatment %in% c("D1","D2"))) %>%
  mutate(Drought = ifelse(Treatment %in% c("D", "HTD", "PSD"), TRUE, FALSE),
         MainTreatment = ifelse(Treatment %in% c("HT", "HTD"), "HT",
                                ifelse(Treatment %in% c("S","PSD"), "PS", "C")))


convertFullTreatment = c(C = "C", D7 = "D", HT = "HT", HTD27 = "HTD",
                         PSD7 = "PSD", Su7 = "PS", Su1 = "S")

photosynthesis_traits = filter(photosynthesis_traits, MT == GT,
                               FullTreatment %in% c("C","D7","HT","HTD27","PSD7","Su7","Su1")) %>%
                          mutate(Treatment = convertFullTreatment[FullTreatment], 
                                 Drought = ifelse(Treatment %in% c("D", "HTD", "PSD"), TRUE, FALSE),
                                 MainTreatment = ifelse(Treatment %in% c("HT", "HTD"), "HT",
                                                        ifelse(Treatment %in% c("PS","PSD"), "PS", 
                                                               ifelse(Treatment == "S", "S", "C"))))

# Select the traits that we will actually use for the linear models
morph_traits = select(morph_traits, Genotype, MainTreatment, Drought, Batch, Timepoint, DW, N, Nmol,
                      RWC, RootLength, Angle, PetioleRatio, ShapeBlade, LeafSize, SLA)

fitness_traits = select(fitness_traits, Genotype, MainTreatment, Drought, Batch, Yield, 
                        Phyllochron, FT, LN)

photosynthesis_traits = select(photosynthesis_traits, MainTreatment, Drought, Batch, 
                               Chl, A, gsw, Rd, LUE)



# Save the data -----------------------------------------------------------

saveRDS(morph_traits, "Intermediate/Morphology/LinearModelTraits.rds")
saveRDS(fitness_traits, "Intermediate/Fitness/LinearModelTraits.rds")
saveRDS(photosynthesis_traits, "Intermediate/Photosynthesis/LinearModelTraits.rds")

