
# 01 Load packages and data ---------------------------------------------------
library(tidyverse)


morphology     = readRDS("Intermediate/Morphology/LinearModelTraits.rds")
photosynthesis = readRDS("Intermediate/Photosynthesis/LinearModelTraits.rds")
fitness        = readRDS("Intermediate/Fitness/LinearModelTraits.rds")
rates          = readRDS("Intermediate/Rates/DataForLinearModel.rds")

# Calculate number of batches and number of replicates per batch ----------

# Angle
n_angle = morphology %>% 
            group_by(Genotype, MainTreatment, Drought) %>%
            filter(!is.na(Angle)) %>%
            summarise(nBatch = length(unique(Batch)),
                      nRepl  = n()/nBatch)

filter(n_angle, Genotype == "An-1")
filter(n_angle, Genotype == "Bay-0")
filter(n_angle, Genotype == "Col-0")
filter(n_angle, Genotype == "Lp2-6")

# DW
n_dw = morphology %>% 
  group_by(Genotype, MainTreatment, Drought) %>%
  filter(!is.na(DW)) %>%
  summarise(nBatch = length(unique(Batch)),
            nRepl  = n()/nBatch)

filter(n_dw, Genotype == "An-1")
filter(n_dw, Genotype == "Bay-0")
filter(n_dw, Genotype == "Col-0")
filter(n_dw, Genotype == "Lp2-6")

# LS
n_LS = morphology %>% 
  group_by(Genotype, MainTreatment, Drought) %>%
  filter(!is.na(LeafSize)) %>%
  summarise(nBatch = length(unique(Batch)),
            nRepl  = n()/nBatch)

filter(n_LS, Genotype == "An-1")
filter(n_LS, Genotype == "Bay-0")
filter(n_LS, Genotype == "Col-0")
filter(n_LS, Genotype == "Lp2-6")

# N
n_N = morphology %>% 
  group_by(Genotype, MainTreatment, Drought) %>%
  filter(!is.na(N)) %>%
  summarise(nBatch = length(unique(Batch)),
            nRepl  = n()/nBatch)

filter(n_N, Genotype == "An-1")
filter(n_N, Genotype == "Bay-0")
filter(n_N, Genotype == "Col-0")
filter(n_N, Genotype == "Lp2-6")

# Nmol
n_Nmol = morphology %>% 
  group_by(Genotype, MainTreatment, Drought) %>%
  filter(!is.na(Nmol)) %>%
  summarise(nBatch = length(unique(Batch)),
            nRepl  = n()/nBatch)

filter(n_Nmol, Genotype == "An-1")
filter(n_Nmol, Genotype == "Bay-0")
filter(n_Nmol, Genotype == "Col-0")
filter(n_Nmol, Genotype == "Lp2-6")

# PetioleRatio
n_PR = morphology %>% 
  group_by(Genotype, MainTreatment, Drought) %>%
  filter(!is.na(PetioleRatio)) %>%
  summarise(nBatch = length(unique(Batch)),
            nRepl  = n()/nBatch)

filter(n_PR, Genotype == "An-1")
filter(n_PR, Genotype == "Bay-0")
filter(n_PR, Genotype == "Col-0")
filter(n_PR, Genotype == "Lp2-6")

# RootLength
n_RL = morphology %>% 
  group_by(Genotype, MainTreatment, Drought) %>%
  filter(!is.na(RootLength)) %>%
  summarise(nBatch = length(unique(Batch)),
            nRepl  = n()/nBatch)

filter(n_RL, Genotype == "An-1")
filter(n_RL, Genotype == "Bay-0")
filter(n_RL, Genotype == "Col-0")
filter(n_RL, Genotype == "Lp2-6")

# ShapeBlade
n_SB = morphology %>% 
  group_by(Genotype, MainTreatment, Drought) %>%
  filter(!is.na(ShapeBlade)) %>%
  summarise(nBatch = length(unique(Batch)),
            nRepl  = n()/nBatch)

filter(n_SB, Genotype == "An-1")
filter(n_SB, Genotype == "Bay-0")
filter(n_SB, Genotype == "Col-0")
filter(n_SB, Genotype == "Lp2-6")

# SLA
n_SLA = morphology %>% 
  group_by(Genotype, MainTreatment, Drought) %>%
  filter(!is.na(SLA)) %>%
  summarise(nBatch = length(unique(Batch)),
            nRepl  = n()/nBatch)

filter(n_SLA, Genotype == "An-1")
filter(n_SLA, Genotype == "Bay-0")
filter(n_SLA, Genotype == "Col-0")
filter(n_SLA, Genotype == "Lp2-6")

# A
n_A = photosynthesis %>% 
  group_by(MainTreatment, Drought) %>%
  filter(!is.na(A)) %>%
  summarise(nBatch = length(unique(Batch)),
            nRepl  = n()/nBatch)
n_A

# Chl
n_Chl = photosynthesis %>% 
  group_by(MainTreatment, Drought) %>%
  filter(!is.na(Chl)) %>%
  summarise(nBatch = length(unique(Batch)),
            nRepl  = n()/nBatch)
n_Chl

# gsw
n_gsw = photosynthesis %>% 
  group_by(MainTreatment, Drought) %>%
  filter(!is.na(gsw)) %>%
  summarise(nBatch = length(unique(Batch)),
            nRepl  = n()/nBatch)
n_gsw

# Fitness
n_fitness = fitness %>% 
  group_by(Genotype, MainTreatment, Drought) %>%
  filter(!is.na(FT)) %>%
  summarise(nBatch = length(unique(Batch)),
            nRepl  = n()/nBatch)

filter(n_fitness, Genotype == "An-1")
filter(n_fitness, Genotype == "Bay-0")
filter(n_fitness, Genotype == "Col-0")
filter(n_fitness, Genotype == "Lp2-6")

# Phyllochron
n_phyllochron = fitness %>% 
  group_by(Genotype, MainTreatment, Drought) %>%
  filter(!is.na(Phyllochron)) %>%
  summarise(nBatch = length(unique(Batch)),
            nRepl  = n()/nBatch)

filter(n_phyllochron, Genotype == "An-1")
filter(n_phyllochron, Genotype == "Bay-0")
filter(n_phyllochron, Genotype == "Col-0")
filter(n_phyllochron, Genotype == "Lp2-6")

# Seed length

SeedSorting = readRDS("Intermediate/Biosorter/SeedSortingAll.rds") %>% 
  select(Block, Sample, Rep, Genotype, Treatment, XP, nSeeds, nNonSeeds, medSize, madSize) %>%
  rename(ID = Sample)

SeedSortingD = filter(SeedSorting, Treatment %in% c("D1", "D2")) %>% 
  mutate(Treatment = "D")
SeedSorting = rbind(SeedSorting, SeedSortingD) %>%
  filter(!(Treatment %in% c("D1", "D2")))

rm(SeedSortingD)

SeedSorting = SeedSorting %>%
  mutate(Drought = ifelse(Treatment %in% c("D", "HTD", "PSD"), TRUE, FALSE),
         MainTreatment = ifelse(Treatment %in% c("HT", "HTD"), "HT",
                                ifelse(Treatment %in% c("S","PSD"), "PS", "C")))

n_SL = SeedSorting %>% 
  group_by(Genotype, MainTreatment, Drought) %>%
  filter(!is.na(medSize)) %>%
  summarise(nBatch = 3,
            nRepl  = n()/nBatch)

filter(n_SL, Genotype == "An-1")
filter(n_SL, Genotype == "Bay-0")
filter(n_SL, Genotype == "Col-0")
filter(n_SL, Genotype == "Lp2-6")

# RGR
n_RGR = rates %>% 
  group_by(Genotype, Treatment) %>%
  filter(!is.na(DW)) %>%
  summarise(nBatch = length(unique(paste0(XP, Block))),
            nRepl  = n()/nBatch)

filter(n_RGR, Genotype == "An-1")
filter(n_RGR, Genotype == "Bay-0")
filter(n_RGR, Genotype == "Col-0")
filter(n_RGR, Genotype == "Lp2-6")

# LAR
n_LAR = rates %>% 
  group_by(Genotype, Treatment) %>%
  filter(!is.na(LeafNumber)) %>%
  summarise(nBatch = length(unique(paste0(XP, Block))),
            nRepl  = n()/nBatch)

filter(n_LAR, Genotype == "An-1")
filter(n_LAR, Genotype == "Bay-0")
filter(n_LAR, Genotype == "Col-0")
filter(n_LAR, Genotype == "Lp2-6")
