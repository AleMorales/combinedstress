# Load packages and data --------------------------------------------------
library(tidyverse)
library(robustbase)
library(furrr)

plan(multiprocess)

# Load the seed sorting traits
if(!file.exists("Intermediate/Biosorter/SeedSortingAll.rds")) 
  source("Scripts/BioSorter/00_SeedSorting.R")
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

# Check number of samples per genotype
with(SeedSorting, table(Genotype, Treatment))

# Remove Lp2-6 and Bur-0
SeedSorting = filter(SeedSorting, !(Genotype %in% c("Lp2-6", "Bur-0")),
                           !is.na(medSize))

# Check number of samples per genotype
with(SeedSorting, table(Genotype, MainTreatment, Drought))

# Calculate robust linear model -------------------------------------------

fit_medSize = lmrob(medSize~Genotype*MainTreatment*Drought - 1, data = SeedSorting,
                    control = lmrob.control(k.max = 1000, setting="KS2014",
                    max.it = 1000, maxit.scale = 1e3, refine.tol = 1e-5))
summary(fit_medSize)
#plot(fit_medSize)

refit = function() {
  new_data = mutate(SeedSorting, medSize = fitted(fit_medSize) + 
                      sample(residuals(fit_medSize), replace = TRUE))
  new_fit = update(fit_medSize, data = new_data)
  coef(new_fit)
}

N = 5e3
medSize_coef = future_map_dfr(1:N, ~refit(), .progress = T, 
                         .options = future_options(packages = "robustbase"))


saveRDS(medSize_coef, file = "Intermediate/Biosorter/LinearModel.rds")

# Calculate values for each genotype and treatment combination ------------

medSize_coef = readRDS(file = "Intermediate/Biosorter/LinearModel.rds")

# Calculate the relative effects of each treatment for each genotype
relative_main_effect_An1 = function(effect, control) {
  effect/control
}
relative_interaction_An1 = function(effect1, effect2, interaction, control) {
  interaction/control
}
relative_main_effect = function(effect, gen_effect, control) {
  (effect + gen_effect)/control
}
relative_interaction = function(effect1, gen_effect1, effect2, gen_effect2, interaction, 
                                gen_interaction, control) {
  (interaction + gen_interaction)/control
}
medSize = medSize_coef %>%
  mutate(# An-1
    An1_HT  = relative_main_effect_An1(MainTreatmentHT, `GenotypeAn-1`),
    An1_PS  = relative_main_effect_An1(MainTreatmentPS, `GenotypeAn-1`),
    An1_D   = relative_main_effect_An1(DroughtTRUE, `GenotypeAn-1`),
    An1_HTD = relative_main_effect_An1(MainTreatmentHT + DroughtTRUE +
                                       `MainTreatmentHT:DroughtTRUE`, `GenotypeAn-1`),
    An1_PSD = relative_main_effect_An1(MainTreatmentPS + DroughtTRUE +
                                       `MainTreatmentPS:DroughtTRUE`, `GenotypeAn-1`),
    An1_HT_D = relative_interaction_An1(MainTreatmentHT,  DroughtTRUE, 
                                       `MainTreatmentHT:DroughtTRUE`, `GenotypeAn-1`),
    An1_PS_D = relative_interaction_An1(MainTreatmentPS,  DroughtTRUE, 
                                       `MainTreatmentPS:DroughtTRUE`, `GenotypeAn-1`),
    # Bay-0
    Bay0_HT  = relative_main_effect(MainTreatmentHT, `GenotypeBay-0:MainTreatmentHT`, `GenotypeBay-0`),
    Bay0_PS  = relative_main_effect(MainTreatmentPS, `GenotypeBay-0:MainTreatmentPS`, `GenotypeBay-0`),
    Bay0_D   = relative_main_effect(DroughtTRUE,     `GenotypeBay-0:DroughtTRUE`, `GenotypeBay-0`),
    
    Bay0_HTD = relative_main_effect(MainTreatmentHT + DroughtTRUE + `MainTreatmentHT:DroughtTRUE`, 
                                     `GenotypeBay-0:MainTreatmentHT` + `GenotypeBay-0:DroughtTRUE` +
                                     `GenotypeBay-0:MainTreatmentHT:DroughtTRUE`, `GenotypeBay-0`),
    Bay0_PSD = relative_main_effect(MainTreatmentPS + DroughtTRUE + `MainTreatmentPS:DroughtTRUE`, 
                                     `GenotypeBay-0:MainTreatmentPS` + `GenotypeBay-0:DroughtTRUE` +
                                     `GenotypeBay-0:MainTreatmentPS:DroughtTRUE`, `GenotypeBay-0`),
    Bay0_HT_D = relative_interaction(MainTreatmentHT,  DroughtTRUE, `MainTreatmentHT:DroughtTRUE`, 
                                    `GenotypeBay-0:MainTreatmentHT`,  `GenotypeBay-0:DroughtTRUE`, 
                                    `GenotypeBay-0:MainTreatmentHT:DroughtTRUE`, `GenotypeBay-0`),
    Bay0_PS_D = relative_interaction(MainTreatmentPS,  DroughtTRUE, `MainTreatmentPS:DroughtTRUE`, 
                                    `GenotypeBay-0:MainTreatmentPS`, `GenotypeBay-0:DroughtTRUE`, 
                                    `GenotypeBay-0:MainTreatmentPS:DroughtTRUE`, `GenotypeBay-0`),
    # Col-0
    Col0_HT  = relative_main_effect(MainTreatmentHT, `GenotypeCol-0:MainTreatmentHT`, `GenotypeCol-0`),
    Col0_PS  = relative_main_effect(MainTreatmentPS, `GenotypeCol-0:MainTreatmentPS`, `GenotypeCol-0`),
    Col0_D   = relative_main_effect(DroughtTRUE,     `GenotypeCol-0:DroughtTRUE`, `GenotypeCol-0`),
    Col0_HTD = relative_main_effect(MainTreatmentHT + DroughtTRUE + `MainTreatmentHT:DroughtTRUE`, 
                                    `GenotypeCol-0:MainTreatmentHT` + `GenotypeCol-0:DroughtTRUE` +
                                      `GenotypeCol-0:MainTreatmentHT:DroughtTRUE`, `GenotypeCol-0`),
    
    Col0_PSD = relative_main_effect(MainTreatmentPS + DroughtTRUE + `MainTreatmentPS:DroughtTRUE`, 
                                    `GenotypeCol-0:MainTreatmentPS` + `GenotypeCol-0:DroughtTRUE` +
                                      `GenotypeCol-0:MainTreatmentPS:DroughtTRUE`, `GenotypeCol-0`),
    Col0_HT_D = relative_interaction(MainTreatmentHT,  DroughtTRUE, `MainTreatmentHT:DroughtTRUE`, 
                                    `GenotypeCol-0:MainTreatmentHT`,  `GenotypeCol-0:DroughtTRUE`, 
                                    `GenotypeCol-0:MainTreatmentHT:DroughtTRUE`, `GenotypeCol-0`),
    Col0_PS_D = relative_interaction(MainTreatmentPS,  DroughtTRUE, `MainTreatmentPS:DroughtTRUE`, 
                                    `GenotypeCol-0:MainTreatmentPS`, `GenotypeCol-0:DroughtTRUE`, 
                                    `GenotypeCol-0:MainTreatmentPS:DroughtTRUE`, `GenotypeCol-0`)) %>%
  dplyr::select(An1_HT,  An1_PS, An1_D, An1_HTD, An1_PSD, An1_HT_D, An1_PS_D,
                Bay0_HT, Bay0_PS, Bay0_D, Bay0_HTD, Bay0_PSD, Bay0_HT_D, Bay0_PS_D,
                Col0_HT, Col0_PS, Col0_D, Col0_HTD, Col0_PSD, Col0_HT_D, Col0_PS_D)

# Create average of genotypes
medSize = mutate(medSize, 
                 Average_HT  = (An1_HT  + Bay0_HT  + Col0_HT)/3,
                 Average_PS  = (An1_PS  + Bay0_PS  + Col0_PS)/3,
                 Average_D   = (An1_D   + Bay0_D   + Col0_D)/3,
                 Average_HTD = (An1_HTD + Bay0_HTD + Col0_HTD)/3,
                 Average_PSD = (An1_PSD + Bay0_PSD + Col0_PSD)/3,
                 Average_HT_D = (An1_HT_D + Bay0_HT_D + Col0_HT_D)/3,
                 Average_PS_D = (An1_PS_D + Bay0_PS_D + Col0_PS_D)/3)

# Split dataset across genotypes, add genotype, rename and rbind it
medSize = map(c("An1", "Bay0", "Col0", "Average"), function(x) {
  out = select(medSize,  contains(x))
  names(out) = c("HT", "PS", "D", "HTD", "PSD", "HT_D", "PS_D")
  mutate(out, Genotype = x)}) %>% 
  do.call("rbind", .)

medSize =  pivot_longer(medSize, -Genotype, names_to = "Coef",
                        values_to = "Value") %>%
            mutate(mu = NA, Trait = "SeedSize") %>% 
            rename(Effect = Coef)


saveRDS(medSize, file = "Intermediate/Biosorter/Effects.rds")

# Calculate average values for each combination and treatment -------------

medSize_coef = readRDS(file = "Intermediate/Biosorter/LinearModel.rds")  %>%
                  mutate(ID = 1:n())

medSize = medSize_coef %>%
  mutate(# An-1
    An1_C   = `GenotypeAn-1`,
    An1_HT  = MainTreatmentHT + `GenotypeAn-1`,
    An1_PS  = MainTreatmentPS + `GenotypeAn-1`,
    An1_D   = DroughtTRUE + `GenotypeAn-1`,
    An1_HTD = MainTreatmentHT + DroughtTRUE + 
              `MainTreatmentHT:DroughtTRUE` + `GenotypeAn-1`,
    An1_PSD = MainTreatmentPS + DroughtTRUE + 
              `MainTreatmentPS:DroughtTRUE` + `GenotypeAn-1`,
    # Bay-0
    Bay0_C   = `GenotypeBay-0`,
    Bay0_HT  = MainTreatmentHT + `GenotypeBay-0:MainTreatmentHT` + `GenotypeBay-0`,
    Bay0_PS  = MainTreatmentPS + `GenotypeBay-0:MainTreatmentPS` + `GenotypeBay-0`,
    Bay0_D   = DroughtTRUE + `GenotypeBay-0:DroughtTRUE` + `GenotypeBay-0`,
    Bay0_HTD = MainTreatmentHT + DroughtTRUE + `MainTreatmentHT:DroughtTRUE` + 
               `GenotypeBay-0:MainTreatmentHT` + `GenotypeBay-0:DroughtTRUE` + 
               `GenotypeBay-0:MainTreatmentHT:DroughtTRUE` + `GenotypeBay-0`,
    Bay0_PSD = MainTreatmentPS + DroughtTRUE + `MainTreatmentPS:DroughtTRUE` + 
               `GenotypeBay-0:MainTreatmentPS` + `GenotypeBay-0:DroughtTRUE` + 
               `GenotypeBay-0:MainTreatmentPS:DroughtTRUE` + `GenotypeBay-0`,
    # Col-0
    Col0_C   = `GenotypeCol-0`,
    Col0_HT  = MainTreatmentHT + `GenotypeCol-0:MainTreatmentHT` + `GenotypeCol-0`,
    Col0_PS  = MainTreatmentPS + `GenotypeCol-0:MainTreatmentPS` + `GenotypeCol-0`,
    Col0_D   = DroughtTRUE + `GenotypeCol-0:DroughtTRUE` + `GenotypeCol-0`,
    Col0_HTD = MainTreatmentHT + DroughtTRUE + `MainTreatmentHT:DroughtTRUE` + 
               `GenotypeCol-0:MainTreatmentHT` + `GenotypeCol-0:DroughtTRUE` + 
               `GenotypeCol-0:MainTreatmentHT:DroughtTRUE` + `GenotypeCol-0`,
    Col0_PSD = MainTreatmentPS + DroughtTRUE + `MainTreatmentPS:DroughtTRUE` + 
               `GenotypeCol-0:MainTreatmentPS` + `GenotypeCol-0:DroughtTRUE` + 
               `GenotypeCol-0:MainTreatmentPS:DroughtTRUE` + `GenotypeCol-0`) %>%
  dplyr::select(ID,
                An1_C,   An1_HT,  An1_PS,  An1_D,  An1_HTD,  An1_PSD,
                Bay0_C, Bay0_HT, Bay0_PS, Bay0_D, Bay0_HTD, Bay0_PSD,
                Col0_C, Col0_HT, Col0_PS, Col0_D, Col0_HTD, Col0_PSD)

# Create average of genotypes
medSize = mutate(medSize, 
                 Average_C = (An1_C + Bay0_C + Col0_C)/3,
                 Average_HT = (An1_HT + Bay0_HT + Col0_HT)/3,
                 Average_PS = (An1_PS + Bay0_PS + Col0_PS)/3,
                 Average_D = (An1_D + Bay0_D + Col0_D)/3,
                 Average_HTD = (An1_HTD + Bay0_HTD + Col0_HTD)/3,
                 Average_PSD = (An1_PSD + Bay0_PSD + Col0_PSD)/3)

# Split dataset accross genotypes, add genotype, rename and rbind it
medSize = map(c("An1", "Bay0", "Col0", "Average"), function(x) {
            out = select(medSize,  ID, contains(x))
            names(out) = c("ID", "C","HT", "PS", "D", "HTD", "PSD")
            mutate(out, Genotype = x)}) %>% 
            do.call("rbind", .)

# Reshape to longest format
medSize = pivot_longer(medSize, c(-Genotype, -ID), 
                       names_to = "Effect", values_to = "SeedLength")

# Save results ------------------------------------------------------------

saveRDS(medSize, file = "Intermediate/Biosorter/DataForCorrelations.rds")


# Calculate mean and stdev of each trait ----------------------------------


SeedLength = readRDS(file = "Intermediate/Biosorter/DataForCorrelations.rds")

sum_SeedLength = group_by(SeedLength, Genotype, Effect) %>%
  summarise(muSeedLength = median(SeedLength), seSeedLength = sd(SeedLength))

