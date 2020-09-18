# Load packages and data --------------------------------------------------
library(tidyverse)

coefs = readRDS("Intermediate/Fitness/LinearFits.rds")

# Calculate relative effects on fitness traits ----------------------------

# Combine into single table, remove the missing treatment combination 
fitness = do.call("rbind", coefs)

# Calculate the relative effects of each treatment for each genotype
relative_main_effect_An1 = function(effect, control, back = exp) {
  (back(effect + control) - back(control))/(back(control))
}
relative_interaction_An1 = function(effect1, effect2, interaction, control, back = exp) {
  (back(effect1 + effect2 + interaction + control) - 
     back(effect1 + effect2 + control))/(back(control))
}
relative_main_effect = function(effect, gen_effect, control, back = exp) {
  (back(effect + gen_effect + control) - back(control))/(back(control))
}
relative_interaction = function(effect1, gen_effect1, effect2, gen_effect2, interaction, gen_interaction, control, back = exp) {
  (back(effect1 + effect2 + interaction + gen_effect1 +gen_effect2 + gen_interaction + control) - 
     back(effect1 + effect2 + gen_effect1 +gen_effect2 + control))/(back(control))
}

fitness = fitness %>%
  group_by(Trait) %>%
  mutate(# An-1
    An1_HT   = relative_main_effect_An1(MainTreatmentHT, `GenotypeAn-1`),
    An1_PS   = relative_main_effect_An1(MainTreatmentPS, `GenotypeAn-1`),
    An1_D    = relative_main_effect_An1(DroughtTRUE, `GenotypeAn-1`),
    An1_HTD  = relative_main_effect_An1(MainTreatmentHT + DroughtTRUE + `MainTreatmentHT:DroughtTRUE`, `GenotypeAn-1`),
    An1_PSD  = relative_main_effect_An1(MainTreatmentPS + DroughtTRUE + `MainTreatmentPS:DroughtTRUE`, `GenotypeAn-1`),
    An1_HT_D = relative_interaction_An1(MainTreatmentHT,  DroughtTRUE, `MainTreatmentHT:DroughtTRUE`, `GenotypeAn-1`),
    An1_PS_D = relative_interaction_An1(MainTreatmentPS,  DroughtTRUE, `MainTreatmentPS:DroughtTRUE`, `GenotypeAn-1`),
    # Bay-0
    Bay0_HT   = relative_main_effect(MainTreatmentHT, `GenotypeBay-0:MainTreatmentHT`, `GenotypeBay-0`),
    Bay0_PS   = relative_main_effect(MainTreatmentPS, `GenotypeBay-0:MainTreatmentPS`, `GenotypeBay-0`),
    Bay0_D    = relative_main_effect(DroughtTRUE,     `GenotypeBay-0:DroughtTRUE`, `GenotypeBay-0`),
    Bay0_HTD = relative_main_effect(MainTreatmentHT + DroughtTRUE + `MainTreatmentHT:DroughtTRUE` +  
                                    `GenotypeBay-0:MainTreatmentHT` + `GenotypeBay-0:DroughtTRUE`, 
                                    `GenotypeBay-0:MainTreatmentHT:DroughtTRUE`, `GenotypeBay-0`),
    Bay0_PSD = relative_main_effect(MainTreatmentPS + DroughtTRUE + `MainTreatmentPS:DroughtTRUE` + 
                                    `GenotypeBay-0:MainTreatmentPS` + `GenotypeBay-0:DroughtTRUE`, 
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
    Col0_HTD = relative_main_effect(MainTreatmentHT + DroughtTRUE + `MainTreatmentHT:DroughtTRUE` +
                                      `GenotypeCol-0:MainTreatmentHT` + `GenotypeCol-0:DroughtTRUE`, 
                                    `GenotypeCol-0:MainTreatmentHT:DroughtTRUE`, `GenotypeCol-0`),
    Col0_PSD = relative_main_effect(MainTreatmentPS + DroughtTRUE + `MainTreatmentPS:DroughtTRUE` + 
                                      `GenotypeCol-0:MainTreatmentPS` + `GenotypeCol-0:DroughtTRUE`, 
                                    `GenotypeCol-0:MainTreatmentPS:DroughtTRUE`, `GenotypeCol-0`),
    Col0_HT_D = relative_interaction(MainTreatmentHT,  DroughtTRUE, `MainTreatmentHT:DroughtTRUE`, 
                                    `GenotypeCol-0:MainTreatmentHT`,  `GenotypeCol-0:DroughtTRUE`, 
                                    `GenotypeCol-0:MainTreatmentHT:DroughtTRUE`, `GenotypeCol-0`),
    Col0_PS_D = relative_interaction(MainTreatmentPS,  DroughtTRUE, `MainTreatmentPS:DroughtTRUE`, 
                                    `GenotypeCol-0:MainTreatmentPS`, `GenotypeCol-0:DroughtTRUE`, 
                                    `GenotypeCol-0:MainTreatmentPS:DroughtTRUE`, `GenotypeCol-0`),
    # Lp2-6
    Lp26_HT  = relative_main_effect(MainTreatmentHT, `GenotypeLp2-6:MainTreatmentHT`, `GenotypeLp2-6`),
    Lp26_PS  = relative_main_effect(MainTreatmentPS, `GenotypeLp2-6:MainTreatmentPS`, `GenotypeLp2-6`),
    Lp26_D   = relative_main_effect(DroughtTRUE,     `GenotypeLp2-6:DroughtTRUE`, `GenotypeLp2-6`),
    Lp26_HTD = relative_main_effect(MainTreatmentHT + DroughtTRUE + `MainTreatmentHT:DroughtTRUE` +  
                                      `GenotypeLp2-6:MainTreatmentHT` + `GenotypeLp2-6:DroughtTRUE`, 
                                    `GenotypeLp2-6:MainTreatmentHT:DroughtTRUE`, `GenotypeLp2-6`),
    Lp26_PSD = relative_main_effect(MainTreatmentPS + DroughtTRUE + `MainTreatmentPS:DroughtTRUE` +
                                      `GenotypeLp2-6:MainTreatmentPS` + `GenotypeLp2-6:DroughtTRUE`, 
                                    `GenotypeLp2-6:MainTreatmentPS:DroughtTRUE`, `GenotypeLp2-6`),
    Lp26_HT_D = relative_interaction(MainTreatmentHT,  DroughtTRUE, `MainTreatmentHT:DroughtTRUE`, 
                                    `GenotypeLp2-6:MainTreatmentHT`,  `GenotypeLp2-6:DroughtTRUE`, 
                                    `GenotypeLp2-6:MainTreatmentHT:DroughtTRUE`, `GenotypeLp2-6`),
    Lp26_PS_D = relative_interaction(MainTreatmentPS,  DroughtTRUE, `MainTreatmentPS:DroughtTRUE`, 
                                    `GenotypeLp2-6:MainTreatmentPS`, `GenotypeLp2-6:DroughtTRUE`, 
                                    `GenotypeLp2-6:MainTreatmentPS:DroughtTRUE`, `GenotypeLp2-6`)) %>%
  dplyr::select(Trait, 
                An1_HT, An1_PS, An1_D, An1_HTD, An1_PSD, An1_HT_D, An1_PS_D,
                Bay0_HT, Bay0_PS, Bay0_D, Bay0_HTD, Bay0_PSD, Bay0_HT_D, Bay0_PS_D,
                Col0_HT, Col0_PS, Col0_D, Col0_HTD, Col0_PSD, Col0_HT_D, Col0_PS_D,
                Lp26_HT, Lp26_PS, Lp26_D, Lp26_HTD, Lp26_PSD, Lp26_HT_D, Lp26_PS_D)

# Create average of genotypes
fitness = mutate(fitness, 
                 Average_HT = (An1_HT + Bay0_HT + Col0_HT + Lp26_HT)/4,
                 Average_PS = (An1_PS + Bay0_PS + Col0_PS + Lp26_PS)/4,
                 Average_D = (An1_D + Bay0_D + Col0_D + Lp26_D)/4,
                 Average_HTD = (An1_HTD + Bay0_HTD + Col0_HTD + Lp26_HTD)/4,
                 Average_PSD = (An1_PSD + Bay0_PSD + Col0_PSD + Lp26_PSD)/4,
                 Average_HT_D = (An1_HT_D + Bay0_HT_D + Col0_HT_D + Lp26_HT_D)/4,
                 Average_PS_D = (An1_PS_D + Bay0_PS_D + Col0_PS_D + Lp26_PS_D)/4)

# Split dataset accross genotypes, add genotype, rename and rbind it
fitness = map(c("An1", "Bay0", "Col0", "Lp26", "Average"), function(x) {
  out = select(fitness, Trait, contains(x))
  names(out) = c("Trait", "HT", "PS", "D", "HTD", "PSD", "HT_D", "PS_D")
  mutate(out, Genotype = x)}) %>% 
  do.call("rbind", .)

# Reshape to long format
fitness = pivot_longer(fitness, c(-Trait, -Genotype), 
                       names_to = "Effect", values_to = "Value")

# Summary statistics
sum_fitness = group_by(fitness, Genotype, Trait, Effect) %>% 
  summarise(mu = median(Value))

saveRDS(fitness, file = "Intermediate/Fitness/Effects.rds")

