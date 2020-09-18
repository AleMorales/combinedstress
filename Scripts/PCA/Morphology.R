# Load packages and data --------------------------------------------------
library(tidyverse)

coefs = readRDS(file = "Intermediate/Morphology/LinearFits.rds")

# Calculate relative effects ----------------------------------------------

coefs = map(coefs, ~mutate(.x, ID = 1:n()))

# Combine into single table, remove the missing treatment combination 
morphology = do.call("rbind", coefs)

# Calculate morphology trait values per genotype & treatment -------------------
morphology = morphology %>%
  group_by(Trait) %>%
  mutate(# An-1
    An1_C   = exp(`GenotypeAn-1`),
    An1_HT  = exp(`GenotypeAn-1` + MainTreatmentHT),
    An1_S   = exp(`GenotypeAn-1` + MainTreatmentS),
    An1_PS  = exp(`GenotypeAn-1` + MainTreatmentPS),
    An1_D   = exp(`GenotypeAn-1` + DroughtTRUE),
    An1_HTD = exp(`GenotypeAn-1` + MainTreatmentHT + DroughtTRUE + `MainTreatmentHT:DroughtTRUE`),
    An1_PSD = exp(`GenotypeAn-1` + MainTreatmentPS + DroughtTRUE + `MainTreatmentPS:DroughtTRUE`),
    # Bay-0
    Bay0_C   = exp(`GenotypeBay-0`),
    Bay0_HT  = exp(`GenotypeBay-0` + MainTreatmentHT + `GenotypeBay-0:MainTreatmentHT`),
    Bay0_S   = exp(`GenotypeBay-0` + MainTreatmentS  + `GenotypeBay-0:MainTreatmentS`),
    Bay0_PS  = exp(`GenotypeBay-0` + MainTreatmentPS + `GenotypeBay-0:MainTreatmentPS`),
    Bay0_D   = exp(`GenotypeBay-0` + DroughtTRUE + `GenotypeBay-0:DroughtTRUE`),
    Bay0_HTD = exp(`GenotypeBay-0` + MainTreatmentHT + `GenotypeBay-0:MainTreatmentHT` +
                     DroughtTRUE + `GenotypeBay-0:DroughtTRUE` +
                     `MainTreatmentHT:DroughtTRUE` + `GenotypeBay-0:MainTreatmentHT:DroughtTRUE`),
    Bay0_PSD = exp(`GenotypeBay-0` + MainTreatmentPS + `GenotypeBay-0:MainTreatmentPS` +
                     DroughtTRUE + `GenotypeBay-0:DroughtTRUE` +
                     `MainTreatmentPS:DroughtTRUE` + `GenotypeBay-0:MainTreatmentPS:DroughtTRUE`),
    # Col-0
    Col0_C   = exp(`GenotypeCol-0`),
    Col0_HT  = exp(`GenotypeCol-0` + MainTreatmentHT + `GenotypeCol-0:MainTreatmentHT`),
    Col0_S   = exp(`GenotypeCol-0` + MainTreatmentS  + `GenotypeCol-0:MainTreatmentS`),
    Col0_PS  = exp(`GenotypeCol-0` + MainTreatmentPS + `GenotypeCol-0:MainTreatmentPS`),
    Col0_D   = exp(`GenotypeCol-0` + DroughtTRUE + `GenotypeCol-0:DroughtTRUE`),
    Col0_HTD = exp(`GenotypeCol-0` + MainTreatmentHT + `GenotypeCol-0:MainTreatmentHT` +
                     DroughtTRUE + `GenotypeCol-0:DroughtTRUE` +
                     `MainTreatmentHT:DroughtTRUE` + `GenotypeCol-0:MainTreatmentHT:DroughtTRUE`),
    Col0_PSD = exp(`GenotypeCol-0` + MainTreatmentPS + `GenotypeCol-0:MainTreatmentPS` +
                     DroughtTRUE + `GenotypeCol-0:DroughtTRUE` +
                     `MainTreatmentPS:DroughtTRUE` + `GenotypeCol-0:MainTreatmentPS:DroughtTRUE`),
    # Lp2-6
    Lp26_C   = exp(`GenotypeLp2-6`),
    Lp26_HT  = exp(`GenotypeLp2-6` + MainTreatmentHT + `GenotypeLp2-6:MainTreatmentHT`),
    Lp26_S   = exp(`GenotypeLp2-6` + MainTreatmentS  + `GenotypeLp2-6:MainTreatmentS`),
    Lp26_PS  = exp(`GenotypeLp2-6` + MainTreatmentPS + `GenotypeLp2-6:MainTreatmentPS`),
    Lp26_D   = exp(`GenotypeLp2-6` + DroughtTRUE + `GenotypeLp2-6:DroughtTRUE`),
    Lp26_HTD = exp(`GenotypeLp2-6` + MainTreatmentHT + `GenotypeLp2-6:MainTreatmentHT` +
                     DroughtTRUE + `GenotypeLp2-6:DroughtTRUE` +
                     `MainTreatmentHT:DroughtTRUE` + `GenotypeLp2-6:MainTreatmentHT:DroughtTRUE`),
    Lp26_PSD = exp(`GenotypeLp2-6` + MainTreatmentPS + `GenotypeLp2-6:MainTreatmentPS` +
                     DroughtTRUE + `GenotypeLp2-6:DroughtTRUE` +
                     `MainTreatmentPS:DroughtTRUE` + `GenotypeLp2-6:MainTreatmentPS:DroughtTRUE`)) %>%
  dplyr::select(Trait, ID, 
                An1_C,  An1_HT,  An1_S, An1_PS, An1_D, An1_HTD, An1_PSD,
                Bay0_C, Bay0_HT, Bay0_S, Bay0_PS, Bay0_D, Bay0_HTD, Bay0_PSD,
                Col0_C, Col0_HT, Col0_S, Col0_PS, Col0_D, Col0_HTD, Col0_PSD,
                Lp26_C, Lp26_HT, Lp26_S, Lp26_PS, Lp26_D, Lp26_HTD, Lp26_PSD)

# Create average of genotypes
morphology = mutate(morphology, 
                 Average_C   = (An1_C   + Bay0_C   + Col0_C   + Lp26_C)/4,
                 Average_HT  = (An1_HT  + Bay0_HT  + Col0_HT  + Lp26_HT)/4,
                 Average_S   = (An1_S   + Bay0_S   + Col0_S   + Lp26_S)/4,
                 Average_PS  = (An1_PS  + Bay0_PS  + Col0_PS  + Lp26_PS)/4,
                 Average_D   = (An1_D   + Bay0_D   + Col0_D   +  Lp26_D)/4,
                 Average_HTD = (An1_HTD + Bay0_HTD + Col0_HTD + Lp26_HTD)/4,
                 Average_PSD = (An1_PSD + Bay0_PSD + Col0_PSD + Lp26_PSD)/4)

# Split dataset accross genotypes, add genotype, rename and rbind it
morphology = map(c("An1", "Bay0", "Col0", "Lp26", "Average"), function(x) {
  out = select(morphology, Trait, ID, contains(x))
  names(out) = c("Trait", "ID", "C", "HT", "S", "PS", "D", "HTD", "PSD")
  mutate(out, Genotype = x)}) %>% 
  do.call("rbind", .)

# Reshape to long format
morphology = pivot_longer(morphology, c(-Trait, -ID, -Genotype), 
                       names_to = "Treatment", values_to = "Value")

# Summary statistics
sum_morphology = group_by(morphology, Genotype, Trait, Treatment) %>% 
  summarise(mu = median(Value))

saveRDS(sum_morphology, file = "Intermediate/Morphology/DataForPCA.rds")
saveRDS(morphology, file = "Intermediate/Morphology/DataForPCA_all.rds")
