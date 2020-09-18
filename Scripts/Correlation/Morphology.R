
# Compute correlations among traits
# Load packages and data --------------------------------------------------
library(tidyverse)

# morphology traits
coefs = readRDS(file = "Intermediate/Morphology/LinearFits.rds")

# Reshape the data as required --------------------------------------------
coefs = map(coefs, ~mutate(.x, ID = 1:n()))

# Combine into single table, remove the missing treatment combination 
morphology = do.call("rbind", coefs) %>% filter(Trait != "RWC")

morphology = morphology %>%
  group_by(Trait) %>%
  mutate(# An-1
    An1_C  = exp(`GenotypeAn-1`),
    An1_HT  = exp(MainTreatmentHT + `GenotypeAn-1`),
    An1_S   = exp(MainTreatmentS + `GenotypeAn-1`),
    An1_PS  = exp(MainTreatmentPS + `GenotypeAn-1`),
    An1_D   = exp(DroughtTRUE + `GenotypeAn-1`),
    An1_HTD = exp(MainTreatmentHT + DroughtTRUE + `MainTreatmentHT:DroughtTRUE` + `GenotypeAn-1`),
    An1_PSD = exp(MainTreatmentPS + DroughtTRUE + `MainTreatmentPS:DroughtTRUE` + `GenotypeAn-1`),
    # Bay-0
    Bay0_C   = exp(`GenotypeBay-0`),
    Bay0_HT  = exp(MainTreatmentHT + `GenotypeBay-0:MainTreatmentHT` + `GenotypeBay-0`),
    Bay0_S   = exp(MainTreatmentS  + `GenotypeBay-0:MainTreatmentS`  + `GenotypeBay-0`),
    Bay0_PS  = exp(MainTreatmentPS + `GenotypeBay-0:MainTreatmentPS` + `GenotypeBay-0`),
    Bay0_D   = exp(DroughtTRUE + `GenotypeBay-0:DroughtTRUE` + `GenotypeBay-0`),
    Bay0_HTD = exp(MainTreatmentHT + DroughtTRUE + `MainTreatmentHT:DroughtTRUE` + 
                     `GenotypeBay-0:MainTreatmentHT` + `GenotypeBay-0:DroughtTRUE` +
                     `GenotypeBay-0:MainTreatmentHT:DroughtTRUE` + `GenotypeBay-0`),
    Bay0_PSD = exp(MainTreatmentPS + DroughtTRUE + `MainTreatmentPS:DroughtTRUE` +
                     `GenotypeBay-0:MainTreatmentPS` + `GenotypeBay-0:DroughtTRUE` + 
                     `GenotypeBay-0:MainTreatmentPS:DroughtTRUE` + `GenotypeBay-0`),
    # Col-0
    Col0_C   = exp(`GenotypeCol-0`),
    Col0_HT  = exp(MainTreatmentHT + `GenotypeCol-0:MainTreatmentHT` + `GenotypeCol-0`),
    Col0_S   = exp(MainTreatmentS  + `GenotypeCol-0:MainTreatmentS`  + `GenotypeCol-0`),
    Col0_PS  = exp(MainTreatmentPS + `GenotypeCol-0:MainTreatmentPS` + `GenotypeCol-0`),
    Col0_D   = exp(DroughtTRUE + `GenotypeCol-0:DroughtTRUE` + `GenotypeCol-0`),
    Col0_HTD = exp(MainTreatmentHT + DroughtTRUE + `MainTreatmentHT:DroughtTRUE` + 
                     `GenotypeCol-0:MainTreatmentHT` + `GenotypeCol-0:DroughtTRUE` +
                     `GenotypeCol-0:MainTreatmentHT:DroughtTRUE` + `GenotypeCol-0`),
    Col0_PSD = exp(MainTreatmentPS + DroughtTRUE + `MainTreatmentPS:DroughtTRUE` +
                     `GenotypeCol-0:MainTreatmentPS` + `GenotypeCol-0:DroughtTRUE` +
                     `GenotypeCol-0:MainTreatmentPS:DroughtTRUE` + `GenotypeCol-0`),
    # Lp2-6
    Lp26_C   = exp(`GenotypeLp2-6`),
    Lp26_HT  = exp(MainTreatmentHT + `GenotypeLp2-6:MainTreatmentHT` + `GenotypeLp2-6`),
    Lp26_S   = exp(MainTreatmentS  + `GenotypeLp2-6:MainTreatmentS`  + `GenotypeLp2-6`),
    Lp26_PS  = exp(MainTreatmentPS + `GenotypeLp2-6:MainTreatmentPS` + `GenotypeLp2-6`),
    Lp26_D   = exp(DroughtTRUE + `GenotypeLp2-6:DroughtTRUE` + `GenotypeLp2-6`),
    Lp26_HTD = exp(MainTreatmentHT + DroughtTRUE + `MainTreatmentHT:DroughtTRUE` + 
                     `GenotypeLp2-6:MainTreatmentHT` + `GenotypeLp2-6:DroughtTRUE` +
                     `GenotypeLp2-6:MainTreatmentHT:DroughtTRUE` + `GenotypeLp2-6`),
    Lp26_PSD = exp(MainTreatmentPS + DroughtTRUE + `MainTreatmentPS:DroughtTRUE` +
                     `GenotypeLp2-6:MainTreatmentPS` + `GenotypeLp2-6:DroughtTRUE` + 
                     `GenotypeLp2-6:MainTreatmentPS:DroughtTRUE` + `GenotypeLp2-6`)) %>%
  dplyr::select(Trait, ID,
                An1_C, An1_HT, An1_S, An1_PS, An1_D, An1_HTD, An1_PSD,
                Bay0_C, Bay0_HT, Bay0_S, Bay0_PS, Bay0_D, Bay0_HTD, Bay0_PSD,
                Col0_C, Col0_HT, Col0_S, Col0_PS, Col0_D, Col0_HTD, Col0_PSD,
                Lp26_C, Lp26_HT, Lp26_S, Lp26_PS, Lp26_D, Lp26_HTD, Lp26_PSD)

# Create average of genotypes
morphology = mutate(morphology, 
                 Average_C   = (An1_C + Bay0_C + Col0_C + Lp26_C)/4,
                 Average_HT  = (An1_HT + Bay0_HT + Col0_HT + Lp26_HT)/4,
                 Average_S   = (An1_S  + Bay0_S  + Col0_S  + Lp26_S)/4,
                 Average_PS  = (An1_PS + Bay0_PS + Col0_PS + Lp26_PS)/4,
                 Average_D   = (An1_D + Bay0_D + Col0_D + Lp26_D)/4,
                 Average_HTD = (An1_HTD + Bay0_HTD + Col0_HTD + Lp26_HTD)/4,
                 Average_PSD = (An1_PSD + Bay0_PSD + Col0_PSD + Lp26_PSD)/4)

# Split dataset accross genotypes, add genotype, rename and rbind it
morphology = map(c("An1", "Bay0", "Col0", "Lp26", "Average"), function(x) {
  out = select(morphology, Trait, ID, contains(x))
  names(out) = c("Trait", "ID", "C","HT", "S", "PS", "D", "HTD", "PSD")
  mutate(out, Genotype = x)}) %>% 
  do.call("rbind", .)

# Reshape to longest format
morphology = pivot_longer(morphology, c(-Trait, -Genotype, -ID), 
                       names_to = "Effect", values_to = "Value")

# Reshape to wider format separating traits
morphology = pivot_wider(morphology, names_from = "Trait", values_from = "Value")


# Save results ------------------------------------------------------------
saveRDS(morphology, file = "Intermediate/Morphology/DataForCorrelations.rds")


# Calculate mean and stdev of each trait ----------------------------------

morphology = readRDS("Intermediate/Morphology/DataForCorrelations.rds")

sum_morphology = group_by(morphology, Genotype, Effect) %>%
    summarise(muDW = median(DW), seDW = sd(DW),
              muRL = median(RootLength), seRL = sd(RootLength),
              muAngle = median(Angle), seAngle = sd(Angle),
              muPR = median(PetioleRatio), sePR = sd(PetioleRatio),
              muBS = median(ShapeBlade), seBS = sd(ShapeBlade),
              muSPA = median(SLA), seSPA = sd(SLA),
              muLS = median(LeafSize), seLS = sd(LeafSize),
              muN = median(N), seN = sd(N),
              muNmol = median(Nmol), seNmol = sd(Nmol))

