
# Compute correlations among traits
# Load packages and data --------------------------------------------------
library(tidyverse)

# rates traits
coefs = readRDS(file = "Intermediate/Rates/LinearFits.rds")

# Reshape the data as required --------------------------------------------
coefs = map(coefs, ~mutate(.x, ID = 1:n()))

# Combine into single table, remove the missing treatment combination 
rates = do.call("rbind", coefs)

rates = rates %>%
  group_by(Trait) %>%
  mutate(# An-1
    An1_C   = Time,
    An1_HT  = Time + `Time:TreatmentHT`,
    An1_S   = Time + `Time:TreatmentS`,
    An1_PS  = Time + `Time:TreatmentPS`,
    An1_D   = Time + `Time:TreatmentD`,
    An1_HTD = Time + `Time:TreatmentHTD`,
    An1_PSD = Time + `Time:TreatmentPSD`,
    # Bay-0
    Bay0_C   = Time + `Time:GenotypeBay-0`,
    Bay0_HT  = Time + `Time:TreatmentHT` + `Time:GenotypeBay-0:TreatmentHT`,
    Bay0_S   = Time + `Time:TreatmentS` + `Time:GenotypeBay-0:TreatmentS`,
    Bay0_PS  = Time + `Time:TreatmentPS` + `Time:GenotypeBay-0:TreatmentPS`,
    Bay0_D   = Time + `Time:TreatmentD` + `Time:GenotypeBay-0:TreatmentD`,
    Bay0_HTD = Time + `Time:TreatmentHTD` + `Time:GenotypeBay-0:TreatmentHTD`,
    Bay0_PSD = Time + `Time:TreatmentPSD` + `Time:GenotypeBay-0:TreatmentPSD`,
    # Col-0
    Col0_C   = Time + `Time:GenotypeCol-0`,
    Col0_HT  = Time + `Time:TreatmentHT` + `Time:GenotypeCol-0:TreatmentHT`,
    Col0_S   = Time + `Time:TreatmentS` + `Time:GenotypeCol-0:TreatmentS`,
    Col0_PS  = Time + `Time:TreatmentPS` + `Time:GenotypeCol-0:TreatmentPS`,
    Col0_D   = Time + `Time:TreatmentD` + `Time:GenotypeCol-0:TreatmentD`,
    Col0_HTD = Time + `Time:TreatmentHTD` + `Time:GenotypeCol-0:TreatmentHTD`,
    Col0_PSD = Time + `Time:TreatmentPSD` + `Time:GenotypeCol-0:TreatmentPSD`,
    # Lp2-6
    Lp26_C  =  Time + `Time:GenotypeLp2-6`,
    Lp26_HT  = Time + `Time:TreatmentHT` + `Time:GenotypeLp2-6:TreatmentHT`,
    Lp26_S   = Time + `Time:TreatmentS` + `Time:GenotypeLp2-6:TreatmentS`,
    Lp26_PS  = Time + `Time:TreatmentPS` + `Time:GenotypeLp2-6:TreatmentPS`,
    Lp26_D   = Time + `Time:TreatmentD` + `Time:GenotypeLp2-6:TreatmentD`,
    Lp26_HTD = Time + `Time:TreatmentHTD` + `Time:GenotypeLp2-6:TreatmentHTD`,
    Lp26_PSD = Time + `Time:TreatmentPSD` + `Time:GenotypeLp2-6:TreatmentPSD`) %>%
  dplyr::select(Trait, ID,
                An1_C, An1_HT, An1_S, An1_PS, An1_D, An1_HTD, An1_PSD,
                Bay0_C, Bay0_HT, Bay0_S, Bay0_PS, Bay0_D, Bay0_HTD, Bay0_PSD,
                Col0_C, Col0_HT, Col0_S, Col0_PS, Col0_D, Col0_HTD, Col0_PSD,
                Lp26_C, Lp26_HT, Lp26_S, Lp26_PS, Lp26_D, Lp26_HTD, Lp26_PSD)

# Create average of genotypes
rates = mutate(rates, 
                    Average_C   = (An1_C + Bay0_C + Col0_C + Lp26_C)/4,
                    Average_HT  = (An1_HT + Bay0_HT + Col0_HT + Lp26_HT)/4,
                    Average_S   = (An1_S  + Bay0_S  + Col0_S  + Lp26_S)/4,
                    Average_PS  = (An1_PS + Bay0_PS + Col0_PS + Lp26_PS)/4,
                    Average_D   = (An1_D + Bay0_D + Col0_D + Lp26_D)/4,
                    Average_HTD = (An1_HTD + Bay0_HTD + Col0_HTD + Lp26_HTD)/4,
                    Average_PSD = (An1_PSD + Bay0_PSD + Col0_PSD + Lp26_PSD)/4)

# Split dataset accross genotypes, add genotype, rename and rbind it
rates = map(c("An1", "Bay0", "Col0", "Lp26", "Average"), function(x) {
  out = select(rates, Trait, ID, contains(x))
  names(out) = c("Trait", "ID", "C","HT", "S", "PS", "D", "HTD", "PSD")
  mutate(out, Genotype = x)}) %>% 
  do.call("rbind", .)

# Reshape to longest format
rates = pivot_longer(rates, c(-Trait, -Genotype, -ID), 
                          names_to = "Effect", values_to = "Value")

# Reshape to wider format separating traits
rates = pivot_wider(rates, names_from = "Trait", values_from = "Value")


# Save results ------------------------------------------------------------
saveRDS(rates, file = "Intermediate/Rates/DataForCorrelations.rds")

# Calculate mean and stdev of each trait ----------------------------------

rates = readRDS("Intermediate/Rates/DataForCorrelations.rds")

sum_rates = group_by(rates, Genotype, Effect) %>%
  summarise(muRGR = median(RGR), seRGR = sd(RGR),
            muLeafRate = median(LeafRate), seLeafRate = sd(LeafRate))

