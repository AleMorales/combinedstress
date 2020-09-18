# Load packages and data --------------------------------------------------
library(tidyverse)

coefs = readRDS(file = "Intermediate/Rates/LinearFits.rds")

# Combine into single table, remove the missing treatment combination 
coefs = map(coefs, ~mutate(.x, ID = 1:n()))

rates = do.call("rbind", coefs)

# Compute the different combinations to get the effects on RGR and LeafRate

rates = rates %>% 
  group_by(Trait) %>%
  mutate(# An-1
    An1_C   = Time,
    An1_HT  = `Time:TreatmentHT` + Time,
    An1_S   = `Time:TreatmentS` + Time,
    An1_PS  = `Time:TreatmentPS` + Time,
    An1_D   = `Time:TreatmentD` + Time,
    An1_HTD = `Time:TreatmentHTD` + Time,
    An1_PSD = `Time:TreatmentPSD` + Time,
    # Bay-0
    Bay0_C   = Time + `Time:GenotypeBay-0`,
    Bay0_HT  = (`Time:TreatmentHT` + `Time:GenotypeBay-0:TreatmentHT`) + (Time + `Time:GenotypeBay-0`),
    Bay0_S    = (`Time:TreatmentS` + `Time:GenotypeBay-0:TreatmentS`) + (Time + `Time:GenotypeBay-0`),
    Bay0_PS  = (`Time:TreatmentPS` + `Time:GenotypeBay-0:TreatmentPS`) + (Time + `Time:GenotypeBay-0`),
    Bay0_D   = (`Time:TreatmentD`  + `Time:GenotypeBay-0:TreatmentD`) + (Time + `Time:GenotypeBay-0`),
    Bay0_HTD = (`Time:TreatmentHTD` + `Time:GenotypeBay-0:TreatmentHTD`) + (Time + `Time:GenotypeBay-0`),
    Bay0_PSD = (`Time:TreatmentPSD` + `Time:GenotypeBay-0:TreatmentPSD`) + (Time + `Time:GenotypeBay-0`),
    # Col-0
    Col0_C   = Time + `Time:GenotypeCol-0`,
    Col0_HT  = (`Time:TreatmentHT` + `Time:GenotypeCol-0:TreatmentHT`) + (Time + `Time:GenotypeCol-0`),
    Col0_S    = (`Time:TreatmentS` + `Time:GenotypeCol-0:TreatmentS`) + (Time + `Time:GenotypeCol-0`),
    Col0_PS  = (`Time:TreatmentPS` + `Time:GenotypeCol-0:TreatmentPS`) + (Time + `Time:GenotypeCol-0`),
    Col0_D   = (`Time:TreatmentD`  + `Time:GenotypeCol-0:TreatmentD`) + (Time + `Time:GenotypeCol-0`),
    Col0_HTD = (`Time:TreatmentHTD` + `Time:GenotypeCol-0:TreatmentHTD`) + (Time + `Time:GenotypeCol-0`),
    Col0_PSD = (`Time:TreatmentPSD` + `Time:GenotypeCol-0:TreatmentPSD`) + (Time + `Time:GenotypeCol-0`),
    # Lp2-6
    Lp26_C   = Time + `Time:GenotypeLp2-6`,
    Lp26_HT  = (`Time:TreatmentHT` + `Time:GenotypeLp2-6:TreatmentHT`) + (Time + `Time:GenotypeLp2-6`),
    Lp26_S    = (`Time:TreatmentS` + `Time:GenotypeLp2-6:TreatmentS`) + (Time + `Time:GenotypeLp2-6`),
    Lp26_PS  = (`Time:TreatmentPS` + `Time:GenotypeLp2-6:TreatmentPS`) + (Time + `Time:GenotypeLp2-6`),
    Lp26_D   = (`Time:TreatmentD`  + `Time:GenotypeLp2-6:TreatmentD`) + (Time + `Time:GenotypeLp2-6`),
    Lp26_HTD = (`Time:TreatmentHTD` + `Time:GenotypeLp2-6:TreatmentHTD`) + (Time + `Time:GenotypeLp2-6`),
    Lp26_PSD = (`Time:TreatmentPSD` + `Time:GenotypeLp2-6:TreatmentPSD`) + (Time + `Time:GenotypeLp2-6`)) %>%
  dplyr::select(Trait, ID, 
                An1_C, An1_HT,  An1_S, An1_PS, An1_D, An1_HTD, An1_PSD,
                Bay0_C, Bay0_HT, Bay0_S, Bay0_PS, Bay0_D, Bay0_HTD, Bay0_PSD,
                Col0_C, Col0_HT, Col0_S, Col0_PS, Col0_D, Col0_HTD, Col0_PSD,
                Lp26_C, Lp26_HT, Lp26_S, Lp26_PS, Lp26_D, Lp26_HTD, Lp26_PSD)


# Create average of genotypes
rates = mutate(rates, 
               Average_C   = (An1_C   + Bay0_C   + Col0_C   + Lp26_C)/4,
               Average_HT  = (An1_HT  + Bay0_HT  + Col0_HT  + Lp26_HT)/4,
               Average_S   = (An1_S   + Bay0_S   + Col0_S   + Lp26_S)/4,
               Average_PS  = (An1_PS  + Bay0_PS  + Col0_PS  + Lp26_PS)/4,
               Average_D   = (An1_D   + Bay0_D   + Col0_D   +  Lp26_D)/4,
               Average_HTD = (An1_HTD + Bay0_HTD + Col0_HTD + Lp26_HTD)/4,
               Average_PSD = (An1_PSD + Bay0_PSD + Col0_PSD + Lp26_PSD)/4)


# Split dataset accross genotypes, add genotype, rename and rbind it
rates = map(c("An1", "Bay0", "Col0", "Lp26", "Average"), function(x) {
  out = select(rates, Trait, ID, contains(x))
  names(out) = c("Trait", "ID", "C", "HT", "S", "PS", "D", "HTD", "PSD")
  mutate(out, Genotype = x)}) %>% 
  do.call("rbind", .)

# Reshape to long format
rates = pivot_longer(rates, c(-Trait, -Genotype, -ID), 
                          names_to = "Treatment", values_to = "Value")

# Summary statistics
sum_rates = group_by(rates, Genotype, Trait, Treatment) %>% 
  summarise(mu = median(Value))


saveRDS(sum_rates, file = "Intermediate/Rates/DataForPCA.rds")

saveRDS(rates, file = "Intermediate/Rates/DataForPCA_all.rds")
