# Load packages and data ------------------------------------------------------
library(tidyverse)

seedlength = readRDS("Intermediate/Biosorter/LinearModel.rds")

seedlength = mutate(seedlength, ID = 1:n())

# Calculate predicted trait values per genotype & treatment -------------------
seedlength = seedlength %>%
  mutate(# An-1
    An1_C   = `GenotypeAn-1`,
    An1_HT  = `GenotypeAn-1` + MainTreatmentHT,
    An1_PS  = `GenotypeAn-1` + MainTreatmentPS,
    An1_D   = `GenotypeAn-1` + DroughtTRUE,
    An1_HTD = `GenotypeAn-1` + MainTreatmentHT + DroughtTRUE + `MainTreatmentHT:DroughtTRUE`,
    An1_PSD = `GenotypeAn-1` + MainTreatmentPS + DroughtTRUE + `MainTreatmentPS:DroughtTRUE`,
    # Bay-0
    Bay0_C   = `GenotypeBay-0`,
    Bay0_HT  = `GenotypeBay-0` + MainTreatmentHT + `GenotypeBay-0:MainTreatmentHT`,
    Bay0_PS  = `GenotypeBay-0` + MainTreatmentPS + `GenotypeBay-0:MainTreatmentPS`,
    Bay0_D   = `GenotypeBay-0` + DroughtTRUE + `GenotypeBay-0:DroughtTRUE`,
    Bay0_HTD = `GenotypeBay-0` + MainTreatmentHT + `GenotypeBay-0:MainTreatmentHT` +
                     DroughtTRUE + `GenotypeBay-0:DroughtTRUE` +
                     `MainTreatmentHT:DroughtTRUE` + `GenotypeBay-0:MainTreatmentHT:DroughtTRUE`,
    Bay0_PSD = `GenotypeBay-0` + MainTreatmentPS + `GenotypeBay-0:MainTreatmentPS` +
                     DroughtTRUE + `GenotypeBay-0:DroughtTRUE` +
                     `MainTreatmentPS:DroughtTRUE` + `GenotypeBay-0:MainTreatmentPS:DroughtTRUE`,
    # Col-0
    Col0_C   = `GenotypeCol-0`,
    Col0_HT  = `GenotypeCol-0` + MainTreatmentHT + `GenotypeCol-0:MainTreatmentHT`,
    Col0_PS  = `GenotypeCol-0` + MainTreatmentPS + `GenotypeCol-0:MainTreatmentPS`,
    Col0_D   = `GenotypeCol-0` + DroughtTRUE + `GenotypeCol-0:DroughtTRUE`,
    Col0_HTD = `GenotypeCol-0` + MainTreatmentHT + `GenotypeCol-0:MainTreatmentHT` +
                     DroughtTRUE + `GenotypeCol-0:DroughtTRUE` +
                     `MainTreatmentHT:DroughtTRUE` + `GenotypeCol-0:MainTreatmentHT:DroughtTRUE`,
    Col0_PSD = `GenotypeCol-0` + MainTreatmentPS + `GenotypeCol-0:MainTreatmentPS` +
                     DroughtTRUE + `GenotypeCol-0:DroughtTRUE` +
                     `MainTreatmentPS:DroughtTRUE` + `GenotypeCol-0:MainTreatmentPS:DroughtTRUE`) %>%
  dplyr::select(ID,
                An1_C,  An1_HT,  An1_PS, An1_D, An1_HTD, An1_PSD,
                Bay0_C, Bay0_HT, Bay0_PS, Bay0_D, Bay0_HTD, Bay0_PSD,
                Col0_C, Col0_HT, Col0_PS, Col0_D, Col0_HTD, Col0_PSD)

# Create average of genotypes
seedlength = mutate(seedlength, 
                 Average_C   = (An1_C   + Bay0_C   + Col0_C)/3,
                 Average_HT  = (An1_HT  + Bay0_HT  + Col0_HT)/3,
                 Average_PS  = (An1_PS  + Bay0_PS  + Col0_PS)/3,
                 Average_D   = (An1_D   + Bay0_D   + Col0_D)/3,
                 Average_HTD = (An1_HTD + Bay0_HTD + Col0_HTD)/3,
                 Average_PSD = (An1_PSD + Bay0_PSD + Col0_PSD)/3)

# Split dataset accross genotypes, add genotype, rename and rbind it
seedlength = map(c("An1", "Bay0", "Col0", "Average"), function(x) {
  out = select(seedlength, ID, contains(x))
  names(out) = c('ID', "C", "HT", "PS", "D", "HTD", "PSD")
  mutate(out, Genotype = x)}) %>% 
  do.call("rbind", .)

# Reshape to long format
seedlength = pivot_longer(seedlength, c(-Genotype, -ID), 
                       names_to = "Treatment", values_to = "Value")

# Summary statistics
sum_seedlength = group_by(seedlength, Genotype, Treatment) %>% 
  summarise(mu = median(Value))

saveRDS(sum_seedlength, file = "Intermediate/Biosorter/DataForPCA.rds")
saveRDS(seedlength, file = "Intermediate/Biosorter/DataForPCA_all.rds")
