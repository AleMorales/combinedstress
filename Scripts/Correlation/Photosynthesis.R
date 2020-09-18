
# Compute correlations among traits
# Load packages and data --------------------------------------------------
library(tidyverse)

# photosynthesis traits
coefs = readRDS(file = "Intermediate/Photosynthesis/LinearFits.rds")

# Reshape the data as required --------------------------------------------
coefs = map(coefs, ~mutate(.x, ID = 1:n()))

# Combine into single table, remove the missing treatment combination 
photosynthesis = do.call("rbind", coefs)

photosynthesis = photosynthesis %>%
  group_by(Trait) %>%
  mutate(
    C  = exp(`(Intercept)`),
    HT  = exp(MainTreatmentHT + `(Intercept)`),
    S   = exp(MainTreatmentS + `(Intercept)`),
    PS  = exp(MainTreatmentPS + `(Intercept)`),
    D   = exp(DroughtTRUE + `(Intercept)`),
    HTD = exp(MainTreatmentHT + DroughtTRUE + `MainTreatmentHT:DroughtTRUE` + `(Intercept)`),
    PSD = exp(MainTreatmentPS + DroughtTRUE + `MainTreatmentPS:DroughtTRUE` + `(Intercept)`)) %>%
  dplyr::select(Trait, ID, C, HT, S, PS, D, HTD, PSD) %>%
  mutate(Genotype = "Average")


# Reshape to longest format
photosynthesis = pivot_longer(photosynthesis, c(-Trait, -Genotype, -ID), 
                          names_to = "Effect", values_to = "Value")

# Reshape to wider format separating traits
photosynthesis = pivot_wider(photosynthesis, names_from = "Trait", values_from = "Value")


# Save results ------------------------------------------------------------
saveRDS(photosynthesis, file = "Intermediate/Photosynthesis/DataForCorrelations.rds")

# Calculate mean and stdev of each trait ----------------------------------

sum_photosynthesis = group_by(photosynthesis, Effect) %>%
  summarise(muRd = median(Rd), seRd = sd(Rd),
            muLUE = median(LUE)*1e2, seLUE = sd(LUE)*1e2,
            muA = median(A), seA = sd(A),
            muChl = median(Chl), seChl = sd(Chl),
            mugsw = median(gsw)*1e3, segsw = sd(gsw)*1e3)

