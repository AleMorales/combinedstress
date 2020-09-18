# Load packages and data ------------------------------------------------------
library(tidyverse)

coefs = readRDS(file = "Intermediate/Photosynthesis/LinearFits.rds")

# Calculate relative effects on fitness traits ----------------------------

coefs = map(coefs, ~mutate(.x, ID = 1:n()))

# Combine into single table, remove the missing treatment combination 
photosynthesis = do.call("rbind", coefs)
photosynthesis = photosynthesis[,-8]

# Calculate morphology trait values per genotype & treatment -------------------
photosynthesis = photosynthesis %>%
  group_by(Trait) %>%
  mutate(
    C   = exp(`(Intercept)`),
    HT  = exp(`(Intercept)` + MainTreatmentHT),
    S   = exp(`(Intercept)` + MainTreatmentS),
    PS  = exp(`(Intercept)` + MainTreatmentPS),
    D   = exp(`(Intercept)` + DroughtTRUE),
    HTD = exp(`(Intercept)` + MainTreatmentHT + DroughtTRUE + `MainTreatmentHT:DroughtTRUE`),
    PSD = exp(`(Intercept)` + MainTreatmentPS + DroughtTRUE + `MainTreatmentPS:DroughtTRUE`))  %>%
  dplyr::select(Trait, ID, C, HT, S, PS, D, HTD, PSD)

# Reshape to long format
photosynthesis = pivot_longer(photosynthesis, c(-Trait, -ID), 
                              names_to = "Treatment", values_to = "Value")

# Summary statistics
sum_photosynthesis = group_by(photosynthesis, Trait, Treatment) %>%
  summarise(mu = median(Value))


saveRDS(sum_photosynthesis, file = "Intermediate/Photosynthesis/DataForPCA.rds")
saveRDS(photosynthesis, file = "Intermediate/Photosynthesis/DataForPCA_all.rds")
