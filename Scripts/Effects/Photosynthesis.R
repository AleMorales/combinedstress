
# Load packages and data ------------------------------------------------------
library(tidyverse)

coefs = readRDS(file = "Intermediate/Photosynthesis/LinearFits.rds")

# Calculate relative effects on fitness traits ----------------------------

# Combine into single table, remove the missing treatment combination 
photosynthesis = do.call("rbind", coefs)
photosynthesis = photosynthesis[,-8]


# Calculate the relative effects of each treatment
relative_main_effect = function(effect, control, back = exp) {
  (back(effect + control) - back(control))/(back(control))
}
relative_interaction = function(effect1, effect2, interaction, control, back = exp) {
  (back(effect1 + effect2 + interaction + control) - 
     back(effect1 + effect2 + control))/(back(control))
}
photosynthesis = photosynthesis %>%
  group_by(Trait) %>%
  mutate(HT   = relative_main_effect(MainTreatmentHT, `(Intercept)`),
         S    = relative_main_effect(MainTreatmentS, `(Intercept)`),
         PS   = relative_main_effect(MainTreatmentPS, `(Intercept)`),
         D    = relative_main_effect(DroughtTRUE, `(Intercept)`),
         HTD  = relative_main_effect(MainTreatmentHT + DroughtTRUE + `MainTreatmentHT:DroughtTRUE`, `(Intercept)`),
         PSD  = relative_main_effect(MainTreatmentPS + DroughtTRUE + `MainTreatmentPS:DroughtTRUE`, `(Intercept)`),
         HT_D = relative_interaction(MainTreatmentHT,  DroughtTRUE, `MainTreatmentHT:DroughtTRUE`, `(Intercept)`),
         PS_D = relative_interaction(MainTreatmentPS,  DroughtTRUE, `MainTreatmentPS:DroughtTRUE`, `(Intercept)`)) %>%
  dplyr::select(Trait, HT, S, PS, D, HTD, PSD, HT_D, PS_D)

# Reshape to long format
photosynthesis = pivot_longer(photosynthesis, -Trait, 
                              names_to = "Effect", values_to = "Value")

# Summary statistics
sum_photosynthesis = group_by(photosynthesis, Trait, Effect) %>%
  summarise(mu = median(Value))

# Add the mean value to each trait-effect for reorder
photosynthesis = merge(photosynthesis, sum_photosynthesis)

photosynthesis$Effect = factor(photosynthesis$Effect, 
                               levels = c("HT", "HTD", "HT_D", "PS",  "PSD", "PS_D", "S", "D"))

saveRDS(photosynthesis, file = "Intermediate/Photosynthesis/Effects.rds")

