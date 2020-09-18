# This scripts fits lmer and rlmer models to the different traits and
# 1. Compare the estimates of fixed effects
# 2. Perform diagnostics on the fits

# 01 Load packages and data -----------------------------------------------
library(tidyverse)
library(lme4)
library(robustlmm)
library(robustbase)

# Load the different traits and identify batches of data
morph_traits = readRDS("Intermediate/Morphology/LinearModelTraits.rds")

fitness_traits = readRDS("Intermediate/Fitness/LinearModelTraits.rds")
photosynthesis_traits = readRDS("Intermediate/Photosynthesis/LinearModelTraits.rds")


# Fit linear and robust linear mixed models to each each trait ------------

# Common functions

compare_effects = function(lmm, rlmm) {
  par(mfrow = c(2,2))
  fl = fixef(lmm)
  frl = fixef(rlmm)
  plot(fl, (frl - fl))
  abline(h = 0)

  rl = ranef(lmm)$Batch[[1]]
  rrl = ranef(rlmm)$Batch[[1]]
  plot(rl, (rrl - rl))
  abline(h = 0)
  
  plot(residuals(lmm), residuals(rlmm))
  abline(a = 0, b = 1)
  
  hist(residuals(lmm))
  print(shapiro.test(residuals(lmm)))
  print(shapiro.test(rl))
}

# DW
dw_data = filter(morph_traits, !is.na(DW))

lmm_DW = lmer(DW~Genotype*MainTreatment*Drought + Genotype*Timepoint - 1 + (1 | Batch),
              data = dw_data)
rlmm_DW = rlmerRcpp(DW~Genotype*MainTreatment*Drought + Genotype*Timepoint - 1 + (1 | Batch),
              data = dw_data)

compare_effects(lmm_DW, rlmm_DW)

#plot(rlmm_DW)

# Yield
y_data = filter(fitness_traits, !is.na(Yield), Yield > 0, Yield < 1000)

lmm_Yield = lmer(Yield~Genotype*MainTreatment*Drought - 1 + (1 | Batch),
              data = y_data)
rlmm_Yield = rlmerRcpp(Yield~Genotype*MainTreatment*Drought - 1 + (1 | Batch),
                    data = y_data)

compare_effects(lmm_Yield, rlmm_Yield)

#plot(rlmm_Yield)

# SLA
sla_data = filter(morph_traits, !is.na(SLA))

lmm_SLA = lmer(SLA~Genotype*MainTreatment*Drought + Genotype*Timepoint - 1 + (1 | Batch),
                 data = sla_data)
rlmm_SLA = rlmerRcpp(SLA~Genotype*MainTreatment*Drought + Genotype*Timepoint - 1 + (1 | Batch),
                       data = sla_data)

compare_effects(lmm_SLA, rlmm_SLA)


# Rd
Rd_data = filter(photosynthesis_traits, !is.na(Rd))

lmm_Rd = lm(Rd~MainTreatment*Drought, data = Rd_data)
rlmm_Rd = lmrob(Rd~MainTreatment*Drought, data = Rd_data)

plot(coef(lmm_Rd), coef(rlmm_Rd) - coef(lmm_Rd))
abline(h = 0)

hist(residuals(lmm_Rd))
shapiro.test(residuals(lmm_Rd))

# Conclusion: We use robust linear mixed model for morphology and fitness to account
# for the batch effect. For photosynthesis, the mixed model is not feasible as 
# batches contain too few data and results in singular fits (and these batches
# were the result of staggered sowing, not independent repetitions). 
# Note that the formulae are different depending on the type of experiment 
# (and not all MainTreatment are present in all 
# the datasets). 
# Note that in the photosynthetic parameters we may not need robust fits after all.