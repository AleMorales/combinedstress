
# 01 Load packages and data ---------------------------------------------------
library(tidyverse)
library(furrr)
library(robustlmm)
library(tidytext)

traitdata = readRDS("Intermediate/Morphology/LinearModelTraits.rds")

plan(multiprocess, workers = 1)

# Robust linear mixed models on fitness traits with residual/ranef bootstrap --

logit = function(x) log(x/(1 - x))
logistic = function(x) 1/(1 + exp(-x))

# Fit robust linear mixed model to the data and return information for bootstrapping
intial_fit = function(trait, data, transformation = log, lwr = 0, upr = Inf) {
  
  data = filter(data, !is.na(!!trait), !!trait > lwr, !!trait < upr) %>% 
    dplyr::select(Genotype, MainTreatment, Drought, Timepoint, Batch, !!trait) %>%
    rename(trait = !!trait) %>% 
    mutate(trait = transformation(trait)) %>%
    arrange(Batch)
  
  fit = rlmerRcpp(trait~Genotype*MainTreatment*Drought + Genotype*Timepoint - 
                  1 + (1 | Batch), data = data)
  
  return(fit)
}

# Create the function for refitting the model to a bootstrap sample
create_refit_model = function(fit) {
  data = fit@frame
  res = residuals(fit)
  ranef = ranef(fit)$Batch[[1]]
  blocks = as.numeric(table(data$Batch))
  refit_model  = function() {
    new_data = mutate(data, trait = predict(fit, re.form = NA) + 
                        sample(res, replace = TRUE) + 
                        rep(sample(ranef, replace = TRUE), blocks))
    new_fit = update(fit, data = new_data)
    fixef(new_fit)
  }
}


N = 5e3

# DW
DW_fit = intial_fit(quo(DW), traitdata)
DW_fun = create_refit_model(DW_fit)
DW_coef = future_map_dfr(1:N, ~DW_fun(), .progress = T, 
                         .options = future_options(packages = "robustlmm"))

# RootLength
RootLength_fit = intial_fit(quo(RootLength), traitdata)
RootLength_fun = create_refit_model(RootLength_fit)
RootLength_coef = future_map_dfr(1:N, ~RootLength_fun(), .progress = T, 
                            .options = future_options(packages = "robustlmm"))

# Angle
Angle_fit = intial_fit(quo(Angle), traitdata)
Angle_fun = create_refit_model(Angle_fit)
Angle_coef = future_map_dfr(1:N, ~Angle_fun(), .progress = T, 
                                  .options = future_options(packages = "robustlmm"))

# RWC
RWC_fit = intial_fit(quo(RWC), traitdata, logit)
RWC_fun = create_refit_model(RWC_fit)
RWC_coef = future_map_dfr(1:N, ~RWC_fun(), .progress = T, 
                            .options = future_options(packages = "robustlmm"))

# PetioleRatio
PetioleRatio_fit = intial_fit(quo(PetioleRatio), traitdata, logit, upr = 0.9)
PetioleRatio_fun = create_refit_model(PetioleRatio_fit)
PetioleRatio_coef = future_map_dfr(1:N, ~PetioleRatio_fun(), .progress = T, 
                          .options = future_options(packages = "robustlmm"))

# ShapeBlade
ShapeBlade_fit = intial_fit(quo(ShapeBlade), traitdata)
ShapeBlade_fun = create_refit_model(ShapeBlade_fit)
ShapeBlade_coef = future_map_dfr(1:N, ~ShapeBlade_fun(), .progress = T, 
                                   .options = future_options(packages = "robustlmm"))

# SLA
SLA_fit = intial_fit(quo(SLA), traitdata)
SLA_fun = create_refit_model(SLA_fit)
SLA_coef = future_map_dfr(1:N, ~SLA_fun(), .progress = T, 
                                 .options = future_options(packages = "robustlmm"))

# LeafSize
LeafSize_fit = intial_fit(quo(LeafSize), traitdata)
LeafSize_fun = create_refit_model(LeafSize_fit)
LeafSize_coef = future_map_dfr(1:N, ~LeafSize_fun(), .progress = T, 
                          .options = future_options(packages = "robustlmm"))

# N
N_fit = intial_fit(quo(N), traitdata)
N_fun = create_refit_model(N_fit)
N_coef = future_map_dfr(1:N, ~N_fun(), .progress = T, 
                               .options = future_options(packages = "robustlmm"))

# Nmol
Nmol_fit = intial_fit(quo(Nmol), traitdata)
Nmol_fun = create_refit_model(Nmol_fit)
Nmol_coef = future_map_dfr(1:N, ~Nmol_fun(), .progress = T, 
                               .options = future_options(packages = "robustlmm"))

# Combine estimations into a single table ---------------------------------

# Add trait name
DW_coef = mutate(DW_coef, Trait = 'DW')
RootLength_coef = mutate(RootLength_coef, Trait = 'RootLength')
Angle_coef = mutate(Angle_coef, Trait = 'Angle')
RWC_coef = mutate(RWC_coef, Trait = 'RWC')
PetioleRatio_coef = mutate(PetioleRatio_coef, Trait = 'PetioleRatio')
ShapeBlade_coef = mutate(ShapeBlade_coef, Trait = 'ShapeBlade')
SLA_coef = mutate(SLA_coef, Trait = 'SLA')
LeafSize_coef = mutate(LeafSize_coef, Trait = 'LeafSize')
N_coef = mutate(N_coef, Trait = 'N')
Nmol_coef = mutate(Nmol_coef, Trait = 'Nmol')


saveRDS(list(DW_coef, RootLength_coef, Angle_coef, RWC_coef, PetioleRatio_coef,
             ShapeBlade_coef, SLA_coef, LeafSize_coef, N_coef, Nmol_coef), 
        file = "Intermediate/Morphology/LinearFits.rds")
