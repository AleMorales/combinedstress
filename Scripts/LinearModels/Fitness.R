
# 01 Load packages and data ---------------------------------------------------
library(tidyverse)
library(furrr)
library(robustlmm)
library(tidytext)

traitdata = readRDS("Intermediate/Fitness/LinearModelTraits.rds")

plan(multiprocess, workers = 8)

# Robust linear mixed models on fitness traits with residual/ranef bootstrap --

# Fit robust linear mixed model to the data and return information for bootstrapping
intial_fit = function(trait, data, transformation = log, lwr = 0, upr = Inf) {
  
  data = filter(data, !is.na(!!trait), !!trait > lwr, !!trait < upr) %>% 
    dplyr::select(Genotype, MainTreatment, Drought, Batch, !!trait) %>%
    rename(trait = !!trait) %>% 
    mutate(trait = transformation(trait)) %>%
    arrange(Batch)
  
  fit = rlmerRcpp(trait~Genotype*MainTreatment*Drought - 1 + (1 | Batch), data = data)

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

# FT
FT_fit = intial_fit(quo(FT), traitdata)
FT_fun = create_refit_model(FT_fit)
FT_coef = future_map_dfr(1:N, ~FT_fun(), .progress = T, 
                         .options = future_options(packages = "robustlmm"))

# Yield
Yield_fit = intial_fit(quo(Yield), traitdata, lwr = 0, upr = 1000)
Yield_fun = create_refit_model(Yield_fit)
Yield_coef = future_map_dfr(1:N, ~Yield_fun(), .progress = T, 
                         .options = future_options(packages = "robustlmm"))

# Phyllochron
Phyllochron_fit = intial_fit(quo(Phyllochron), traitdata)
Phyllochron_fun = create_refit_model(Phyllochron_fit)
Phyllochron_coef = future_map_dfr(1:N, ~Phyllochron_fun(), .progress = T, 
                            .options = future_options(packages = "robustlmm"))

# LN
LN_fit = intial_fit(quo(LN), traitdata)
LN_fun = create_refit_model(LN_fit)
LN_coef = future_map_dfr(1:N, ~LN_fun(), .progress = T, 
                                  .options = future_options(packages = "robustlmm"))


# Combine estimations into a single table ---------------------------------

# Add trait name
FT_coef = mutate(FT_coef, Trait = 'FT')
Yield_coef = mutate(Yield_coef, Trait = 'Yield')
Phyllochron_coef = mutate(Phyllochron_coef, Trait = 'Phyllochron')
LN_coef = mutate(LN_coef, Trait = 'LN')

saveRDS(list(FT_coef, Yield_coef, Phyllochron_coef), 
        file = "Intermediate/Fitness/LinearFits.rds")



coefs = readRDS("Intermediate/Fitness/LinearFits.rds")
coefs = c(coefs, LN_coef)
saveRDS(coefs, 
        file = "Intermediate/Fitness/LinearFits.rds")
