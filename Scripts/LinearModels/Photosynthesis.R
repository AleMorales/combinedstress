
# Load packages and data ------------------------------------------------------
library(tidyverse)
library(furrr)
library(robustbase)

traitdata = readRDS("Intermediate/Photosynthesis/LinearModelTraits.rds")

plan(multiprocess)

# Robust linear models on photosynthetic traits with residual bootstrap -------

logit = function(x) log(x/(1 - x))
logistic = function(x) 1/(1 + exp(-x))

# Fit robust linear model to data and return information for bootstrapping
intial_fit = function(trait, data, transformation = log) {
  data = filter(data, !is.na(!!trait)) %>% 
    dplyr::select(MainTreatment, Drought, !!trait) %>%
    rename(trait = !!trait) %>% 
    mutate(trait = transformation(trait))
  
  fit = lmrob(trait~MainTreatment*Drought, data = data,
              control = lmrob.control(k.max = 1000, setting="KS2014",
                      max.it = 1000, maxit.scale = 1e3, refine.tol = 1e-5))
  
  rcoef = coef(fit)
  
  res = residuals(fit)

  return(list(data = data, fit = fit, res = res, rcoef = rcoef))
}

# Create the function for refitting the model to a bootstrap sample
create_refit_model = function(trait_list) {
  data = trait_list$data
  fit = trait_list$fit
  res = trait_list$res
  refit_model  = function(x) {
    new_data = mutate(data, trait = predict(fit) + sample(res, replace = TRUE))
    new_fit = update(fit, data = new_data)
    coef(new_fit)
  }
}

N = 5e3

# Rd
Rd_list = intial_fit(quo(Rd), traitdata)
refit_Rd = create_refit_model(Rd_list)
Rd_coef = future_map_dfr(1:N, refit_Rd, .progress = T, 
                         .options = future_options(packages = "robustbase"))

# LUE
LUE_list = intial_fit(quo(LUE), traitdata)
refit_LUE = create_refit_model(LUE_list)
LUE_coef = future_map_dfr(1:N, refit_LUE, .progress = T, 
                         .options = future_options(packages = "robustbase"))

# A
A_list = intial_fit(quo(A), traitdata)
refit_A = create_refit_model(A_list)
A_coef = future_map_dfr(1:N, refit_A, .progress = T, 
                          .options = future_options(packages = "robustbase"))

# Chl
Chl_list = intial_fit(quo(Chl), traitdata)
refit_Chl = create_refit_model(Chl_list)
Chl_coef = future_map_dfr(1:N, refit_Chl, .progress = T, 
                          .options = future_options(packages = "robustbase"))

# gsw
gsw_list = intial_fit(quo(gsw), traitdata)
refit_gsw = create_refit_model(gsw_list)
gsw_coef = future_map_dfr(1:N, refit_gsw, .progress = T, 
                          .options = future_options(packages = "robustbase"))


# Combine estimations into a single table ---------------------------------

# Add trait name
Rd_coef = mutate(Rd_coef, Trait = 'Rd')
LUE_coef = mutate(LUE_coef, Trait = 'LUE')
A_coef = mutate(A_coef, Trait = 'A')
Chl_coef = mutate(Chl_coef, Trait = 'Chl')
gsw_coef = mutate(gsw_coef, Trait = 'gsw')

saveRDS(list(Rd_coef, LUE_coef, A_coef, Chl_coef, gsw_coef), 
        file = "Intermediate/Photosynthesis/LinearFits.rds")

